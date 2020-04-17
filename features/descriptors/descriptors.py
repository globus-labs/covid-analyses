import argparse
import os
import parsl
import pickle
import time
from pathlib import Path
from parsl.app.app import python_app

@python_app
def compute_descriptors(smiles=None, smiles_file=None, start_index=0, batch_size=0, out_file=None, bad_file=None, save_csv=False, ignore_3D=False):
    from mordred import Calculator, descriptors
    from rdkit import Chem
    import numpy as np
    import pickle
    import csv

    if smiles_file: 
        with open(smiles_file) as current:
            current.seek(start_index)
            smiles = [current.readline() for i in range(batch_size)]

    bad = []
    csv_results = []
    results = {}
    mol = None

    calc = Calculator(descriptors, ignore_3D=ignore_3D)

    for smile_tuple in smiles:
        s = smile_tuple.split(',')
        smile = s[2].rstrip()
        identifier = s[1].rstrip()
        
        try:
            mol = Chem.MolFromSmiles(smile)

            if mol is None or len(mol.GetAtoms()) > 100:
                bad.append((smile, identifier))
            else:
                descs = calc(mol)

                # Mods from Xuefeng
                descs.fill_missing("nan")
                
                descriptor_array = np.array(descs).flatten().astype(np.float32)
                if save_csv: 
                    data = [smile, identifier] + list(descriptor_array)
                    data = ['' if str(i) == 'nan' else i for i in data]
                    csv_results.append(data) 
                else: 
                    results[smile] = ([identifier], np.array(descs).flatten().astype(np.float32))
        except:
            bad.append((smile, identifier))

    if save_csv: 
        with open(out_file, 'w', newline='') as o_file:
            writer = csv.writer(o_file, delimiter=',', quoting=csv.QUOTE_MINIMAL)
            writer.writerows(csv_results)
        if bad_file and len(bad) > 0:
            with open(bad_file, 'w', newline='') as b_file:
                b_writer = csv.writer(b_file, delimiter=',', quoting=csv.QUOTE_MINIMAL)
                for b in bad:
                    b_writer.writerow(b)
    else:
        with open(out_file, 'wb') as o_file:
            pickle.dump(results, o_file, protocol=pickle.HIGHEST_PROTOCOL)

        if bad_file and len(bad) > 0:
            with open(bad_file, 'wb') as b_file:
                pickle.dump(bad, b_file, protocol=pickle.HIGHEST_PROTOCOL)
         
    return out_file

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-v", "--version",
                        help="Print Endpoint version information")
    parser.add_argument("-d", "--debug", action='store_true',
                        help="Enables debug logging")
    parser.add_argument("-n", "--num_smiles", default=0, type=int,
                        help="Number of smiles to load and run. Default=10000, if set to 0 the entire file will be used")
    parser.add_argument("-i", "--input_file", required=True,
                        help="File path to the input smiles csv file")
    parser.add_argument("-bo", "--bad_output_dir", default=None,
                        help="Output directory for bad smile list")
    parser.add_argument("-b", "--batch_size", default=0, type=int,
                        help="Size of the batch of smiles to send to each node for processing. Batch size of 0 will process all smiles in a single batch. Default=0.")
    parser.add_argument("-o", "--output_dir", default="outputs",
                        help="Output directory. Default : outputs")
    parser.add_argument("-off", "--offset", default=0, type=int,
                        help="Offset for starting processing from separate files. Default=0")
    parser.add_argument("-csv", "--csv", default=False, action='store_true',
                        help="Save output as CSV. Default is to save as Pickle.")
    parser.add_argument("-ig", "--ignore3d", default=False, action='store_true',
                        help="Ignore 3D descriptors. True: 1613 descriptors; False: 1826 descriptors. Default False.")
    parser.add_argument("-c", "--config", default="local",
                        help="Parsl config defining the target compute resource to use. Default: local")
    args = parser.parse_args()

    if args.config == "local":
        from parsl.configs.htex_local import config
        from parsl.configs.htex_local import config
        config.executors[0].label = "Foo"
        config.executors[0].max_workers = 5
    elif args.config == "theta":
        from theta import config
    elif args.config == "comet":
        from comet import config
    elif args.config == "cvd":
        from configs.theta_cvd import config

    parsl.load(config)

    if args.debug:
        parsl.set_stream_logger()

    output_dir = args.output_dir
    input_file = args.input_file
    bad_dir = args.bad_output_dir
    bad_file = None
    chunksize = int(args.batch_size)
    ignore_3D = args.ignore3d
    save_csv = args.csv
    extension = "csv" if save_csv else 'pkl'

    with open(input_file) as f:
        if int(args.num_smiles) > 0:
            smiles = f.readlines()[:int(args.num_smiles)]
        else:
            smiles = f.readlines()

    try:
        os.makedirs(output_dir)
    except:
        pass
    try:
        os.makedirs(bad_dir)
    except:
        pass

    num_smiles = len(smiles)
    print(f"[Main] Processing {input_file}")
    print(f"[Main] Ignore3D: {ignore_3D}, SaveCSV: {save_csv}")
    print(f"[Main] Chunksize: {chunksize}, Smiles: {num_smiles}")
    print(f"[Main] Output: {output_dir}")
    batch_futures = {}
   
    f = Path(input_file).stem 

    start = time.time()

    if chunksize == 0 or chunksize > len(smiles):
        output_file = os.path.join(output_dir, '%s.%s' % (f, extension))
        if bad_dir:
            bad_file = os.path.join(bad_dir, "%s.%s" % (f, extension))
        batch_futures[0] = compute_descript(smiles, output_file, bad_file, save_csv, ignore_3D)
    else:
        for i in range(0, len(smiles), chunksize):
            start_num = i + args.offset
            output_file = os.path.join(output_dir, "%s-%s-%s.%s" % (f, start_num, start_num+chunksize, extension))
            if bad_dir:
                bad_file = os.path.join(bad_dir, "%s-%s-%s.%s" % (f, start_num, start_num+chunksize, extension))
            batch_futures[(i,i+chunksize)] = compute_descript(smiles[i:i+chunksize], output_file, bad_file, save_csv, ignore_3D)
    
    print("[Main] Waiting for {} futures...".format(len(batch_futures)))
    for i in batch_futures:
        try:
            x = batch_futures[i].result()
            print(f"Chunk {i} is done with {x}")
        except Exception as e:
            print(f"Exception : {e}")
            print(f"Chunk {i} failed")

    print("[Main] All done! %s" % (time.time() - start))
