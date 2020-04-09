import argparse
import os
import parsl
import pickle
import time
from pathlib import Path
from parsl.app.app import python_app

@python_app
def compute_descript(smiles, out_file, bad_file=None):
    from mordred import Calculator, descriptors
    from rdkit import Chem
    import numpy as np
    import pickle
    import csv
    bad = []
    results = {}
    mol = None

    calc = Calculator(descriptors, ignore_3D=False) 

    for smile_tuple in smiles:
        s = smile_tuple.split(',')
        smile = s[2].rstrip()
        identifier = s[1]
        
        try:
            mol = Chem.MolFromSmiles(smile)

            if mol is None or len(mol.GetAtoms()) > 100:
                print("Error processing mol")
                bad.append(smile_tuple)
            else:
                descs = calc(mol)

                # Mods from Xuefeng
                descs.fill_missing("nan")
            
                data = np.array(descs).flatten().astype(np.float32)
                results[smile] = ([identifier], data)
        except:
            bad.append(smile_tuple)

    with open(out_file, 'wb') as output_file:
        pickle.dump(results, output_file, protocol=pickle.HIGHEST_PROTOCOL)

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
                        help="Size of the batch of smiles to send to each node for processing. Default=4, should be 10K")
    parser.add_argument("-o", "--output_dir", default="outputs",
                        help="Output directory. Default : outputs")
    parser.add_argument("-off", "--offset", default=0, type=int,
                        help="Offset for starting processing from separate files. Default=0")
    parser.add_argument("-c", "--config", default="local",
                        help="Parsl config defining the target compute resource to use. Default: local")
    args = parser.parse_args()

    if args.config == "local":
        from parsl.configs.htex_local import config
        from parsl.configs.htex_local import config
        config.executors[0].label = "Foo"
        config.executors[0].max_workers = 1
    elif args.config == "theta":
        from theta import config
    elif args.config == "comet":
        from comet import config

    # config.retries = 2
    parsl.load(config)

    if args.debug:
        parsl.set_stream_logger()

    output_dir = args.output_dir
    input_file = args.input_file
    bad_dir = args.bad_output_dir
    bad_file = None
    chunksize = int(args.batch_size)

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

    chunksize = int(args.batch_size)
    num_smiles = len(smiles)
    print(f"[Main] Chunksize : {chunksize}, Smiles: {num_smiles}")
    batch_futures = {}
   
    f = Path(input_file).stem 

    start = time.time()

    if chunksize == 0 or chunksize > len(smiles):
        print("Processing file: %s" % input_file)
        output_file = os.path.join(output_dir, f)
        if bad_dir:
            bad_file = os.path.join(bad_dir, f)
        batch_futures[0] = compute_descript(smiles, output_file, bad_file)  
    else:
        for i in range(0, len(smiles), chunksize):
            start_num = i + args.offset
            output_file = os.path.join(output_dir, "%s-%s-%s.pkl" % (f, start_num, start_num+chunksize))
            if bad_dir:
                bad_file = os.path.join(bad_dir, "%s-%s-%s.pkl" % (f, start_num, start_num+chunksize))
            batch_futures[(i,i+chunksize)] = compute_descript(smiles[i:i+chunksize], output_file, bad_file)
    
    print("[Main] Waiting for {} futures...".format(len(batch_futures)))
    for i in batch_futures:
        try:
            x = batch_futures[i].result()
            print(f"Chunk {i} is done with {x}")
        except Exception as e:
            print(f"Exception : {e}")
            print(f"Chunk {i} failed")

    print("All done! %s" % (time.time() - start))
