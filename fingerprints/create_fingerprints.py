import parsl
import os
from parsl.app.app import python_app, bash_app
import argparse
import re
import time
from pathlib import Path

@python_app
def create_fingerprints(smiles, out_file, bad_file=None, save_csv=False):
    import os
    import logging
    import pickle
    from rdkit import Chem
    from rdkit.Chem import AllChem
    import csv

    sep = ","
    results = []
    bad = []
    for s in smiles:
        #mol = s.split(sep)i
        mol = None
        try:
            mol = s.split(sep)
            molecule = mol[2].rstrip()
            identifier = mol[1]
            fp = AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(molecule), 2, nBits=2048)
            if save_csv:
                results.append((molecule, identifier, fp.ToBase64())) 
            else:
                results.append((molecule, identifier, fp))
        except:
             fp = None
             bad.append(mol)

    if save_csv: 
        with open(out_file, 'w') as output_file:
            writer = csv.writer(output_file, delimiter=',') # quoting=csv.QUOTE_MINIMAL)
            writer.writerows(results)
        if bad_file and len(bad) > 0:
            with open(bad_file, 'w') as b_file:
                b_writer = csv.writer(b_file, delimiter=',') # quoting=csv.QUOTE_MINIMAL)
                b_writer.writerows(bad)
    else:
        with open(out_file, 'wb') as output_file:
            pickle.dump(results, output_file, protocol=pickle.HIGHEST_PROTOCOL)
    
        if bad_file and len(bad) > 0:
            with open(bad_file, 'wb') as b_file:
                pickle.dump(bad, b_file, protocol=pickle.HIGHEST_PROTOCOL)

    return out_file


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--debug", action='store_true',
                        help="Enables debug logging")
    parser.add_argument("-i", "--input_file", default=None,
                        help="input directory of smiles to process")
    parser.add_argument("-o", "--output_dir", default="outputs",
                        help="Output directory. Default : outputs")
    parser.add_argument("-bo", "--bad_output_dir", default=None,
                        help="Output directory for bad smile list")
    parser.add_argument("-b", "--batch_size", default=0,
                        help="Batch size. Default: 0")
    parser.add_argument("-n", "--num_smiles", default=0, 
                        help="Number of smiles to process (for testing)")
    parser.add_argument("-off", "--offset", default=0, type=int,
                        help="Offset to start numbering")
    parser.add_argument("-csv", "--csv", default=False, action='store_true',
                        help="Save output as CSV. Default is to save as Pickle.")
    parser.add_argument("-c", "--config", default="local",
                        help="Parsl config defining the target compute resource to use. Default: local")
    args = parser.parse_args()

 
    if args.config == "local":
        from parsl.configs.htex_local import config
        from parsl.configs.htex_local import config
        config.executors[0].label = "Foo"
        config.executors[0].max_workers = 2
    elif args.config == "theta":
        from configs.theta import config
    elif args.config == "cvd":
        from configs.theta_cvd import config

    parsl.load(config)
    
    output_dir = args.output_dir
    input_file = args.input_file
    bad_dir = args.bad_output_dir
    bad_file = None
    chunksize = int(args.batch_size)

    save_csv = args.csv
    extension = "csv" if save_csv else 'pkl'

    if not os.path.isdir(output_dir):
        try:
            os.makedirs(output_dir, exist_ok = True)
            print("Output directory '%s' created" % output_dir)
        except OSError as error:
            print("Output directory '%s' can not be created" % output_dir)
            exit()
    if bad_dir and not os.path.isdir(bad_dir):
        try:
            os.makedirs(bad_dir, exist_ok = True)
            print("bad output directory '%s' created" % bad_dir)
        except OSError as error:
            print("bad output directory '%s' can not be created" % bad_dir)
            exit()

    start = time.time()
    batch_futures = {}
    # Process each file in the dataset
    #for f in os.listdir(input_dir):i

    #input_file = os.path.join(input_dir, f)
    with open(input_file) as ff:
        if int(args.num_smiles) > 0:
            smiles = ff.readlines()[:int(args.num_smiles)]
        else:
            smiles = ff.readlines()

    f = Path(input_file).stem
    batch_futures = {}
    if chunksize == 0 or chunksize > len(smiles):
        print("Processing file: %s" % input_file)
        output_file = os.path.join(output_dir, '%s.%s' % (f, extension))
        if bad_dir:
            bad_file = os.path.join(bad_dir, "%s.%s" % (f, extension))
        batch_futures[0] = create_fingerprints(smiles, output_file, bad_file, save_csv)
    else: 
        for i in range(0, len(smiles), chunksize):
            start_num = i + args.offset # start num for naming file
            output_file = os.path.join(output_dir, "%s-%s-%s.%s" % (f, start_num, start_num+chunksize, extension))
            if bad_dir: 
                bad_file = os.path.join(bad_dir, "%s-%s-%s.%s" % (f, start_num, start_num+chunksize, extension))
            batch_futures[(i,i+chunksize)] = create_fingerprints(smiles[i:i+chunksize], output_file, bad_file, save_csv)

    
    print("[Main] Waiting for {} futures...".format(len(batch_futures)))
    for i in batch_futures:
        try:
            x = batch_futures[i].result()
            print(f"Chunk {i} is done with {x}")
        except Exception as e:
            print(f"Exception : {e}")
            print(f"Chunk {i} failed")

    print("All done! %s" % (time.time() - start))
    



