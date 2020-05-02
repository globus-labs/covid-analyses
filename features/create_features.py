import parsl
import os
import argparse
import re
import time
from pathlib import Path
import math
from fingerprints.fingerprints import compute_fingerprints
from images.images import compute_images
from descriptors.descriptors import compute_descriptors
from neural_fingerprints.neural_fingerprints import compute_neural_fingerprints

def generate_batch(filename, start=0, batchsize=10, max_batches=10):
    counter = 0
    if max_batches == 0:
        max_batches = 999999999

    with open(filename) as current:
        while current and counter < max_batches:
            yield current.tell()
            counter += 1
            for i in range(batchsize):
                x = current.readline()
                if not x:
                    current = None
                    break

def create_file_names(output_dir, name, start, finish, bad_dir, extension):
     output_file = os.path.join(output_dir, "%s-%s-%s.%s" % (name, start, finish, extension))
     bad_file = None
     if bad_dir:
         bad_file = os.path.join(bad_dir, "%s-%s-%s.%s" % (name, start, finish, extension))
     return (output_file, bad_file)

def create_dirs(dir_name):
    if dir_name: 
        try:
            os.makedirs(dir_name, exist_ok = True)
            print("Directory '%s' created" % dir_name)
        except OSError as error:
            print("Directory '%s' can not be created" % dir_name)
            exit()

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--debug", action='store_true', help="Enables debug logging")
    parser.add_argument("-t", "--type", default=None, help="Feature type: fingerprints, neural_fingerprints, descriptors, images", required=True)
    parser.add_argument("-i", "--input_file", default=None, help="input directory of smiles to process", required=True)
    parser.add_argument("-o", "--output_dir", default="outputs", help="Output directory. Default : outputs", required=True)
    parser.add_argument("-bo", "--bad_output_dir", default=None, help="Output directory for bad smile list")
    parser.add_argument("-b", "--batch_size", default=0, help="Batch size. Default: 0", type=int)
    parser.add_argument("-n", "--num_smiles", default=0, help="Number of smiles to process (for testing)")
    parser.add_argument("-csv", "--csv", default=False, action='store_true',  help="Save output as CSV. Default is to save as Pickle.")
    parser.add_argument("-gz", "--gzip", default=False, action='store_true',  help="Gzip output files. Default False.")
    parser.add_argument("-ow", "--overwrite", default=False, action='store_true',  help="Overwrite output files. Default False.")
    parser.add_argument("-wr", "--worker_read", default=False, action='store_true', help="Read smiles file on worker. Default: False.")
    parser.add_argument("-w", "--walltime", default=0, help="Walltime limit for running a batch (in seconds). 0 is no walltime limit. Default: 0.", type=int)
    parser.add_argument("-ig", "--ignore3d", default=False, action='store_true', help="Ignore 3D descriptors. True: 1613 descriptors; False: 1826 descriptors. Default False.")
    parser.add_argument("-m", "--model_file", default=None,  help="Model file for neural network fingerprints. Default None.")

    parser.add_argument("-c", "--config", default="local", help="Parsl config defining the target compute resource to use. Default: local")
    args = parser.parse_args()

 
    if args.config == "local":
        from parsl.configs.htex_local import config
        from parsl.configs.htex_local import config
        config.executors[0].label = "Foo"
        config.executors[0].max_workers = 5
    elif args.config == "theta":
        from configs.theta import config
    elif args.config == "cvd":
        from configs.theta_cvd import config

    parsl.load(config)
    
    output_dir = args.output_dir
    input_file = args.input_file
    bad_dir = args.bad_output_dir
    bad_file = None
    batch_size = int(args.batch_size)
    num_smiles = int(args.num_smiles)
    save_csv = args.csv
    save_gzip = args.gzip
    overwrite = args.overwrite
    worker_read = args.worker_read
    compute_type = args.type

    extension = "csv" if save_csv else 'pkl'
    extension = "%s.gz" % extension if save_gzip else extension

    extra_args = {}   
    if compute_type.startswith('fingerprint'):
        process_function = compute_fingerprints
    elif compute_type.startswith('image'):
        process_function = compute_images 
        if save_csv:
            print("No CSV support for images")
            exit()
    elif compute_type.startswith('descriptor'):
        process_function = compute_descriptors
        extra_args = {"ignore_3D" : args.ignore3d}
    elif compute_type.startswith('neural'):
        process_function = compute_neural_fingerprints
        extra_args = {"model_file" : args.model_file}
    else: 
        print("Type must be one of fingerprints, descriptors, images")
        exit()

    if worker_read and batch_size == 0:
        print("Must set a batch size for worker read")
        exit()

    if args.walltime > 0: 
        extra_args = {'walltime' : args.walltime}

    print(f"[Main] Computing {compute_type}") 
    print(f"[Main] Processing {input_file}")
    print(f"[Main] SaveCSV: {save_csv}, Gzip: {save_gzip}")
    print(f"[Main] WorkerRead: {worker_read}, Overwrite: {overwrite}")
    print(f"[Main] BatchSize: {batch_size}, Smiles: {num_smiles}")
    print(f"[Main] Output: {output_dir}")

    if extra_args:
        print(f"[Main] Extra args: {extra_args}")
   
    create_dirs(output_dir)
    create_dirs(bad_dir)

    f = Path(input_file).stem
 
    start = time.time() 
    batch_futures = {}

    if worker_read:
        # send file offset locations to each worker
        batch_generator = generate_batch(input_file, start=0, batchsize=batch_size, max_batches=0)
        i = 0
        for batch_index in batch_generator:
             print("Starting batch index: %s %s %s" % (batch_index, i, batch_size))
             output_file, bad_file = create_file_names(output_dir, f, i, i+batch_size, bad_dir, extension)

             batch_futures[(i,i+batch_size)] = process_function(smiles_file=input_file, start_index=batch_index, 
                                                      batch_size=batch_size, out_file=output_file, bad_file=bad_file, 
                                                      save_csv=save_csv, overwrite=overwrite, save_gzip=save_gzip, **extra_args)

             i += batch_size
    else:
        # Read in all the smiles
        with open(input_file) as ff:
            if int(args.num_smiles) > 0:
                smiles = ff.readlines()[:int(num_smiles)]           
            else: 
                smiles = ff.readlines()
        
        if batch_size == 0 or batch_size > len(smiles):
            print("Processing file: %s" % input_file)
            output_file, bad_file = create_file_names(output_dir, f, 0, len(smiles), bad_dir, extension) 
            batch_futures[0] = process_function(smiles=smiles, out_file=output_file, bad_file=bad_file, save_csv=save_csv, 
                                                overwrite=overwrite, save_gzip=save_gzip, **extra_args)
        else:
            print("Processing %s batches" % (math.ceil(len(smiles)/batch_size)))
            for i in range(0, len(smiles), batch_size):
                output_file, bad_file = create_file_names(output_dir, f, i, i+batch_size, bad_dir, extension)
                batch_futures[(i,i+batch_size)] = process_function(smiles=smiles[i:i+batch_size], out_file=output_file, bad_file=bad_file, 
                                                                   save_csv=save_csv, overwrite=overwrite, save_gzip=save_gzip, **extra_args)

 
    print("[Main] Waiting for {} futures...".format(len(batch_futures)))
    for i in batch_futures:
        try:
            x = batch_futures[i].result()
            print(f"Chunk {i} is done with {x}")
        except Exception as e:
            print(f"Exception : {e}")
            print(f"Chunk {i} failed")

    print("All done! %s" % (time.time() - start))
    



