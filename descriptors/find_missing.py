# from candle_apps.candle_apps import run_intranode
# from candle_apps.candle_apps import ModelInferer
import pandas as pd
import argparse
import os

import parsl
import pickle



if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-v", "--version",
                        help="Print Endpoint version information")
    parser.add_argument("-d", "--debug", action='store_true',
                        help="Enables debug logging")
    parser.add_argument("-n", "--num_smiles", default=8,
                        help="Number of smiles to load and run. Default=10000, if set to 0 the entire file will be used")
    parser.add_argument("-s", "--smile_file", default="train.csv",
                        help="File path to the smiles csv file")
    parser.add_argument("-b", "--batch_size", default="10000",
                        help="Size of the batch of smiles to send to each node for processing. Default=4, should be 10K")
    parser.add_argument("-o", "--outdir", default="outputs",
                        help="Output directory. Default : outputs")
    parser.add_argument("-c", "--config", default="local",
                        help="Parsl config defining the target compute resource to use. Default: local")
    args = parser.parse_args()

    """
    if args.config == "local":
        from parsl.configs.htex_local import config
        from parsl.configs.htex_local import config
        config.executors[0].label = "Foo"
        config.executors[0].max_workers = 1
    elif args.config == "theta":
        from theta import config
    elif args.config == "comet":
        from comet import config
    # Most of the app that hit the timeout will complete if retried.
    # but for this demo, I'm not setting retries.
    # config.retries = 2
    parsl.load(config)
    parsl_runner = parsl.python_app(run_local)
    if args.debug:
        parsl.set_stream_logger()
    if args.num_smiles == "0":
        #print("[Main] Loading all data available")
        # smiles = pd.read_csv("train.csv").iloc[:,0].tolist()
        smiles = pd.read_csv(args.smile_file, error_bad_lines=False) # .iloc[:,0].tolist()
    else:
        print(f"[Main] Loading {args.num_smiles} smiles from file")
        smiles = pd.read_csv(args.smile_file,  error_bad_lines=False, nrows=int(args.num_smiles)) # .iloc[:,0].tolist()
    """

    with open(args.smile_file) as f:
        if int(args.num_smiles) > 0:
            smiles = f.readlines()[:int(args.num_smiles)]
        else:
            smiles = f.readlines()

    try:
        os.makedirs(args.outdir)
    except:
        pass
    #models_to_test = [ModelInferer(1613), ModelInferer(1613)]
    models_to_test = [1, 2]

    chunksize = int(args.batch_size)
    print(f"[Main] Chunksize : {chunksize}")
    batch_futures = {}

    for parentdir in os.listdir(args.outdir):

        outdir = f"{args.outdir}/{parentdir}"
        print(outdir)
        sourcefile = parentdir.replace('_descriptors', '.smi')

        i = 10000
        
        paths = ["{}/{}".format(outdir, p) for p in os.listdir(outdir)]
        # print(paths[:10])
        for i in range(0, 60*(10**6) , chunksize):
            out_file=f"{outdir}/{sourcefile}.chunk-{i}-{i+chunksize}.pkl"
            # if not os.path.exists(out_file):
            if out_file not in paths:
                print("Missing result file: {}".format(out_file))

    print("All done!")
