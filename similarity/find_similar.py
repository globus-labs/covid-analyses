# from candle_apps.candle_apps import run_intranode
# from candle_apps.candle_apps import ModelInferer
#from candle_apps.candle_node_local import run_local
#import pandas as pd
import os
import glob
import argparse
import traceback
import parsl
import pickle
import gc
from itertools import islice
from rdkit import Chem
from rdkit.Chem import AllChem
from operator import itemgetter
#from targets import target_smiles
import logging

def set_file_logger(filename: str, name: str = 'runner', level: int = logging.DEBUG, format_string = None):
    """Add a stream log handler.

    Args:
        - filename (string): Name of the file to write logs to
        - name (string): Logger name
        - level (logging.LEVEL): Set the logging level.
        - format_string (string): Set the format string

    Returns:
       -  None
    """
    if format_string is None:
        format_string = "%(asctime)s.%(msecs)03d %(name)s:%(lineno)d [%(levelname)s]  %(message)s"

    logger = logging.getLogger(name)
    logger.setLevel(logging.DEBUG)
    handler = logging.FileHandler(filename)
    handler.setLevel(level)
    formatter = logging.Formatter(format_string, datefmt='%Y-%m-%d %H:%M:%S')
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    return logger


@parsl.python_app
def process_one_target(file, targets, top_n_matches, outfile=None):
    import rdkit
    import base64
    from rdkit import RDLogger
    from rdkit import DataStructs
    from pstats import SortKey
    import pickle
    from operator import itemgetter
    import csv
    
    target_results = {}

    read = pickle.load( open(file, 'rb') )

    fingerprint_set = []
    for sm, identifier, fp in read:
        # Next TRY here as some do not convert
        try:
            bv = DataStructs.ExplicitBitVect(base64.b64decode(fp))
        except:
            bv = None
        fingerprint_set += [(sm, bv, identifier)]

    for smile_target in targets:
        best_so_far = [('', 0.0) for index in range(top_n_matches)]
        bit_target = targets[smile_target]
            
        # Find scores for non-None fingerprints in fingerprint set
        scores = []
        for (smile, fingerprint, identifier) in fingerprint_set:
            try:
                score = DataStructs.TanimotoSimilarity(fingerprint, bit_target)
                scores += [(smile, score, identifier)]
            except:
                pass

        sorted_scores = sorted(scores, key=itemgetter(1))
        new_list = sorted_scores[-top_n_matches:]
        target_results[smile_target] = new_list

    if outfile:
        with open(outfile, 'wb') as f:
            pickle.dump(target_results, f)
        return outfile
    else:
        return(target_results)


def launch_slice(all_pickle_files, target_subset, prefix, wait=True):

    target_futures = []

    count  = 0
    for pickle_file in all_pickle_files:
        count += 1 
        if not pickle_file.endswith('.pkl'):
            logger.warning(f"Ignoring {pickle_file} not smile file")
            continue

        logger.info("Processing pickle_file: {} Target:{}".format(pickle_file, 
                                                                  prefix))

        outfile =  "{}/{}.{}.pkl".format(args.outdir, 
                                         os.path.basename(pickle_file).strip('.pkl'),
                                         prefix)                                  
        if os.path.exists(outfile) and os.stat(outfile).st_size > 0:
            logger.debug(f"{outfile} already exists")
            continue

        x = process_one_target(pickle_file,
                               target_subset,
                               top_n_matches = int(args.top_n_matches),
                               outfile=outfile
                           )
        target_futures.append(x)
        # Wait for all futures
        logger.info(f"Waiting for all futures for {prefix}")

    count = 0
    
    if wait is True:
        outfiles = [fu.result() for fu in target_futures]
        print(f"Output files : {outfiles}")
        print(f"{prefix} done")
    else:
        outfiles = target_futures

    return outfiles

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-v", "--version",
                        help="Print Endpoint version information")
    parser.add_argument("-d", "--debug", action='store_true',
                        help="Enables debug logging")
    parser.add_argument("-n", "--name", required=True,
                        help="Name to use in csv")
    parser.add_argument("-s", "--smile_glob", default=".",
                        help="Glob pattern that points to all .smi smile files")

    parser.add_argument("--target_glob", required=True,
                        help="Target glob")
    
    parser.add_argument("--top_n_targets", default="1000",
                        help="Top N items to take from the target source csv file")

    parser.add_argument("--top_n_matches", default="100",
                        help="Top N matches to the target to calculate")

    parser.add_argument("-o", "--outdir", default="outputs",
                        help="Output directory. Default : outputs")
    parser.add_argument("-c", "--config", default="local",
                        help="Parsl config defining the target compute resource to use. Default: local")
    parser.add_argument("-l", "--logfile", default="runner.log",
                        help="log file path ")
    parser.add_argument("-p", "--pickle", action='store_true', default=False, 
                        help="Output data in pickled format. Default: False")
    args = parser.parse_args()

    #for smile_dir in glob.glob(args.smile_dir):
    #    print(smile_dir)
    #exit(0)

    logger = set_file_logger(args.logfile)

    if args.config == "local":
        from parsl.configs.htex_local import config
        from parsl.configs.htex_local import config
        config.executors[0].label = "Login"
        # config.executors[0].max_workers = 4
    elif args.config == "theta":
        from theta import config
    elif args.config == "frontera":
        from frontera import config
    elif args.config == "frontera_small":
        from frontera_small import config
    elif args.config == "theta_test":
        from theta_test import config
    elif args.config == "comet":
        from comet import config

    # Most of the app that hit the timeout will complete if retried.
    # but for this demo, I'm not setting retries.
    # config.retries = 2
    # print("*"*10, "Skip loading config", "*"*10)
    parsl.load(config)


    if args.debug:
        parsl.set_stream_logger()

    all_targets = {}
    for target_file in glob.glob(args.target_glob):        
        target_smiles = None
        target_name = os.path.basename(target_file).strip('_dock.csv').strip('top.7.5k.ml.')
        print("Target name : ", target_name)

        all_targets[target_name] = {}

        with open(target_file) as f:
            smiles = f.readlines()[:int(args.top_n_targets)]
            target_smiles = [smile.strip().split(',')[-1] for smile in smiles]

            for smile in target_smiles:
                try:
                    bit_target = AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(smile), 2, nBits=2048)
                    all_targets[target_name][smile] = bit_target
                except Exception as e:
                    logger.exception(f"Caught error processing {smile}")
                    print(f"{smile} failed to compute fingerprint")
                        
    logger.info("Targets computed")

    for source in all_targets:
        print(all_targets[source])
        print("Processing {}.count_{}.{}".format(args.name,
                                                 len(all_targets[source]),
                                                 source))
    # exit()

    os.makedirs(args.outdir, exist_ok=True)

    all_pickle_files = glob.glob(args.smile_glob)

    all_outfiles = []
    wait = False
    for target_name in all_targets:
        print("Launching against {}".format(target_name))
        outfiles = launch_slice(all_pickle_files, all_targets[target_name], prefix=target_name, wait=wait)                    
        all_outfiles.extend(outfiles)
        print("Garbage collect : ", gc.collect())
        #break


    for fu in all_outfiles:
        outfile = None
        if wait is False:
            try:
                outfile = fu.result()
            except:
                pass
        else:
            outfile = fu
        print("Outfile : {}".format(outfile))
    
    print("All done!")

