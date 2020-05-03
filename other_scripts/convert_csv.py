import pickle
import csv
import sys
from parsl.app.app import python_app
import argparse
import parsl
from pathlib import Path
import os
import time

@python_app
def convert_file(dataset_name, input_file, output_file, bad_file):
    import csv
    import pickle

    pickle_none = pickle.dumps(None)
    csv_results = []
    bad = []
    pickle_list = pickle.load(open(input_file, 'rb'))
    for s, d in pickle_list.items():
        if d[1] == pickle_none:
           bad.append((dataset_name, ''.join(d[0]), s))
        else:
            data = [dataset_name, ''.join(d[0]), s] + list(d[1])
            data = ['' if str(i) == 'nan' else i for i in data]
            csv_results.append(data)
        #count += 1

    with open(output_file, 'w', newline='') as o_file:
        writer = csv.writer(o_file, delimiter=',', quoting=csv.QUOTE_MINIMAL)
        writer.writerows(csv_results)
    with open(bad_file, 'w', newline='') as b_file:
        writer = csv.writer(b_file, delimiter=',', quoting=csv.QUOTE_MINIMAL)
        writer.writerows(bad)

    return output_file


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--dataset_name", default=None, help="Dataset name")
    parser.add_argument("-i", "--input_dir", default=None, help="input directory of smiles to process")
    parser.add_argument("-o", "--output_dir", default="outputs", help="Output directory. Default : outputs")
    parser.add_argument("-ow", "--overwrite", default=False, action='store_true',  help="Overwrite output files. Default False.")
    parser.add_argument("-c", "--config", default="local", help="Parsl config defining the target compute resource to use. Default: local")
    parser.add_argument("-s", "--subdirs", default=False, action='store_true', help="Use subdirs. Default: False")
    args = parser.parse_args()


    if args.config == "local":
        from parsl.configs.htex_local import config
        from parsl.configs.htex_local import config
        config.executors[0].label = "Foo"
        config.executors[0].max_workers = 4
    elif args.config == "cvd":
        from configs.theta_cvd import config

    parsl.load(config)


    try:
        os.makedirs(args.output_dir)
    except:
        pass

    try:
        os.makedirs(os.path.join(args.output_dir, 'missing'))
    except:
        pass

    batch_futures = []
    start = time.time()
    if args.subdirs:
        for d in os.listdir(args.input_dir):
            for f in os.listdir(os.path.join(args.input_dir, d)):
                try:
                    os.makedirs(os.path.join(args.output_dir, d))
                    os.makedirs(os.path.join(args.output_dir, d, 'missing'))
                except:
                    pass
                if f.endswith('.pkl'):
                    input_file = os.path.join(args.input_dir, d, f)
                    output_file = os.path.join(args.output_dir, d, '%s%s' % (Path(f).stem, '.csv'))
        
                    bad_file = os.path.join(args.output_dir, d, 'missing', '%s%s' % (Path(f).stem, '.csv'))
                    #print(input_file)
                    #print(output_file)
                    batch_futures.append(convert_file(args.dataset_name, input_file, output_file, bad_file))
    else:
        for f in os.listdir(args.input_dir):
            if f.endswith('.pkl'):
                input_file = os.path.join(args.input_dir, f)
                output_file = os.path.join(args.output_dir, '%s%s' % (Path(f).stem, '.csv'))
                bad_file = os.path.join(args.output_dir, 'missing', '%s%s' % (Path(f).stem, '.csv'))
                print(input_file)
                #print(output_file)

                batch_futures.append(convert_file(args.dataset_name, input_file, output_file, bad_file))

    print("[Main] Waiting for {} futures...".format(len(batch_futures)))
    for i in batch_futures:
        try:
            x = i.result()
            print(f"Chunk {i} is done with {x}")
        except Exception as e:
            print(f"Exception : {e}")
            print(f"Chunk {i} failed")

    print("All done! %s" % (time.time() - start))

