import os
import argparse
import parsl
from parsl.app.app import python_app, bash_app
import subprocess


split_prefix = "split-"


def split(input_path, filename, batch_size):
    cmd = f"cd {input_path};split -d -l {batch_size} {filename} {split_prefix}{filename.replace('.csv', '')}"
    subprocess.check_output(cmd, shell=True)


def combine(input_path, filename):
    cmd = f"cd {input_path};cat "
    for f in os.listdir(input_path):
        if f.startswith(f"new_{split_prefix}"):
            cmd += f + " "
    cmd += f" > new_can_{filename}"
    subprocess.check_output(cmd, shell=True)
        

@python_app
def can(input_path, filename, recan_path):
    import subprocess
    cmd = f"bash {recan_path} {input_path} {filename}"
    try:
        res = subprocess.check_output(cmd, shell=True)
    except subprocess.CalledProcessError as e:
        return e.returncode, e.output
    return cmd


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_file", required=True,
                        help="File path to the input smiles csv file")
    parser.add_argument("-b", "--batch_size", default=500000, type=int,
                        help="Size of each split. Default=0.")
    parser.add_argument("-r", "--recan_path", default=f"{os.getcwd()}/recan-new.sh",
                        help="File path to the input smiles csv file")
    parser.add_argument("-c", "--config", default="local",
                        help="Parsl config defining the target compute resource to use. Default: local")
    args = parser.parse_args()

    if args.config == "local":
        from parsl.configs.htex_local import config
        from parsl.configs.htex_local import config
        config.executors[0].label = "Foo"
        config.executors[0].max_workers = 4
    elif args.config == "theta":
        from theta import config
    elif args.config == "comet":
        from comet import config
    elif args.config == "cvd":
        from configs.theta_cvd import config

    parsl.load(config)

    input_file = args.input_file
    batch_size = args.batch_size

    basepath = os.path.dirname(input_file)
    filename = os.path.basename(input_file)

    print(f"Basepath: {basepath}. filename: {filename}")
    
    # Split the smiles
    split(basepath, filename, args.batch_size)

    # Can with python apps
    can_futures = []
    for f in os.listdir(basepath):
        if not f.startswith(split_prefix):
            continue
        can_future = can(basepath, f, args.recan_path)
        can_futures.append(can_future)
    
    for fut in can_futures:
        print(fut.result())
    print("Finished waiting for all futures, combining")

    # Combine afterwards, can be commented out
    combine(basepath, filename)
    
    

