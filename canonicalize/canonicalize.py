import subprocess
import argparse
import math
import sys
import os
from subprocess import PIPE
import pandas as pd
import numpy as np

def prepare_files(work_dir, filename, smile, iden, header, delim, size, remove):
    """Prepare the file for can'ing. Because the files can be in excess of 100GB
     we use command line tools and subprocess.

    This uses the following process:
    1. Unzip the file if it is zipped.
    2. Use awk to separate out the smile and id col.
    3. Sort -u the file to remove duplicates
    4. Remove the header line if set
    5. Split the file into 1M line chunks
    6. Return the list of filenames"""

    preped_files = []
    print('starting prep')

    # check for compressed

    # load the data into pandas
    print('reading file')
    df = pd.read_csv(filename, sep=delim, header=None)

    # remove the header
    if header:
        df = df[1:]
    df = df.fillna('-1')
    #print(df[4])

    print(smile)
    if smile is not None or iden is not None:
        print('trimming')
        to_keep = pd.DataFrame()
        if iden is not None:
            try:
                # remove .0 from int ids
                to_keep[iden] = df[iden].astype(np.int64)
            except:
                pass

        if smile is not None:
            to_keep[smile] = df[smile]  
        print(to_keep)
        df = to_keep

    print(df)
    # remove anything we need to, e.g., quotes or </value>
    if remove:
        df[0] = df[0].str.replace(remove, '')
    # drop duplicates
    print('dropping dupes')
    df.drop_duplicates(subset=None, inplace=True)

    print(df.shape)
    print(df.shape[0])
    print('writing output')
    # write the file back out as 1M line chunks
    num_files = math.ceil(df.shape[0] / size)
    print(f'num files: {num_files}')

    for id, df_i in enumerate(np.array_split(df, num_files)):
        out_file = f"{work_dir}/smiles-{id}.smi"
        print(out_file)
        df_i.to_csv(out_file, sep="\t", index=False, header=False)
        preped_files.append(out_file)

    print(preped_files)
    return preped_files


def can_files(work_dir, to_process):
    """Run openbabel on each of the files and return the list of files processed"""

    os.chdir((work_dir))

    canned_files = []

    # deal with it being a string file, like a list
    if not isinstance(to_process, list):
        to_process = [to_process]

    for fn in to_process:
        filename = fn.replace(work_dir + "/", "")
        print("Running open babel")
        cmd = f"obabel {filename} -O can_{filename} -ocan -e"
        canned_files.append(f"{work_dir}/can_{filename}")
        sp = subprocess.run(cmd, stdout=PIPE, stderr=PIPE,
                            shell=True, executable='/bin/bash')
    print(canned_files)
    return canned_files

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--filename", default=None,
                        help="The smile file to process")
    parser.add_argument("-o", "--output_dir", default=None,
                        help="The output directory to store canonicalized smiles")
    parser.add_argument("-s", "--smile_col", type=int, default=None,
                        help="The column to find the smile file in")
    parser.add_argument("-i", "--identifier_col", type=int, default=None,
                        help="The column to find the identifier in. None if none exist.")
    parser.add_argument("-n", "--size", default="500000000",
                        help="Num smiles per file")
    parser.add_argument("--header",  action='store_true',
                        help="Whether the file has a header that need removing")
    parser.add_argument("-d", "--delimiter", default="\t",
                        help="How smiles are separated")
    parser.add_argument("-r", "--remove", default=None,
                        help="A string to remove from the smile column (e.g., quotes or <value>")
    args = parser.parse_args()

    num_smiles = int(args.size)

    work_dir = args.output_dir

    if work_dir is None:
        work_dir = "/".join(args.filename.split("/")[:-1])

    if not os.path.isdir(work_dir):
        try: 
            os.makedirs(work_dir, exist_ok = True) 
            print("Output directory '%s' created" % work_dir) 
        except OSError as error: 
            print("Output directory '%s' can not be created" % work_dir) 
            exit()

    os.chdir(work_dir)

    to_process = prepare_files(work_dir, args.filename, args.smile_col, args.identifier_col,
                                        args.header, args.delimiter, num_smiles, args.remove)

    caned_files = can_files(work_dir, to_process)

