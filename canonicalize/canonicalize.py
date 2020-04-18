import subprocess
import argparse
import math
import sys
import os
from subprocess import PIPE
import pandas as pd
import numpy as np

def prepare_files(work_dir, filename, smile, iden, label, header, delim, size, remove, out_delim):
    """Prepare the file for can'ing. Because the files can be in excess of 100GB
     we use command line tools and subprocess.

    This uses the following process:
    - Read the file
    - Remove the header if specificed
    - Remove duplicates
    - Split the file into 'size' line chunks
    - Return the list of filenames"""

    preped_files = []
    print('starting prep')

    # load the data into pandas
    df = pd.read_csv(filename, sep=delim, header=None)

    # remove the header
    if header:
        df = df[1:]
    df = df.fillna('-1')

    if smile is not None or iden is not None:
        to_keep = pd.DataFrame()
        if smile is not None:
            to_keep[smile] = df[smile]
        if iden is not None:
            try:
                # remove .0 from int ids
                to_keep[iden] = df[iden].astype(np.int64)
            except:
                pass
        if label is not None:
            to_keep[label] = df[label] 
        df = to_keep
    print(df)
    # Switch order of columns. OpenBabel wants SMILES first.
    #if iden and iden < smile:
    #    print("Switching order of SMILES and identifier")
    #    df = df[[2,1]]
    print(df)
    # remove anything we need to, e.g., quotes or </value>
    if remove:
        df[0] = df[0].str.replace(remove, '')
    # drop duplicates
    df.drop_duplicates(subset=None, inplace=True)

    print(df.shape)
    print('writing output')
    # write the file back out as size chunks
    num_files = 1
    if size > 0:
        num_files = math.ceil(df.shape[0] / size)
    print(f'num files: {num_files}')

    for id, df_i in enumerate(np.array_split(df, num_files)):
        out_file = f"{work_dir}/smiles-{id}.smi"
        print(out_file)
        df_i.to_csv(out_file, sep=out_delim, index=False, header=False)
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
    parser.add_argument("-f", "--input_file", default=None, required=True,
                        help="The smile file to process")
    parser.add_argument("-o", "--output_dir", default=None, required=True,
                        help="The output directory to store canonicalized smiles")
    parser.add_argument("-s", "--smile_col", type=int, default=None,
                        help="The column to find the smile file in")
    parser.add_argument("-id", "--identifier_col", type=int, default=None,
                        help="The column to find the identifier in. None if none exist.")
    parser.add_argument("-ld", "--label_col", type=int, default=None,
                        help="The column to find the label in. None if none exist.")
    parser.add_argument("-b", "--batch_size", default="0",
                        help="Num smiles per file. Set to 0 for all smiles in file")
    parser.add_argument("--header",  action='store_true',
                        help="Whether the file has a header that need removing")
    parser.add_argument("-d", "--delimiter", default="\t",
                        help="How smiles are separated")
    parser.add_argument("-od", "--out_delimiter", default="\t",
                        help="How output smiles are separated")
    parser.add_argument("-r", "--remove", default=None,
                        help="A string to remove from the smile column (e.g., quotes or <value>")
    args = parser.parse_args()

    num_smiles = int(args.batch_size)

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

    to_process = prepare_files(work_dir, args.input_file, args.smile_col, args.identifier_col, args.label_col,
                                        args.header, args.delimiter, num_smiles, args.remove, args.out_delimiter)

    caned_files = can_files(work_dir, to_process)

