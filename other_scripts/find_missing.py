import sys
import argparse
import os
import time
import re
import shutil

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_dir", required=True, default=None, help="input directory of smiles to process")
    parser.add_argument("-b", "--batch_size",  type=int, default=10000, help="batch size")
    parser.add_argument("-pm", "--print_missing",  action='store_true', default=False, help="print missing files")
    parser.add_argument("-e", "--extension",  default='csv', help="file extension")
    args = parser.parse_args()

    input_dir = args.input_dir
    batch_size = args.batch_size

    found = []
    
    for f in os.listdir(os.path.join(input_dir)):
        if f.endswith(args.extension):
           result = re.findall('(-[\d]+-[\d]+)', str(f))
           _, s, e = result[-1].split('-')
           found.append(int(s))
      
    found.sort()
    max_num = found[-1]
    missing = []
    for i in range(0, int(max_num/batch_size)):
        if (i*batch_size) not in found: 
            missing.append((i*batch_size))

    print("Total files: %s" % len(found))
    print("Highest chunk: %s" % max_num)
    print("Missing: %s" % len(missing))

    if args.print_missing:
        print(missing)
