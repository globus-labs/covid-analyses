import sys
import argparse
import os
import time
import re
import shutil

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_dir", required=True, default=None, help="input directory of smiles to process")
    parser.add_argument("-o", "--output_dir", required=True, default="outputs", help="Output directory. Default : outputs")
    parser.add_argument("-d", "--dataset", required=True, default="outputs", help="Output directory. Default : outputs")
    parser.add_argument("-off", "--offset", required=True, help="Offset")
    args = parser.parse_args()
    
    offset  = int(args.offset)

    out_dir = args.output_dir
    missing_dir = (os.path.join(args.output_dir, 'missing'))

    try: 
        os.makedirs(missing_dir, exist_ok = True)
    except: 
        pass
  

def find_top(input_dir):
    for d in os.listdir(input_dir):
        chunk = int(d[-2:])
        top = 0
        for f in os.listdir(os.path.join(input_dir, d)):
            try:
                result = re.findall('(-[\d]+-[\d]+)', str(f))
                _, s, e = result[-1].split('-')
                if int(e) > top: 
                    top = int(e)
            except: 
                if f != "missing":
                    print(f)
        print("%s: %s  (offset: %s)"  % (d, top, chunk*offset))
 

def get_to_name(dataset, f, offset):
    try:
        result = re.findall('(-[\d]+-[\d]+)', str(f))
        _, s, e = result[-1].split('-')
        return "%s-%s-%s.csv" % (dataset, int(s) + offset, int(e) + offset)
    except:
        return None

def copy_files(dataset, input_dir, output_dir, copy):
    for d in os.listdir(input_dir):
        chunk = int(d[-2:])
        chunk_offset = chunk *  offset
        for f in os.listdir(os.path.join(input_dir, d)):
            if f == 'missing':
                print ("missing")
                for m in os.listdir(os.path.join(input_dir, d, f)):
                    in_file = os.path.join(input_dir, d, f, m)
                    o_file = get_to_name(dataset, m, chunk_offset)
                    if o_file:
                        out_file = os.path.join(output_dir, f, o_file) 
                        print ("%s --> %s" % (in_file, out_file))
                        if copy: 
                            shutil.copyfile(in_file, out_file)
            else:
                o_file = get_to_name(dataset, f, chunk_offset) 
                in_file = os.path.join(input_dir, d, f) 
                if o_file:
                    out_file = os.path.join(output_dir, o_file)
                    print ("%s --> %s" % (in_file, out_file))
                    if copy:
                        shutil.copyfile(in_file, out_file)

copy_files(args.dataset, args.input_dir, args.output_dir, True)
#find_top(args.input_dir)
