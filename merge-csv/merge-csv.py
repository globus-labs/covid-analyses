import argparse
import os
import pickle
import csv
from pathlib import Path

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-v", "--version",
                        help="Print Endpoint version information")
    parser.add_argument("-if", "--input_file", 
                        help="File path to the input smiles file")
    parser.add_argument("-id", "--input_dir",
                        help="File path to the input smiles csv file")
    parser.add_argument("-od", "--output_dir", default="outputs",
                        help="Output directory. Default : outputs")
    parser.add_argument("-of", "--output_file", default=None,
                        help="Output file. ")
    parser.add_argument("-l", "--label", default="",
                        help="Label to include in CSV")
    args = parser.parse_args()

    try:
        os.makedirs(args.output_dir)
    except:
        pass

    if args.input_file:
        with open(args.input_file) as f:
            smiles = f.readlines()   
    else: 
        smiles = []
        for f in os.listdir(args.input_dir):
            input_file = os.path.join(args.input_dir, f)
            with open(input_file) as of:
                smiles.extend(of.readlines())
    print("Total smiles: %s" % len(smiles) )
    
    f = Path(input_file).stem

    output_file = args.output_file
    if args.output_file is None:
        output_file = os.path.join(args.output_dir, "%s.csv" % f) 
   
    print("Output file: %s" % output_file) 
    with open(output_file, 'w') as out_csv:
        writer = csv.writer(out_csv, delimiter=",")
        
        for i in range(0, len(smiles)):
             s = smiles[i].split('\t')
             identifier = ''
             if len(s) > 1:
                 identifier = s[1].rstrip()
             writer.writerow((args.label, identifier, s[0].rstrip()))
     
    print("All done!")

