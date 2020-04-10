# Canonicalize

The canonicalize tool uses Open Babel to convert input SMILES into canonical SMILES. 
The tool can be configured to select a particular column to be used for SMILES and molecule 
identifiers. 

## Example 
The script can be run as follows to convert a tab-separated input file into one or more files with canonical SMILES in the output directory.
```
$ python canonicalize.py --input_file <input_file> --output_dir <output_dir>
```

## Usage
```
usage: canonicalize.py [-h] -f INPUT_FILE -o OUTPUT_DIR [-s SMILE_COL]
                       [-id IDENTIFIER_COL] [-b BATCH_SIZE] [--header]
                       [-d DELIMITER] [-r REMOVE]

optional arguments:
  -h, --help            show this help message and exit
  -f INPUT_FILE, --input_file INPUT_FILE
                        The smile file to process
  -o OUTPUT_DIR, --output_dir OUTPUT_DIR
                        The output directory to store canonicalized smiles
  -s SMILE_COL, --smile_col SMILE_COL
                        The column to find the smile file in
  -id IDENTIFIER_COL, --identifier_col IDENTIFIER_COL
                        The column to find the identifier in. None if none
                        exist.
  -b BATCH_SIZE, --batch_size BATCH_SIZE
                        Num smiles per file. Set to 0 for all smiles in file
  --header              Whether the file has a header that need removing
  -d DELIMITER, --delimiter DELIMITER
                        How smiles are separated
  -r REMOVE, --remove REMOVE
                        A string to remove from the smile column (e.g., quotes
                        or <value>
```
