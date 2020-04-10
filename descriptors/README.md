#  Descriptors
The create descriptors tool uses mordred to compute molecular descriptors.  It takes an input CSV file formatted with lines containing <SOURCE, IDENTIFIER, SMILE> and creates one or more output Python Pickle files of the format: 

```
{ SMILE:  (
    [IDENTIFIERS], 
    NumPy Array[descriptor1, descriptor2, .., descriptorN]
    ),  
..
}
```

# Example

The following example will create a single Pickle output file in the output directory for an input CSV file.
```
$ python create_descriptors.py --input_file <INPUT> --output_dir <OUTPUT> 

```
# Usage

```
usage: create_descriptors.py [-h] [-v VERSION] [-d] [-n NUM_SMILES] -i
                             INPUT_FILE [-bo BAD_OUTPUT_DIR] [-b BATCH_SIZE]
                             [-o OUTPUT_DIR] [-off OFFSET] [-c CONFIG]

optional arguments:
  -h, --help            show this help message and exit
  -v VERSION, --version VERSION
                        Print Endpoint version information
  -d, --debug           Enables debug logging
  -n NUM_SMILES, --num_smiles NUM_SMILES
                        Number of smiles to load and run. Default=10000, if
                        set to 0 the entire file will be used
  -i INPUT_FILE, --input_file INPUT_FILE
                        File path to the input smiles csv file
  -bo BAD_OUTPUT_DIR, --bad_output_dir BAD_OUTPUT_DIR
                        Output directory for bad smile list
  -b BATCH_SIZE, --batch_size BATCH_SIZE
                        Size of the batch of smiles to send to each node for
                        processing. Batch size of 0 will process all smiles in
                        a single batch. Default=0.
  -o OUTPUT_DIR, --output_dir OUTPUT_DIR
                        Output directory. Default : outputs
  -off OFFSET, --offset OFFSET
                        Offset for starting processing from separate files.
                        Default=0
  -c CONFIG, --config CONFIG
                        Parsl config defining the target compute resource to
                        use. Default: local
 ```
