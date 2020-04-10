# Create Descriptors
This script uses mordred to compute molecular descriptors. 


## Running the workflow

Usage::
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
                        processing. Default=4, should be 10K
  -o OUTPUT_DIR, --output_dir OUTPUT_DIR
                        Output directory. Default : outputs
  -off OFFSET, --offset OFFSET
                        Offset for starting processing from separate files.
                        Default=0
  -c CONFIG, --config CONFIG
                        Parsl config defining the target compute resource to
                        use. Default: local
