# Images

The create images tool will create 128*128 2D images of a molecule using RDKit. 
It consumes a CSV file formatted as <SOURCE,IDENTIFIER,SMILE> and will create one or more output Pickle files 
containg SMILES and Python Image Library (PIL) images. Specifically:

```
[(SOURCE, IDENTIFIER, SMILES, <PIL_IMAGE>), .. ]
```

# Example
The following example will create a single Pickle output file in the output directory with 2D molecular images
for each input SMILES string from a CSV file.

```
$ python create_images.py --input_file <INPUT_FILE>  --output_dir <OUTPUT_DIRECTORY>
```

# Usage

```
usage: create_images.py [-h] [-d] [-i INPUT_FILE] [-o OUTPUT_DIR]
                        [-bo BAD_OUTPUT_DIR] [-b BATCH_SIZE] [-n NUM_SMILES]
                        [-off OFFSET] [-c CONFIG]

optional arguments:
  -h, --help            show this help message and exit
  -d, --debug           Enables debug logging
  -i INPUT_FILE, --input_file INPUT_FILE
                        input file of smiles to process
  -o OUTPUT_DIR, --output_dir OUTPUT_DIR
                        Output directory. Default : outputs
  -bo BAD_OUTPUT_DIR, --bad_output_dir BAD_OUTPUT_DIR
                        Output directory for bad smile list
  -b BATCH_SIZE, --batch_size BATCH_SIZE
                        Batch size. Default: 0
  -n NUM_SMILES, --num_smiles NUM_SMILES
                        Number of smiles to process (for testing)
  -off OFFSET, --offset OFFSET
                        Offset to start numbering
  -c CONFIG, --config CONFIG
                        Parsl config defining the target compute resource to
                        use. Default: local
```
