## Using `create_feature.py`

The `create_features.py` program runs a supplied function over many molecules, assumed to be in a CSV file with format `<source>, <id>, <smiles>`.

For example:

```
python create_features.py 
       -t inchi    # generate InChIs, InChiKeys, etc.
       -wr         # Have workers read data, rather than master send to workers (recommended)
       -csv        # Output in CSV format
       -b 1000000  # Block size
       -i /projects/CVD_Research/datasets/release2/canonical_smiles/ENA/ENA.csv # Input
       -o O        # Location for output files
       -bo BO      # Location for bad output files
       -c theta    # Use `theta` configuration
```

The following functions (`-t` values) are currently supported:
* **descriptor**: Computes descriptor
* **druggable**: Determines (True/False) whether molecule passes OpenEye BlockBuster filter 
* **fingerprint**: Computes fingerprint
* **inchi**: Computes InChI, InChiKey, etc. Has timeout.
* **neural_fingerprint**: Details TBD.

## To add another function, say `foo`: 

1) Create directory features/foo
2) Create file `foo/foo.py`, e.g., as a copy of inchis/inchi.py (if you want a timeout)
3) Rename the function `compute_inchi' within your new `foo/foo.py` to `compute_foo`.
4) Rename the function `compute_one_inchi` within your new foo/foo.py to `compute_one_foo`, and update it to perform whatever function you want to apply to each molecule. 
5) Add code to `create_features.py` to load the new module, `foo/foo.py`, and call `compute_foo`. You're done.
