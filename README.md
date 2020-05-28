# Molecular Feature Extraction Pipeline Tools

This repository contains a set of parallel and scalable tools for processing raw molecule data and producing a set of features for use in machine learning pipelines. Specifically, the tools can be used together in a pipeline to convert the SMILES representation for each molecule to a canonical SMILES to allow for de-duplication and consistency across data sources. Next, for each molecule, three different types of features are computed: 1) molecular fingerprints that encode the structure of molecules; 2) 2D and 3D molecular descriptors; and 3) 2D images of the molecular structure. These features are being used as input to various machine learning and deep learning models that will be used to predict important characteristics of candidate molecules including docking scores, toxicity, and more.

The tools use the Parsl parallel programming library to enable these processes to be run on multiple cores on a single node or across many nodes concurrently. 


## Installation

The tools are written in Python and have several dependencies that can be installed from Conda. 

The follow instructions outline how to set up a basic Conda environment to run these tools. 

```
conda create --name covid_py3.7 python=3.7
```

The tools use Parsl to parallelize execution. To get the latest Parsl install from GitHub

```
pip install git+https://github.com/Parsl/parsl.git
```

If you're using Parsl on Theta you wll need to install psutil seperately. 

```
conda install psutil
pip install parsl
```

Many of the features are computed using mordred, OpenEye, and RDKit.

```
conda install -c rdkit -c mordred-descriptor 
conda install -c openeye openeye-toolkits
```

# Example Workflow

The following pipeline can be applied to any column-formatted input dataset. At a minimum it requires a column of SMILES. The pipeline can also track an identifier for each molecule.

## Preparation
For the following example we’ll work with the PubChem dataset that is available here: https://2019-ncov.e.globus.org/databases/PubChem/smiles.pubchem.txt.gz

To reduce computation time let’s first unzip and take the first 10M molecules. 

```
$ gunzip smiles.pubchem.txt.gz

$ Head -n 10000000 smiles.pubchem.txt >> pubchem.txt
```

## Canonicalize smiles

First we will convert SMILES to a canonical form using Open Babel. To simplify processing with latter stages we will set the batch size to 1M (resulting in 10 files).

```
$ python canonicalize.py --input_file /data/smiles/pubchem.txt --output_dir /data/pubchem/canned/ --batch_size 1000000
```

## Compute features

In this part of the pipeline we will compute descriptor, fingerprints, and image features for use with other processing tools.

### Molecular descriptors 

First we will use mordred to compute molecular descriptors. By default this step will create ~1800 descriptors for each molecule.  Note: we set num_smiles here to 10000 to limit computation cost, in production this should be set to 0 to do the entire dataset. We set a batch size to control parallelism. Computation for each batch will be parallelized where possible.

```
$ python create_features.py --type descriptors --input_file /data/pubchem/csv/pubchem.csv  --output_dir /data/pubchem/descriptors/  --num_smiles 1000 --batch_size 100 --csv
```

### Molecular fingerprints

We compute fingerprints for each of the molecules using RDKit to create representative bit vectors. 

```
$ python create_features.py --type fingerprints --input_file /data/pubchem/csv/pubchem.csv  --output_dir /data/pubchem/fingerprints/  --num_smiles 1000 --batch_size 100 --csv
```

### 2D Molecular images

We compute 2D images of each molecule using RDKit.

```
$ python create_features.py --type images  --input_file /data/pubchem/csv/pubchem.csv  --output_dir /data/pubchem/images/  --num_smiles 1000 --batch_size 100
```

Images are saved in pickle files as PNGs. 

```python
import pickle
p = pickle.load(open(‘/data/pubchem/images/pubchem-0-100.pkl’, 'rb'))

for i in range (0,5):
    p[i][3].save('mol-%s.png' % i)
    
print(p[:5])

[('PC', '', 'CC(=O)OC(CC(=O)[O-])C[N+](C)(C)C', <PIL.PngImagePlugin.PngImageFile image mode=RGB size=128x128 at 0x7F7E1D5259D0>), 
('PC', '', 'CC(=O)OC(CC(=O)O)C[N+](C)(C)C', <PIL.PngImagePlugin.PngImageFile image mode=RGB size=128x128 at 0x7F7E1D4B8810>), 
('PC', '', 'C1=CC(C(C(=C1)C(=O)O)O)O', <PIL.PngImagePlugin.PngImageFile image mode=RGB size=128x128 at 0x7F7E1D058E90>), 
('PC', '', 'CC(CN)O', <PIL.PngImagePlugin.PngImageFile image mode=RGB size=128x128 at 0x7F7E1D058F90>), 
('PC', '', 'C(C(=O)COP(=O)(O)O)N', <PIL.PngImagePlugin.PngImageFile image mode=RGB size=128x128 at 0x7F7E1D05F090>)]
```

### International Chemical Identifier (InChI)

We compute InChIs for each of the molecules using RDKit. 

```
$ python create_features.py --type inchis --input_file /data/pubchem/csv/pubchem.csv  --output_dir /data/pubchem/inchis/  --num_smiles 1000 --batch_size 100 --csv
```

### Conformers

We compute conformers using OpenEye.

```
$ python create_features.py --type conformers --input_file /data/pubchem/csv/pubchem.csv  --output_dir /data/pubchem/conformers/  --num_smiles 1000 --batch_size 100 --license <PATH_TO_OE_LICENSE>
```

### Druggable

We compute druggable molecules using OpenEye's BlockBuster filter.
```
$ python create_features.py --type druggable  --input_file /data/pubchem/csv/pubchem.csv  --output_dir /data/pubchem/conformers/  --num_smiles 1000 --batch_size 100--csv --license <PATH_TO_OE_LICENSE>
```

