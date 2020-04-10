# covid-analyses

This repository contains a set of parallel and scalable tools for processing raw molecule data and producing a set of features for use in machine learning pipelines. Specifically, the tools can be used together in a pipeline to convert molecule SMILES into canonicalized SMILES and compute molecular fingerprints, descriptors, and images. 

The tools use the Parsl parallel programming library to enable these processes to be run on multiple cores on a single node or across many nodes concurrently. 


## Installation

The tools are written in Python and have several dependencies that can be installed from Conda. 

The follow instructions outline how to set up a basic Conda environment to run these tools. 

```
conda create --name candle_py3.7 python=3.7

conda install -c rdkit -c mordred-descriptor 
mordred conda install psutil
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

# Example Workflow

The following pipeline can be applied to any column-formatted input dataset. At a minimum we require a column of SMILES, the pipeline can also track an identifier for each molecule.

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
$ python canonicalize.py --input_file ~/data/test/smiles/pubchem.txt --output_dir ~/data/test/pubchem/canned/ --batch_size 1000000
```

## Merge CSV

In many cases, datasets are deposited as a collection of files. For later stages of the workflow it is easier to manage them from a single input file. In this stage we merge many input files into a single CSV. 

```
$ python merge-csv.py --input_dir ~/data/test/pubchem/canned/ --output_file ~/data/test/pubchem/csv/pubchem.csv  --label PC
```

## Compute features

In this part of the pipeline we will compute descriptor, fingerprints, and image features for use with other processing tools.

### Molecular descriptors 

First we will use mordred to compute molecular descriptors. By default this step will create ~1800 descriptors for each molecule.  Note: we set num_smiles here to 10000 to limit computation cost, in production this should be set to 0 to do the entire dataset. We set a batch size to control parallelism. Computation for each batch will be parallelized where possible.

```
$ python create_descriptors.py --input_file ~/data/test/pubchem/csv/pubchem.csv  --output_dir ~/data/test/pubchem/descriptors/  --bad_output_dir ~/data/test/pubchem/descriptors/missing/  --num_smiles 1000 --batch_size 100
```

Before moving on, lets quickly check that the descriptors have been created correctly: 

```python
import pickle
p = pickle.load(open('data/test/pubchem/descriptors/pubchem-0-100.pkl', 'rb'))
print(p)

{'CC(=O)OC(CC(=O)[O-])C[N+](C)(C)C': ([''],
  array([10.106703 , 10.090993 ,  1.       , ..., 62.       ,  8.145833 ,
          3.0416667], dtype=float32)),
 'CC(=O)OC(CC(=O)O)C[N+](C)(C)C': ([''],
  array([10.106703 , 10.090993 ,  1.       , ..., 62.       ,  8.145833 ,
          3.0416667], dtype=float32)),
 'C1=CC(C(C(=C1)C(=O)O)O)O': ([''],
  array([ 8.094414 ,  7.861189 ,  1.       , ..., 59.       ,  5.1944447,
          2.5      ], dtype=float32)),
 'CC(CN)O': ([''],
  array([ 3.0472066,  3.305183 ,  0.       , ..., 14.       ,  3.3611112,
          1.3333334], dtype=float32)),
 'C(C(=O)COP(=O)(O)O)N': ([''],
  array([ 6.9501066,  7.141538 ,  2.       , ..., 41.       ,  5.923611 ,
          2.2916667], dtype=float32)),
```


### Molecular fingerprints

Next we will compute fingerprints for each of the molecules using RDKit to create representative bit vectors. 

```
$ python create_fingerprints.py --input_file ~/data/test/pubchem/csv/pubchem.csv  --output_dir ~/data/test/pubchem/fingerprints/  --bad_output_dir ~/data/test/pubchem/fingerprints/missing/  --num_smiles 1000 --batch_size 100
```

We can now check that the fingerprints have been created

```python
import pickle
p = pickle.load(open('/home/chard/data/test/pubchem/fingerprints/pubchem-0-100.pkl', 'rb'))
print(p[:5])

[('CC(=O)OC(CC(=O)[O-])C[N+](C)(C)C', 
'', 
<rdkit.DataStructs.cDataStructs.ExplicitBitVect object at 0x7f9d3e7a3970>), ('CC(=O)OC(CC(=O)O)C[N+](C)(C)C', 
'', 
<rdkit.DataStructs.cDataStructs.ExplicitBitVect object at 0x7f9d3e7a3c30>), ('C1=CC(C(C(=C1)C(=O)O)O)O', 
'', 
<rdkit.DataStructs.cDataStructs.ExplicitBitVect object at 0x7f9d3e7a38f0>), 
('CC(CN)O', 
'', 
<rdkit.DataStructs.cDataStructs.ExplicitBitVect object at 0x7f9d3e3ca0b0>), 
('C(C(=O)COP(=O)(O)O)N',
 '', 
<rdkit.DataStructs.cDataStructs.ExplicitBitVect object at 0x7f9d3e7a3af0>)]
```

### 2D Molecular images

Finally, we will compute 2D images of each molecule using RDKit.

```
$ python create_images.py --input_file ~/data/test/pubchem/csv/pubchem.csv  --output_dir ~/data/test/pubchem/images/  --bad_output_dir ~/data/test/pubchem/images/missing/  --num_smiles 1000 --batch_size 100
```

Now we will check the images by saving them as PNGs. 

```python
import pickle
p = pickle.load(open(‘~/data/test/pubchem/images/pubchem-0-100.pkll’, 'rb'))

for i in range (0,5):
    p[i][3].save('mol-%s.png' % i)
```
