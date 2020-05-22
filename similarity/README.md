# Similarity search tool

The program `find_similar.py` is a tool for using molecular fingerprints to find similar molecules. It uses the Parsl parallel programming library to enable these processes to be run on multiple cores on a single node or across many nodes concurrently. 

See [the COVID Analyses tools](https://github.com/globus-labs/covid-analyses/README.md) for installation instructions.

## Description

A **molecular fingerprint** is a mapping of a molecule to a bit vector (in our case, 2048 bits). Given fingerprints for two molecules, one can then use the [Tanimoto distance](https://en.wikipedia.org/wiki/Jaccard_index) between those two fingerprints to obtain an estimate of the [chemical similarity](https://en.wikipedia.org/wiki/Chemical_similarity) between the molecules. 

The basic idea of this tool is to take a set of *NT* **target molecules** and *NS* **source molecule fingerprints** contained in *NF* files, and:
* compute a fingerprint for each target molecule
* compute the similarity between each target molecule and the source molecules in the *NF* files (a total of *NTxNS* similarities)
* for each target, select the *N* source molecules from each of the *NF* files with the highest similarity scores
* write out the top-scoring molecules: a total of *NTxNFxN* molecules 

The primary inputs to the program are:
1. a set of one or more files containing target molecules, specified as SMILES strings (the parameter ``--target_glob``);
1. a set of one or more files containing fingerprints for source molecules (the parameter `--smile_glob`);
1. the value *N* (the parameter `--top_n_matches`)

The outputs are a set of files, one each per XXX, containing lines of the form:
```
XXXX
```

Notes:
* The fact that the program computes top *N* matches separately for each source file is anachronistic: it should return just the top *N* matches overall.
* The *NTxNFxN* matches output can contain duplicates if two or more targets have high similarities to the same source molecule, and/or if the same molecule appears in two or more source files. 
* We have computed, and distribute via the [nCov-Group Data Repository](https://2019-ncovgroup.github.io/data/), fingerprints for close to 4 billion small molecules. The code used to compute those fingerprints is in the [the COVID Analyses tools](https://github.com/globus-labs/covid-analyses/README.md).
* The tool also takes an argument `--top_n_targets` to specify the number of targets to take from each supplied target file.

## Running the tool

To see details on arguments, run:
```
python find_similar.py --help
```

Here is an example run:
```
python find_similar.py --smile_glob "/myhome/fingerprints/*" --target_glob "/myhome/targets/*" -o "/myhome/outputs" -n FOO
```


## A further note
* The approach taken here of computing *NTxNS* similarities is precise but inefficient. We might want instead to create an index to enable faster, if less precise, searches. For an example of how to do this, see [this paper](http://www.davidanastasiu.net/pdf/papers/2016-AnastasiuK-DSAA-tapnn.pdf.).
