import parsl
import os
from parsl.app.app import python_app, bash_app
import argparse
import re
import time
from pathlib import Path


@python_app
def compute_neural_fingerprints(smiles=None, smiles_file=None, start_index=0, batch_size=0, out_file=None, bad_file=None, save_csv=False, overwrite=False, save_gzip=False, model_file=None, max_degree=6):
    import os
    import logging
    import pickle
    from rdkit import Chem
    from rdkit.Chem import AllChem
    import csv
    import numpy as np
    import torch
    from NeuralGraph.model import QSAR
    from NeuralGraph.util import dev, enlarge_weights

    if not csv:
        raise Exception("Neural FPs only support CSV output")

    if not overwrite and  os.path.exists(out_file):
        raise Exception("File exists: %s" % out_file)

    #Load the model
    if model_file is not None:
        model_file = Path(model_file)
        if model_file.exists() and model_file.is_file():
            net = torch.load(model_file, map_location=dev)
        else:
            raise FileNotFoundError
    else: # random large weights
        net = QSAR(hid_dim=128, n_class=1, max_degree=6)
        enlarge_weights(net, -1e4, 1e4)

    net = net.to(dev)

    # read smiles
    if smiles_file: 
        with open(smiles_file) as current:
            current.seek(start_index)
            smiles = [current.readline() for i in range(batch_size)]
   
    sep = ","
    bad = []
    results = []
    
    for s in smiles:
        try: 
            mol_tuple = s.split(',')
            dataset = mol_tuple[0].rstrip()
            identifier = mol_tuple[1].rstrip()
            sml = mol_tuple[2].rstrip()
       
            good = True
            mol = Chem.MolFromSmiles(sml)
            if mol:
                atoms = mol.GetAtoms()
                for atom in atoms:
                    if atom.GetDegree() >= max_degree:
                        bad.append(mol_tuple)
                        good = False
            else:
                bad.append(mol_tuple)
                good = False
            if good:
                fp = np.concatenate(net.calc_nfp([sml]))
                #fp = net.calc_nfp([sml])
                fp_ = ':'.join("{:.7f}".format(x) for x in fp)
                results.append((dataset, identifier, sml, fp_))
        except: 
            bad.append(s)

    with open(out_file, 'w') as output_file:
        writer = csv.writer(output_file, delimiter=',', quoting=csv.QUOTE_MINIMAL)
        writer.writerows(results)
    if bad_file and len(bad) > 0:
         with open(bad_file, 'w') as b_file:
             b_writer = csv.writer(b_file, delimiter=',') # quoting=csv.QUOTE_MINIMAL)
             b_writer.writerows(bad)

    return out_file

