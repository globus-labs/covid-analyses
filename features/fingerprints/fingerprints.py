import parsl
import os
from parsl.app.app import python_app, bash_app
import argparse
import re
import time
from pathlib import Path


@python_app
def compute_fingerprints(smiles=None, smiles_file=None, start_index=0, batch_size=0, out_file=None, bad_file=None, save_csv=False):
    import os
    import logging
    import pickle
    from rdkit import Chem
    from rdkit.Chem import AllChem
    import csv


    if smiles_file: 
        with open(smiles_file) as current:
            current.seek(start_index)
            smiles = [current.readline() for i in range(batch_size)]

    sep = ","
    results = []
    bad = []
    for s in smiles:
        #mol = s.split(sep)i
        mol = None
        try:
            mol = s.split(sep)
            molecule = mol[2].rstrip()
            identifier = mol[1].rstrip()
            dataset = mol[0].rstrip()
            fp = AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(molecule), 2, nBits=2048)
            if save_csv:
                results.append((dataset, identifier, molecule, fp.ToBase64())) 
            else:
                results.append((molecule, identifier, fp))
        except:
             fp = None
             bad.append(mol)

    if save_csv: 
        with open(out_file, 'w') as output_file:
            writer = csv.writer(output_file, delimiter=',') # quoting=csv.QUOTE_MINIMAL)
            writer.writerows(results)
        if bad_file and len(bad) > 0:
            with open(bad_file, 'w') as b_file:
                b_writer = csv.writer(b_file, delimiter=',') # quoting=csv.QUOTE_MINIMAL)
                b_writer.writerows(bad)
    else:
        with open(out_file, 'wb') as output_file:
            pickle.dump(results, output_file, protocol=pickle.HIGHEST_PROTOCOL)
    
        if bad_file and len(bad) > 0:
            with open(bad_file, 'wb') as b_file:
                pickle.dump(bad, b_file, protocol=pickle.HIGHEST_PROTOCOL)

    return out_file

