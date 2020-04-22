import argparse
import os
import parsl
import pickle
import time
from pathlib import Path
from parsl.app.app import python_app

@python_app
def compute_descriptors(smiles=None, smiles_file=None, start_index=0, batch_size=0, out_file=None, bad_file=None, save_csv=False, overwrite=False, save_gzip=False, ignore_3D=False):
    from mordred import Calculator, descriptors
    from rdkit import Chem
    import numpy as np
    import pickle
    import csv
    import os
    import gzip

    if save_gzip:
        raise Exception("GZip not supported yet")

    if not overwrite and  os.path.exists(out_file):
        raise Exception("File exists: %s" % out_file)

    if smiles_file: 
        with open(smiles_file) as current:
            current.seek(start_index)
            smiles = [current.readline() for i in range(batch_size)]

    bad = []
    csv_results = []
    results = {}
    mol = None

    calc = Calculator(descriptors, ignore_3D=ignore_3D)

    for smile_tuple in smiles:
        try:
            s = smile_tuple.split(',')
            smile = s[2].rstrip()
            identifier = s[1].rstrip()
        
            mol = Chem.MolFromSmiles(smile)

            if mol is None or len(mol.GetAtoms()) > 100:
                bad.append((smile, identifier))
            else:
                descs = calc(mol)

                # Mods from Xuefeng
                descs.fill_missing("nan")
                
                descriptor_array = np.array(descs).flatten().astype(np.float32)
                if save_csv: 
                    data = [smile, identifier] + list(descriptor_array)
                    data = ['' if str(i) == 'nan' else i for i in data]
                    csv_results.append(data) 
                else: 
                    results[smile] = ([identifier], np.array(descs).flatten().astype(np.float32))
        except:
            bad.append(s)

    if save_csv: 
        with open(out_file, 'w', newline='') as o_file:
            writer = csv.writer(o_file, delimiter=',', quoting=csv.QUOTE_MINIMAL)
            writer.writerows(csv_results)
        if bad_file and len(bad) > 0:
            with open(bad_file, 'w', newline='') as b_file:
                b_writer = csv.writer(b_file, delimiter=',', quoting=csv.QUOTE_MINIMAL)
                for b in bad:
                    b_writer.writerow(b)
    else:
        with open(out_file, 'wb') as o_file:
            pickle.dump(results, o_file, protocol=pickle.HIGHEST_PROTOCOL)

        if bad_file and len(bad) > 0:
            with open(bad_file, 'wb') as b_file:
                pickle.dump(bad, b_file, protocol=pickle.HIGHEST_PROTOCOL)
         
    return out_file
