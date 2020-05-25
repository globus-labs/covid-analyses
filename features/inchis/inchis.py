import parsl
import os
from parsl.app.app import python_app, bash_app
import argparse
import re
import time
from pathlib import Path


@python_app
def compute_inchis(smiles=None, smiles_file=None, start_index=0, batch_size=0, out_file=None, bad_file=None, save_csv=False,  overwrite=False, save_gzip=False):
    import os
    import logging
    import pickle
    from rdkit import Chem
    from rdkit.Chem import AllChem, inchi
    import csv
    import os
    import gzip
    import signal
    import hashlib

    def compute_inchi(f, mol):
        ini = inchi.MolToInchi(Chem.MolFromSmiles(mol))
        #f.write(' inchi %s\n'%ini)
        #f.flush()
        inikey = inchi.InchiToInchiKey(ini)
        #f.write(' inchikey %s\n'%inikey)
        #f.flush()
        hsh = hashlib.md5(ini.encode('utf-8')).hexdigest()
        return((hsh, inikey, ini))

    def alarm_handler(signum, frame):
        print("ALARM signal received")
        raise Exception('alarm')

    #f = open("XXX.txt", "w")

    #f.write('COMPUTE_INCHIS XXX\n')
    #f.flush()
    #f.write('FILE1 %s\n'%smiles_file)
    #f.flush()

    if save_gzip: 
        raise Exception("GZip not supported yet")
 
    if not overwrite and  os.path.exists(out_file):
        raise Exception("File exists: %s" % out_file)

    if smiles_file: 
        #f.write('FILE %s\n'%smiles_file)
        #f.flush()
        with open(smiles_file) as current:
            current.seek(start_index)
            smiles = [current.readline() for i in range(batch_size)]
            #f.write('SMILES %s'%smiles)
            #f.flush()
   
    #f.write('SMILES %d\n'%len(smiles))
    #f.flush()

    if len(smiles) == 0: 
        return ""
    
    timeout = 120
    sep = ","
    results = []
    bad = []
    #cnt = 0
    signal.signal(signal.SIGALRM, alarm_handler)
    for s in smiles:
        s = s.rstrip()
        #f.write('\nNext %d: %s\n'%(cnt,s))
        #cnt += 1
        #f.flush()
        mol = s.split(sep)
        if len(mol) < 3:
            bad.append(mol+['length_error'])
            continue
        molecule = mol[2].rstrip()
        identifier = mol[1].rstrip()
        dataset = mol[0].rstrip()
        #f.write(' mol = %s\n'%molecule)
        #f.flush()
        error = False
        error_type = 'unknown'
        signal.alarm(timeout)
        smiles_hash = hashlib.md5(molecule.encode('utf-8')).hexdigest()
        #f.write(' smiles_hash = %s\n'%smiles_hash)
        #f.flush()
        try:
            (ini_hash, inikey, ini) = compute_inchi(f, molecule)
            #f.write(' success')
            #f.flush()
        except Exception as e:
            #f.write(' error: %s'%e.args)
            #f.flush()
            if e.args=='alarm':
                error_type = 'alarm'
            error = True
        signal.alarm(0)

        if error:
            #f.write(' bad XXXXXX')
            #f.flush()
            bad.append(mol+[smiles_hash,error_type])
        elif save_csv:
            results.append((dataset, identifier, smiles_hash, ini_hash, inikey, ini)) 
        else:
            results.append((identifier, smiles_hash, ini_hash, inikey, ini))
    
    if save_csv: 
        with open(out_file, 'w') as output_file:
            writer = csv.writer(output_file, delimiter=',')
            writer.writerows(results)
        if bad_file and len(bad) > 0:
            with open(bad_file, 'w') as b_file:
                b_writer = csv.writer(b_file, delimiter=',')
                b_writer.writerows(bad)
    else:
        with open(out_file, 'wb') as output_file:
            pickle.dump(results, output_file, protocol=pickle.HIGHEST_PROTOCOL)
    
        if bad_file and len(bad) > 0:
            with open(bad_file, 'wb') as b_file:
                pickle.dump(bad, b_file, protocol=pickle.HIGHEST_PROTOCOL)

    return out_file

