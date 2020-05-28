import parsl
import os
from parsl.app.app import python_app, bash_app
import argparse
import re
import time
from pathlib import Path

@python_app
def compute_inchis(smiles=None, smiles_file=None, start_index=0, batch_size=0, out_file=None, bad_file=None, overwrite=False, save_gzip=False, save_csv=False, timeout=0 ):
    import os
    import logging
    from rdkit import Chem
    from rdkit.Chem import inchi
    import csv
    import os
    import signal
    import hashlib
    # According to https://github.com/rdkit/rdkit/issues/2683
    from rdkit import RDLogger
    RDLogger.DisableLog('rdApp.*')

    # Take a SMILES as input and return data to be written to output
    def compute_one_inchi(mol):
        molfrom = Chem.MolFromSmiles(mol)
        if molfrom==None:
            return((None, None, None))
        ini = inchi.MolToInchi(molfrom)
        inikey = inchi.InchiToInchiKey(ini)
        hsh = hashlib.md5(ini.encode('utf-8')).hexdigest()
        return((hsh, inikey, ini))

    def alarm_handler(signum, frame):
        print("ALARM signal received")
        raise Exception('alarm')

    if save_gzip or not save_csv: 
        raise Exception("GZip and PKL not supported yet")
 
    if not overwrite and  os.path.exists(out_file):
        raise Exception("File exists: %s" % out_file)

    if smiles_file: 
        with open(smiles_file) as current:
            current.seek(start_index)
            smiles = [current.readline() for i in range(batch_size)]
   
    if len(smiles) == 0: 
        return ""
    
    sep = ","
    results = []
    bad = []
    signal.signal(signal.SIGALRM, alarm_handler)
    for s in smiles:
        s1 = s.rstrip()
        mol = s1.split(sep)
        if len(mol) < 3:
            continue # Could be 'break'?
        molecule = mol[2].rstrip()
        identifier = mol[1].rstrip()
        dataset = mol[0].rstrip()
        error = False
        error_type = 'unknown'
        signal.alarm(timeout)
        smiles_hash = hashlib.md5(molecule.encode('utf-8')).hexdigest()
        try:
            (ini_hash, inikey, ini) = compute_one_inchi(molecule)
        except Exception as e:
            if e.args=='alarm':
                error_type = 'alarm'
            error = True
        signal.alarm(0)

        if error:
            bad.append(mol+[smiles_hash,error_type])
        else:
            results.append((dataset, identifier, smiles_hash, ini_hash, inikey, ini)) 
    
    with open(out_file, 'w') as output_file:
        writer = csv.writer(output_file, delimiter=',')
        writer.writerows(results)
    if bad_file and len(bad) > 0:
        with open(bad_file, 'w') as b_file:
            b_writer = csv.writer(b_file, delimiter=',')
            b_writer.writerows(bad)

    return out_file

