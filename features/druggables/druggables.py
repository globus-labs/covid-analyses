import parsl
import os
from parsl.app.app import python_app, bash_app
import argparse
import time
from pathlib import Path

# conda install -c openeye openeye-toolkits

@python_app
def compute_druggables(smiles=None, smiles_file=None, start_index=0, batch_size=0, out_file=None, bad_file=None, save_csv=False,  overwrite=False, save_gzip=False, license=None, timeout=0, max_failures=2):
    import csv
    import os
    from openeye import oechem
    from openeye import oeomega
    from openeye import oemolprop
    import signal

    os.environ['OE_LICENSE'] = license    
    
    if save_gzip or save_csv: 
        raise Exception("GZip and CSV not supported")
 
    if not overwrite and  os.path.exists(out_file):
        raise Exception("File exists: %s" % out_file)

    if smiles_file: 
        with open(smiles_file) as current:
            current.seek(start_index)
            smiles = [current.readline() for i in range(batch_size)]
    
    if len(smiles) == 0: 
        return ""

    sep = "," 
    mols = []
    full_mols = []
    for s in smiles:
        mol = s.split(sep)
        if len(mol) > 1:
            mols.append(mol[2].rstrip())
            full_mols.append(mol)
    
    # put all mols in a string as we're told it is faster to process this way
    in_smiles = "\n".join(mols)

    ims = oechem.oemolistream()
    ims.SetFormat(oechem.OEFormat_SMI)
    ims.openstring(in_smiles)

    # Turn off logging except errors
    oechem.OEThrow.SetLevel(5)

    filt = oemolprop.OEFilter(oemolprop.OEFilterType_BlockBuster)

    bad = []
    results = []
    for (m, mol) in zip(ims.GetOEMols(), full_mols):
        if len(mol) > 1: 
            results.append((mol[0], mol[1], filt(m)))

    with open(out_file, 'w') as output_file:
        writer = csv.writer(output_file, delimiter=',') # quoting=csv.QUOTE_MINIMAL)
        writer.writerows(results)

    return out_file

