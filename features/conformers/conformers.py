import parsl
import os
from parsl.app.app import python_app, bash_app
import argparse
import time
from pathlib import Path

# conda install -c openeye openeye-toolkits

@python_app
def compute_conformers(smiles=None, smiles_file=None, start_index=0, batch_size=0, out_file=None, bad_file=None, save_csv=False,  overwrite=False, save_gzip=False, license=None, timeout=0, max_failures=2):
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

    # function to compute enantiomers
    # separated out so we can use an alarm to timeout
    def get_enan(omega, enan):
        enan = oechem.OEMol(enan)
        ret = omega.Build(enan)
        return enan, ret


    def alarm_handler(signum, frame):
        #print("ALARM signal received")
        raise Exception()

    #ofs = oechem.oemolostream()
    #ofs.SetFormat(oechem.OEFormat_OEB)
    #if not ofs.open(out_file):
    #    oechem.OEThrow.Fatal("Unable to open %s for writing" % out_file)

    sep = ","
    bad = []
    mols = []
    for s in smiles:
        mol = s.split(sep)
        if len(mol) > 1: 
            mols.append(mol[2].rstrip())

    # put all mols in a string as we're told it is faster to process this way
    in_smiles = "\n".join(mols)

    ims = oechem.oemolistream()
    ims.SetFormat(oechem.OEFormat_SMI)
    ims.openstring(in_smiles)

    # Turn off logging except errors
    oechem.OEThrow.SetLevel(5)

    filt = oemolprop.OEFilter(oemolprop.OEFilterType_BlockBuster)
    #ofs.open(out_file)
    signal.signal(signal.SIGALRM, alarm_handler)
    oe_results = []
    for mol in ims.GetOEMols():
        if filt(mol):
            oemols = []
            ret_code = None
            omegaOpts = oeomega.OEOmegaOptions(oeomega.OEOmegaSampling_FastROCS)
            omega = oeomega.OEOmega(omegaOpts)
            
            failures = 0
            for enantiomer in oeomega.OEFlipper(mol.GetActive(), 6, True):
                if max_failures > 0 and failures >= max_failures: 
                    break
                if len(oemols) >= 10:
                    break

                ret_code = None
                error = False
                signal.alarm(timeout)
                try: 
                    enantiomer, ret_code = get_enan(omega, enantiomer)
                except:
                    print("Timeout %s" % out_file) 
                    failures += 1
                    error = True
                signal.alarm(0)

                if not error and ret_code == oeomega.OEOmegaReturnCode_Success:
                    halfMol = oechem.OEMol(mol, oechem.OEMCMolType_HalfFloatCartesian)
                    oemols.append(halfMol)
                #else:
                    #oechem.OEThrow.Warning("%s: %s" %
                    #    (enantiomer.GetTitle(), oeomega.OEGetOmegaError(ret_code)))
                oe_results.append(oemols)

    ofs = oechem.oemolostream()
    ofs.SetFormat(oechem.OEFormat_OEB)
    ofs.open(out_file)
    for r in oe_results: 
        for res in r:
            oechem.OEWriteMolecule(ofs, res)

    ofs.close()
 
    return out_file

