import parsl
import os
from parsl.app.app import python_app, bash_app
import argparse
import re
import time
from pathlib import Path

@python_app
def compute_images(smiles=None, smiles_file=None, start_index=0, batch_size=0, out_file=None, bad_file=None, save_csv=False,  overwrite=False, save_gzip=False, molSize=(128, 128), kekulize=True, mol_name='', mol_computed=False):
    from rdkit import Chem
    from PIL import Image, ImageDraw, ImageFont
    import io
    import cairosvg
    from PIL import Image
    from rdkit.Chem import rdDepictor
    from rdkit.Chem.Draw import rdMolDraw2D
    from rdkit.Chem import RDConfig
    import pickle
    import os
    import gzip

    if not overwrite and  os.path.exists(out_file):
         raise Exception("File exists: %s" % out_file)

    if smiles_file:
        with open(smiles_file) as current:
            current.seek(start_index)
            smiles = [current.readline() for i in range(batch_size)]
    sep = ","
    results = []
    bad = []
    for s in smiles:
        mol = s.split(sep)
        mc = None
        try:
            molecule = mol[2].rstrip()
            if not mol_computed:
                molecule = Chem.MolFromSmiles(molecule)
            mc = Chem.Mol(molecule.ToBinary())
        except:
            image = None
            bad.append(mol)

        if mc: 
            if kekulize:
                try:
                    Chem.Kekulize(mc)
                except:
                    mc = Chem.Mol(molecule.ToBinary())
            try:
                if not mc.GetNumConformers():
                    rdDepictor.Compute2DCoords(mc)
                drawer = rdMolDraw2D.MolDraw2DSVG(molSize[0], molSize[1])
                drawer.DrawMolecule(mc)
                drawer.FinishDrawing()
                svg = drawer.GetDrawingText()
                image = Image.open(io.BytesIO(cairosvg.svg2png(bytestring=svg, parent_width=100, parent_height=100,
                                                   scale=1)))
                image.convert('RGB')
            except:
                image = None
                bad.append(mol)
        # only save the result if we have an image
        if image: 
            results.append((mol[0], mol[1], mol[2].rstrip(), image))
    
    if save_gzip: 
         with gzip.open(out_file, 'wb') as output_file:
             pickle.dump(results, output_file, protocol=pickle.HIGHEST_PROTOCOL)
         if bad_file and len(bad) > 0:
             with gzip.open(bad_file, 'wb') as b_file:
                 pickle.dump(bad, b_file, protocol=pickle.HIGHEST_PROTOCOL)
    else:
        with open(out_file, 'wb') as output_file:
            pickle.dump(results, output_file, protocol=pickle.HIGHEST_PROTOCOL)

        if bad_file and len(bad) > 0: 
            with open(bad_file, 'wb') as b_file:
                pickle.dump(bad, b_file, protocol=pickle.HIGHEST_PROTOCOL)
    return out_file
