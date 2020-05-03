import parsl
import os
from parsl.app.app import python_app, bash_app
import argparse
import re
import time
from pathlib import Path


def generate_batch(filename, start=0, batchsize=10, max_batches=10):

    counter = 0
    if max_batches == 0:
        max_batches = 999999999

    with open(filename) as current:

        while current and counter < max_batches:
            yield current.tell()
            counter += 1

            for i in range(batchsize):
                x = current.readline()
                if not x:
                    current = None
                    break

@python_app
def smiles_to_images(smiles=None, smiles_file=None, start_index=0, batch_size=0, out_file=None, bad_file=None, molSize=(128, 128), kekulize=True, mol_name='', mol_computed=True):
    from rdkit import Chem
    from PIL import Image, ImageDraw, ImageFont
    import io
    import cairosvg
    from PIL import Image
    from rdkit.Chem import rdDepictor
    from rdkit.Chem.Draw import rdMolDraw2D
    from rdkit.Chem import RDConfig
    import pickle


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
            molecule = mol[2]
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

        results.append((mol[0], mol[1], mol[2].rstrip(), image))

    with open(out_file, 'wb') as output_file:
        pickle.dump(results, output_file, protocol=pickle.HIGHEST_PROTOCOL)

    if bad_file and len(bad) > 0: 
        with open(bad_file, 'wb') as b_file:
            pickle.dump(bad, b_file, protocol=pickle.HIGHEST_PROTOCOL)
    return out_file

def create_file_names(output_dir, name, start, finish, bad_dir, extension):
     output_file = os.path.join(output_dir, "%s-%s-%s.%s" % (name, start, finish, extension))
     bad_file = None
     if bad_dir:
         bad_file = os.path.join(bad_dir, "%s-%s-%s.%s" % (name, start, finish, extension))
     return (output_file, bad_file)

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--debug", action='store_true',
                        help="Enables debug logging")
    parser.add_argument("-i", "--input_file", default=None,
                        help="input file of smiles to process")
    parser.add_argument("-o", "--output_dir", default="outputs",
                        help="Output directory. Default : outputs")
    parser.add_argument("-bo", "--bad_output_dir", default=None,
                        help="Output directory for bad smile list")
    parser.add_argument("-b", "--batch_size", default=0,
                        help="Batch size. Default: 0")
    parser.add_argument("-n", "--num_smiles", default=0, 
                        help="Number of smiles to process (for testing)")
    parser.add_argument("-off", "--offset", default=0, type=int,
                        help="Offset to start numbering")
    parser.add_argument("-wr", "--worker_read", default=False, action='store_true',
                        help="Read smiles file on worker. Default: False.")
    parser.add_argument("-c", "--config", default="local",
                        help="Parsl config defining the target compute resource to use. Default: local")
    args = parser.parse_args()

 
    if args.config == "local":
        from parsl.configs.htex_local import config
        from parsl.configs.htex_local import config
        config.executors[0].label = "Foo"
        config.executors[0].max_workers = 2
    elif args.config == "theta":
        from configs.theta import config
    elif args.config == "cvd":
        from configs.theta_cvd import config

    parsl.load(config)
    
    output_dir = args.output_dir
    input_file = args.input_file
    bad_dir = args.bad_output_dir
    bad_file = None
    chunksize = int(args.batch_size)
    worker_read = args.worker_read

    if not os.path.isdir(output_dir):
        try:
            os.makedirs(output_dir, exist_ok = True)
            print("Output directory '%s' created" % output_dir)
        except OSError as error:
            print("Output directory '%s' can not be created" % output_dir)
            exit()
    if bad_dir and not os.path.isdir(bad_dir):
        try:
            os.makedirs(bad_dir, exist_ok = True)
            print("bad output directory '%s' created" % bad_dir)
        except OSError as error:
            print("bad output directory '%s' can not be created" % bad_dir)
            exit()

    start = time.time()
    batch_futures = {}
    # Process each file in the dataset
    #for f in os.listdir(input_dir):
        #input_file = os.path.join(input_dir, f)

    f = Path(input_file).stem

    with open(input_file) as ff:
        if int(args.num_smiles) > 0:
            smiles = ff.readlines()[:int(args.num_smiles)]
        else:
            smiles = ff.readlines()

    batch_futures = {}

    if worker_read:
        batch_generator = generate_batch(input_file, start=0, batchsize=int(args.batch_size), max_batches=0)
        i = 0
        for batch_index in batch_generator:
             print("starting batch index: %s %s %s" % (batch_index, i, chunksize))
             output_file, bad_file = create_file_names(output_dir, f, i, i+chunksize, bad_dir, extension)
             batch_futures[(i,i+chunksize)] = create_fingerprints(smiles_file=input_file, start_index=batch_index,
                                                      batch_size=int(args.batch_size), out_file=output_file, bad_file=bad_file, save_csv=save_csv)
             i += chunksize
    else:

        if chunksize == 0 or chunksize > len(smiles):
            print("Processing file: %s" % input_file)
            output_file, bad_file = create_file_names(output_dir, f, 0, len(smiles), bad_dir, 'pkl')
            batch_futures[0] = smiles_to_images(smiles, output_file, bad_file, mol_computed=False)
        else: 
            for i in range(0, len(smiles), chunksize):
                start_num = i + args.offset
                output_file, bad_file = create_file_names(output_dir, f, i, i+chunksize, bad_dir, extension)
                output_file = os.path.join(output_dir, "%s-%s-%s.pkl" % (f, start_num, start_num+chunksize))
                if bad_dir: 
                    bad_file = os.path.join(bad_dir, "%s-%s-%s.pkl" % (f, start_num, start_num+chunksize))
                batch_futures[(i,i+chunksize)] = smiles_to_images(smiles[i:i+chunksize], output_file, bad_file, mol_computed=False)

    
    print("[Main] Waiting for {} futures...".format(len(batch_futures)))
    for i in batch_futures:
        try:
            x = batch_futures[i].result()
            print(f"Chunk {i} is done with {x}")
        except Exception as e:
            print(f"Exception : {e}")
            print(f"Chunk {i} failed")

    total_time = time.time() - start
    print("All done! %s" % total_time)
    

