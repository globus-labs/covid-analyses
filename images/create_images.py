import parsl
import os
from parsl.app.app import python_app, bash_app
import argparse
import re

@python_app
def smiles_to_image(in_file, out_file, start=None, end=None, molSize=(128, 128), kekulize=True, mol_name='', mol_computed=True):
    from rdkit import Chem
    from PIL import Image, ImageDraw, ImageFont
    import io
    import cairosvg
    from PIL import Image
    from rdkit.Chem import rdDepictor
    from rdkit.Chem.Draw import rdMolDraw2D
    from rdkit.Chem import RDConfig
    import pickle

    p = pickle.load(open(in_file, "rb"))

    if start is None:
        start = 0
    if end is None:
        end = len(p)
    smiles = p[start:end]


    results = []
    for mol_tuple in smiles:
        try:
            mol = mol_tuple[0]
            if not mol_computed:
                mol = Chem.MolFromSmiles(mol)
            mc = Chem.Mol(mol.ToBinary())
        except:
            image = None
        if kekulize:
            try:
                Chem.Kekulize(mc)
            except:
                mc = Chem.Mol(mol.ToBinary())
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
        results.append((mol_tuple[0], mol_tuple[1], image))
    with open(out_file, 'wb') as output_file:
        pickle.dump(results, output_file, protocol=pickle.HIGHEST_PROTOCOL)
    return out_file


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--debug", action='store_true',
                        help="Enables debug logging")
    parser.add_argument("-n", "--input_dir", default=None,
                        help="input directory of smiles to process")
    parser.add_argument("-b", "--batch_size", default=0, type=int,
                        help="Size of the batch of smiles to send to each node for processing. 0 means follow batch size of existing files. Default : 0")
    parser.add_argument("-o", "--output_dir", default="outputs",
                        help="Output directory. Default : outputs")
    parser.add_argument("-c", "--config", default="local",
                        help="Parsl config defining the target compute resource to use. Default: local")
    args = parser.parse_args()

 
    if args.config == "local":
        from parsl.configs.htex_local import config
        from parsl.configs.htex_local import config
        config.executors[0].label = "Foo"
        config.executors[0].max_workers = 2
    elif args.config == "theta":
        from theta import config
    elif args.config == "cvd":
        from theta_cvd import config

    parsl.load(config)
    
    batch_size = args.batch_size
    output_dir = args.output_dir
    input_dir = args.input_dir

    if not os.path.isdir(output_dir):
        try:
            os.makedirs(output_dir, exist_ok = True)
            print("Output directory '%s' created" % output_dir)
        except OSError as error:
            print("Output directory '%s' can not be created" % output_dir)
            exit()

    results = []

    for f in os.listdir(input_dir):
        input_file = os.path.join(input_dir, f)
        if args.batch_size > 0:
            try:
                result = re.findall('([\d]+-[\d]+)', f)
                file_range = result[-1].split('-')
                assert(len(file_range) == 2)
                span = (int(file_range[0]), int(file_range[1]))
                span_difference = span[1] - span[0]
                pos = re.search('([\d]+-[\d]+)', f)
                base_name = f[:pos.start(0)]
            except:
                print("Can't find range in file: %s" % f)
        print("Starting %s" % f)

        # define ranges for batch-sized chunks
        #print("%s %s" % (batch_size, span_difference))
        if batch_size == 0 or batch_size > span_difference:
             print("Processing entire file: %s" % input_file)
             output_file = os.path.join(args.output_dir, f)
             results.append(smiles_to_image(input_file, output_file, mol_computed = False))
        else:
            print("Processing %s batches: %s" % (int(span_difference/batch_size), input_file))
            for i in range (0, int(span_difference/batch_size)):
                start = span[0] + (i * batch_size)
                end = span[0] + ((i +1) * batch_size)
                output_file = os.path.join(output_dir, '%s%s_%s.pkl' % (base_name, start, end))
                results.append(smiles_to_image(input_file, output_file, i*batch_size, (i+1)*batch_size,mol_computed = False))
            if span_difference % batch_size > 0:
                pass
                results.append(smiles_to_image(input_file, output_file, (i+1)*batch_size, mol_computed = False)) 
    for r in results:
        print(r.result())



    #Identify the files we want to process
    #files = []
    #for i in range (50,100):
    #    if i == 0: 
    #        files.append((0, i+1 * 1000000))
    #    else:
    #        files.append((i*1000000, (i+1)*1000000))
    #file_names = []
    #for f in files:
    #    file_names.append((f[0], "/projects/CVD_Research/datasets/Fingerprints-64/pubchem_canonical.smi/pubchem_canonical.chunk-%s-%s.pkl" % (f[0], f[1])))

    #results = []
    # Go through the input files and chunk into 10K files

    # import time

    #start = time.time()

    #for f in file_names:
    #    print("Starting %s" % f[1])
    #    input_file = f[1]
    
        # define ranges for 10K chunks
    #    chunks = []
    #    for i in range (0,100):
    #        if i == 0: 
    #            chunks.append((0, (i+1)* 10000))
    #        else:
    #            chunks.append((i*10000+1, (i+1)*10000))
#
#        for c in chunks:
#            output_file = '/projects/CVD_Research/datasets/images/pubchem_canonical.smi/pubchem_images_%s_%s.pkl' % (c[0]+f[0], c[1]+f[0])
#            results.append(smiles_to_image(input_file, output_file, c[0], c[1],mol_computed = False))
#
#    for r in results: 
#        print(r.result())
#
#    total_time = time.time() - start
#
#    print("Completed %s from %s files in %s" % (len(results), len(file_names), total_time)) 



