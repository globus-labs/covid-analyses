from parsl.app.app import python_app

#may need full paths to work easly with parsl!
#module_dir is location to add to path if needed to specify location of impress_md folder

@python_app
def dock_one_smile(smile, target_oeb_file, license_file, module_dir=None):
    import os
    import sys

    if module_dir is not None:
        sys.path.append(module_dir)

    os.environ['OE_LICENSE'] = license_file
    from openeye import oechem, oedocking
    from impress_md import interface_functions

    docker, _ = interface_functions.get_receptor(target_oeb_file, use_hybrid=True,
                                                        high_resolution=True)
    score, _, _ = interface_functions.RunDocking_(smile,
                                                         dock_obj=docker,
                                                         pos=0,
                                                         name="",
                                                         target_name="",
                                                         force_flipper=True)
    return score

if __name__ == "__main__":
    import os
    import parsl
    from parsl.data_provider.files import File

    from parsl.configs.htex_local import config
    parsl.load(config)

#example smiles - asprin
    smiles_to_dock = ["CC(=O)OC1=CC=CC=C1C(=O)O", ]


    #set current directory as location of license file, impress_md folder
    module_dir = os.getcwd()
    #need your license file full path
    license_file = os.path.join(module_dir,"oe_license.txt")

#example oeb file, you will need one of these
    receptor_files_to_dock = [os.path.join(module_dir,"site1.oeb"),]

    futures = []
    label_list = []
    for receptor_file in receptor_files_to_dock:
        for smile in smiles_to_dock:
            futures.append(dock_one_smile(smile,receptor_file,license_file,module_dir))
            label_list.append("%s - %s"%(smile,receptor_file))


    score_list = [f.result() for f in futures]
    for i in range(len(score_list)):
        print("%s: %f"%(label_list[i],score_list[i]))
