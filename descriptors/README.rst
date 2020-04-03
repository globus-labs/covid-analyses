Installation
------------


* Set up a conda env ::
    conda create --name candle_py3.7 python=3.7

* Install requirements ::
    conda install -c rdkit -c mordred-descriptor mordred

    conda install psutil
* Install parsl from master:
   
    pip install git+https://github.com/Parsl/parsl.git
    

Installing parsl on Theta
-------------------------

On most systems install parsl via pip works::
    pip install parsl

However on some machines, like Theta, the pip module for `psutil` simply won't install properly.
In this situation, it's easier to install the `psutil` module from conda and then do the parsl install::
    conda install psutil
    pip install parsl


Running the workflow
--------------------

To run on theta, try::
     python3 test.py --num_smiles 10000 --batch_size 1000 -c theta -s <SMILE_FILE> -i <SMILE_OUTPUT_DIR>

For eg::

     python3 test.py -s 2019q3-4_Enamine_REAL_02.smi -b 10000 -n 0 -o 2019q3-4_Enamine_REAL_02_descriptors -c theta

Note::
  Expect the command to take a while to start running since loading the whole mass of python libs necessary from
  the shared-fs will be slow. But, once loaded, it should fly


Usage::
    usage: test.py [-h] [-v VERSION] [-d] [-n NUM_SMILES] [-s SMILE_FILE]
    [-b BATCH_SIZE] [-c CONFIG]

    optional arguments:
    -h, --help            show this help message and exit
    -v VERSION, --version VERSION
    Print Endpoint version information
    -d, --debug           Enables debug logging
    -n NUM_SMILES, --num_smiles NUM_SMILES
    Number of smiles to load and run. Default=10000, if
    set to 0 the entire file will be used
    -s SMILE_FILE, --smile_file SMILE_FILE
    File path to the smiles csv file
    -b BATCH_SIZE, --batch_size BATCH_SIZE
    Size of the batch of smiles to send to each node for
    processing. Default=4, should be 10K
    -c CONFIG, --config CONFIG
    Parsl config defining the target compute resource to
    use. Default: local


Notes:
------

Following are some of the changes to the original code in `main.py` that needed some tweaking:

The `ModelInferer` class has basically been stripped to avoid a deserialization error. Instead we just ship out the
model itself, on which we call `model.infer(...)`.

We do explicit pickle based serialization/deserialization to avoid a bug in the ipyparallel.serialize library that
Parsl uses underneath. I'm working on a cleaner fix.



Notebooks :

`initial_timing_study.ipynb` looks at the cost of serializing and dispatching tasks to a remote worker.

`resilient_workflow.ipynb` looks at creating batches of tasks and handling failures.

Notes::
  Train contains some smiles to play aroudn with. The models are fake here
  since loading up the real models would take some time to do.

  This interface is exactly the same however.
