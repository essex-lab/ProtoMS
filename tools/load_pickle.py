"""Simple script to load a pickle of calculation results and import
some useful libraries to investigate the data.

To run:
python -i $PROTOMSHOME/tools/load_pickle.py FILENAME

The -i flag given before the script name tells python to drop into an
interactive session after executing the script. The variable 'data' in
the interactive session will contain the loaded object.

This script also be run with ipython for more powerful features:
ipython -i $PROTOMSHOME/tools/load_pickle.py -- FILENAME

Loaded libraries:
pickle
numpy as np
matplotlib.pyplot as plt
all objects from free_energy_base are imported in the default namespace

"""

import pickle
import numpy as np
import os
import sys
import matplotlib
from protomslib import *
from protomslib import simulationobjects
from protomslib.free_energy import *
# from free_energy_base import *
# from calc_ti_decomposed import TI_decomposed
# from gcmc_free_energy_base import GCMCMBAR, GCMCResult, GCMCPMF

if "DISPLAY" not in os.environ or os.environ["DISPLAY"] == "":
    matplotlib.use('Agg')
import matplotlib.pyplot as plt

with open(sys.argv[1], 'rb') as f:
    data = pickle.load(f)
