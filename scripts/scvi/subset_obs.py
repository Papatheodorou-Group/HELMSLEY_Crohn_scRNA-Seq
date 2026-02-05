# source ~/.bashrc
# conda activate sklearn-env

from optparse import OptionParser
import os
import sys

###################
# PARSE OPTIONS
###################

usage = "python %prog [options] <in_obj> <metadata> <out_obj>\n"
arg1 = "\nin_obj: input anndata object"
arg2 = "\nmetadata: comma-separated file with obs fields to retain"
arg3 = "\nout_obj: output anndata object"
print("")

parser = OptionParser(usage=usage + arg1 + arg2 + arg3)
(options, args) = parser.parse_args()

if len(args) < 3:
	sys.stderr.write("ERROR: <in_obj>, <metadata> and <out_obj> are required\nTry --help for help\n")
	exit()
	
IN_OBJ = args[0]
METADATA = args[1]
OUT_OBJ = args[2]

###################
# PROCESSING
###################

import scanpy as sc
import pandas as pd
import numpy as np
from scipy import sparse
import pickle

adata = sc.read_h5ad(IN_OBJ)
meta = METADATA.split(",")
adata.obs = adata.obs[meta]
adata.write_h5ad(OUT_OBJ)


exit()

