
from optparse import OptionParser
import os
import sys

###################
# PARSE OPTIONS
###################

usage = "python %prog [options] <in_obj> <out>\n"
arg1 = "\nin_obj: input anndata object"
arg2 = "\nout: output data frame"
arg3 = "\nregion: region code for subsetting (optional)"
print("")

parser = OptionParser(usage=usage + arg1 + arg2)
(options, args) = parser.parse_args()

if len(args) < 2:
	sys.stderr.write("ERROR: <in_obj> and <out> are required\nTry --help for help\n")
	exit()
	
OBJ=args[0]
OUT=args[1]

###################
# EXECUTE
###################

import scanpy as sc
import pandas as pd

adata = sc.read_h5ad(OBJ)
df = adata.obs
if len(args) > 2:
	REGION=args[2]
	df = df[df['Region code'] == REGION]
df.to_csv(OUT, sep = "\t")

exit()

