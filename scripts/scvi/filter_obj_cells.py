from optparse import OptionParser
import os
import sys

###################
# PARSE OPTIONS
###################

usage = "python %prog [options] <in_q_obj> <cells> <out_q_obj>\n"
arg1 = "\nin_q_obj: input h5ad object"
arg2 = "\ncells: file with input cells to select"
arg3 = "\nout_q_obj: output h5ad object"
print("")

parser = OptionParser(usage=usage + arg1 + arg2 + arg3)
(options, args) = parser.parse_args()

if len(args) < 3:
	sys.stderr.write("ERROR: <in_q_obj>, <cells> and <out_q_obj> are required\nTry --help for help\n")
	exit()
	
IN_Q_OBJ=args[0]
CELLS=args[1]
OUT_Q_OBJ=args[2]

###################
# EXECUTE
###################

import scanpy as sc
import pandas as pd

adata_query = sc.read_h5ad(IN_Q_OBJ)
cells = pd.read_csv(CELLS, sep = "\t", header = None)
c = cells[0].tolist()
adata_query = adata_query[c,:]
adata_query.write_h5ad(OUT_Q_OBJ)


exit()

