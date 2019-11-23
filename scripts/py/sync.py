
# Script to synchronize files in metadata with outputs (removes outputs not in the meta table), but can also do the opposite.
import pandas as pd
import os
import numpy as np
import shutil
import sys
import re


if (len(sys.argv) < 2):
    print("You should pass o to sync outputs, m to sync metadata or mo to do both")
    sys.exit()

in_command = sys.argv[1]
meta_table_path ="/home/murillor/projects/greatapes_sims/meta/sims/rand_jid.txt"
out_path = "/home/murillor/projects/greatapes_sims/output/"
trash_path = "/home/murillor/projects/greatapes_sims/tmp/"

meta = pd.read_csv(meta_table_path, sep="\t", header=None, names=["rand", "jid", "info"])

out_files = np.asarray(os.listdir(out_path))
# for each outfile, check if rand is present in meta
out_remove = np.logical_not([any(r in f for r in meta.rand) for f in out_files])
all_rands = np.array([re.findall(r'RAND_(.+).trees', o)[0] for o in out_files])
out_rands = np.unique(all_rands)
if ("o" in in_command):
    for f in out_files[out_remove]:
        shutil.move(out_path+f, trash_path)
        print("moving "+f+" to "+trash_path)

if("m" in in_command):
    meta[np.logical_not(meta.rand.isin(out_rands))].to_csv(trash_path+"removed_rand_jid.txt", sep="\t", header=False, index=False)
    meta[meta.rand.isin(out_rands)].to_csv(meta_table_path, sep="\t", header=False, index=False)
    print("removed",len(meta[np.logical_not(meta.rand.isin(out_rands))]),"rows not in output")

