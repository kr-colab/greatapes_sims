
# Script to synchronize files in metadata with outputs (removes outputs not in the meta table), but can also do the opposite.
import pandas as pd
import os
import numpy as np

meta_table_path ="/home/murillor/projects/greatapes_sims/meta/sims/rand_jid.txt"
out_path = "/home/murillor/projects/greatapes_sims/output/"
trash_path = "/home/murillor/projects/greatapes_sims/tmp/"

meta = pd.read_csv(meta_table_path, sep="\t", header=None, names=["rand", "jid", "info"])

out_files = np.asarray(os.listdir(out_path))
# for each outfile, check if rand is present in meta
to_remove = np.logical_not([any(r in f for r in meta.rand) for f in out_files])


