import pandas as pd
import shutil
import os
import glob

base_path = "../../output/"
meta = pd.read_csv("../../output/sims_info.tsv", header=None, sep="\t")
header = pd.read_csv("../../output/header_sims_info.tsv", sep="\t")
meta.columns = header.columns

for rid in meta.rand_id.unique():
    dir_path = f"{base_path}{rid}"
    print(dir_path)
    files = glob.glob(f"{base_path}*{rid}*")
    if not os.path.exists(dir_path):
        os.makedirs("../../output/"+rid)
    for f in files:
        mv = shutil.move(f, dir_path)
        print("Moved: "+mv)
