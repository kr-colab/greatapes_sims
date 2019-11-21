import subprocess
import os
from sim_funcs import *
import numpy as np
import pandas as pd

print("DO NOT FORGET TO LOAD THE RIGHT CONDA ENV")

slurm=True
table_path ="/home/murillor/projects/greatapes_sims/meta/param_sims_spp.csv"
meta_path = "/home/murillor/projects/greatapes_sims/meta/sims/"
out_path = "/home/murillor/projects/greatapes_sims/output/"

script_path = "/home/murillor/projects/greatapes_sims/scripts/slim/neutral_burnin_greatapes.slim"
prefix = "neut_burn_unmut"
#prefix_mut = "neut_burn_mut"
mutate=0

params = pd.read_csv(table_path)
anc_list = params.real_pop_size_anc.unique()
anc_list = np.append(anc_list, 10000)
anc_list = np.array([1000])
for ancN in anc_list:
    rand = id_generator()
    var_names = ["ancN", "mu", "r", "L", "RAND"]
    values = [str(ancN),"0","1.5e-8","200000000",rand]

    write_sim_sh(var_names, values, prefix, out_path ,script_path, meta_path, rand, time = "15-00:00:00", mem = "128G")
    if (slurm):
        cmd = "sbatch "+rand+".sh"
    else:
        cmd = "sh "+rand+".sh"
    out = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE)

    if (out.returncode != 0):
        print("sh failed")
        print(out)
    else:
        jid = out.stdout.decode().strip().split(" ")[-1]
        write_meta(meta_path, var_names, values, rand, jid)

    os.remove(rand+".sh")
