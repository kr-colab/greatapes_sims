import subprocess
import os
import pandas as pd
from sim_funcs import *

print("DO NOT FORGET TO LOAD THE RIGHT CONDA ENVIRONMENT")

slurm=True
table_path ="/home/murillor/projects/greatapes_sims/meta/param_sims_spp.csv"
meta_path = "/home/murillor/projects/greatapes_sims/meta/sims/"
reps = 10
out_path = "/home/murillor/projects/greatapes_sims/output/"

script_path = "/home/murillor/projects/greatapes_sims/scripts/slim/neutral_split_greatapes.slim"
burnin_path = "/home/murillor/projects/greatapes_sims/output/neut_burn_unmut_RAND_R3MVZ1633Q.trees"
prefix = "neut_split"
prefix_mut = "neut_split_mut"
mutate=1e-9

tb = pd.read_csv(table_path)
gens = tb.gens.unique()

for i in range(len(gens)):
    for j in range(reps):
        rand = id_generator()
        var_names = ["N1","N2", "path_burnin", "gens", "mu", "r", "L", "RAND"]
        values = ["10000","10000", burnin_path, str(gens[i]), "0","1.5e-8","50818468",rand]

        write_sim_sh(var_names, values, prefix, out_path ,script_path, meta_path, rand, mutate, prefix_mut)
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
