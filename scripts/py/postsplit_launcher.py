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
mutate=1.66e-8 #rate fromhttps://www.sciencedirect.com/science/article/pii/S0002929715004085?via%3Dihub

params = pd.read_csv(table_path)
params.drop_duplicates(['real_gens', 'real_pop_size_anc', 'real_pop_size1', 'real_pop_size2'], inplace=True) #making sure to not do more sims than necessary
params
gens = tb.gens.unique()

for i in range(len(gens)):
    for j in range(reps):
        rand = id_generator()
        var_names = ["N1","N2", "path_burnin", "gens", "mu", "r", "L", "RAND"]
        values = ["10000","10000", burnin_path, str(gens[i]), "0","1.5e-8","50818468",rand]

        write_sim_sh(var_names, values, prefix, out_path ,script_path, meta_path, rand, mutate, prefix_mut, time="48:00:00")
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
