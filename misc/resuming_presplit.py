import subprocess
import os
import math
from sim_funcs import *
import numpy as np
import pandas as pd

print("DO NOT FORGET TO LOAD THE RIGHT CONDA ENV")

slurm=True
table_path ="/home/murillor/projects/greatapes_sims/meta/param_sims_spp.csv"
meta_path = "/home/murillor/projects/greatapes_sims/meta/sims/"
out_path = "/home/murillor/projects/greatapes_sims/output/"
rec_file = "/home/murillor/projects/greatapes_sims/meta/chr12_rec_rate_hg18.tsv"
ex_file = "/home/murillor/projects/greatapes_sims/meta/chr12_exons_hg18.tsv"

script_path = "/home/murillor/projects/greatapes_sims/scripts/slim/resuming_sel_presplit5N_greatapes.slim"
prefix = "neut_presplit"
#prefix_mut = "neut_burn_mut"
mutate=0

params = pd.read_csv(table_path)
anc_list = params.real_pop_size_anc.unique()
anc_list = np.append(anc_list, [1000, 10000])
anc_list = anc_list[anc_list>50000]
ancN=anc_list[0]
sim_ints=[500]
for simint in sim_ints:
    mem = math.ceil((ancN/600))+2
    days = math.ceil(ancN/30000)
    rand = id_generator()
    var_names = ["ancN", "mu", "recfile", "exonfile", "L", "RAND", "posprop", "poscoef", "delprop", "delcoef","siminterval"]
    values = [str(ancN),"0",rec_file,ex_file,"132000000",rand, "0", "0", "0", "0", str(simint)]
    write_sim_sh(var_names, values, prefix, out_path ,script_path, meta_path, rand, time = str(days)+"-00:00:00", mem = str(mem)+"G")
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
