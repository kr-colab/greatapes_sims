import subprocess
import os
import pandas as pd
from sim_funcs import *

print("DO NOT FORGET TO LOAD THE RIGHT CONDA ENVIRONMENT")

slurm=True
table_path ="/home/murillor/projects/greatapes_sims/meta/param_sims_spp.csv"
meta_path = "/home/murillor/projects/greatapes_sims/meta/sims/"
out_path = "/home/murillor/projects/greatapes_sims/output/"
rec_file = "/home/murillor/projects/greatapes_sims/meta/chr12_rec_rate_hg18.tsv"
ex_file = "/home/murillor/projects/greatapes_sims/meta/chr12_exons_hg18.tsv"

sims_table = pd.read_csv("/home/murillor/projects/greatapes_sims/meta/sims/rand_jid.txt", sep="\t", header=None, names=["rand", "jobid", "args"])
reps = 1
out_path = "/home/murillor/projects/greatapes_sims/output/"

script_path = "/home/murillor/projects/greatapes_sims/scripts/slim/neutral_split_greatapes.slim"
burnin_prefix = "neut_burn_unmut"
burnin_path = "/home/murillor/projects/greatapes_sims/output/" #neut_burn_unmut_RAND_bla.trees

prefix = "neut_split"
prefix_mut = "neut_split_mut"
mut_rate=1.66e-8 #rate fromhttps://www.sciencedirect.com/science/article/pii/S0002929715004085?via%3Dihub

params = pd.read_csv(table_path)
params.drop_duplicates(['real_gens', 'real_pop_size_anc', 'real_pop_size1', 'real_pop_size2'], inplace=True) #making sure to not do more sims than necessary
params

for i in range(len(params)):
    #figure out which burnin tree refers to that ancN
    ancN = str(params.loc[i,]["real_pop_size_anc"])
    is_ancN = sims_table.args.str.contains("ancN="+ancN)
    assert sum(is_ancN <= 1), "More than one burnin with this ancN:"+ancN
    burnin_rand = sims_table[is_ancN].rand[0]
    N1=params.loc[i,]["real_pop_size1"]
    N2=params.loc[i,]["real_pop_size2"]
    gens=params.loc[i,]["real_gens"]
    Ntot = N1+N2
    for j in range(reps):
        rand = id_generator()
        var_names = ["N1","N2", "path_burnin", "gens", "mu", "r", "L", "RAND"]
        values = [str(N1),str(N2), burnin_path+burnin_prefix+"_RAND_"+burnin_rand+".trees", str(int(gens)), "0","1.5e-8","200000000",rand]
        write_sim_sh(var_names, values, prefix, out_path ,script_path, meta_path, rand, mut_rate, prefix_mut,time = str(days)+"-00:00:00", mem = str(mem)+"G")
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
