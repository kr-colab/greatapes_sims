import subprocess
import os
from sim_funcs import *

meta_path = "/home/murillor/projects/greatapes_sims/meta/sims/"
out_path = "/home/murillor/projects/greatapes_sims/output/"

script_path = "/home/murillor/projects/greatapes_sims/scripts/slim/neutral_burnin_greatapes.slim"
name = "neut_burn_unmut"
rand = id_generator()

var_names = ["path", "prefix", "ancN", "mu", "r", "L", "RAND"]
values = [out_path,name,"10000","0","1.5e-8","50818468",rand]

write_slim_sbatch(name, var_names, values, script_path, meta_path, rand)
cmd = "sbatch "+rand+".srun"
out = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE)

if (out.returncode != 0):
    print("sbatch failed")
    print(out)
else:
    jid = out.stdout.decode().strip().split(" ")[-1]
    write_meta(meta_path, var_names, values, rand, jid)

os.remove(rand+".srun")
