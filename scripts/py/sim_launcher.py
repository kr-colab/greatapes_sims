import subprocess
import os
from sim_funcs import *

slurm=True
meta_path = "/home/murillor/projects/greatapes_sims/meta/sims/"
out_path = "/home/murillor/projects/greatapes_sims/output/"

script_path = "/home/murillor/projects/greatapes_sims/scripts/slim/neutral_burnin_greatapes.slim"
prefix = "neut_burn_unmut"
prefix_mut = "neut_burn_mut"
mutate=1
rand = id_generator()

var_names = ["ancN", "mu", "r", "L", "RAND"]
values = ["10000","0","1.5e-8","50818468",rand]

write_sim_sh(var_names, values, prefix, out_path ,script_path, meta_path, rand, 1e-9, prefix_mut, slurm)
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
