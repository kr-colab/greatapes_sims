mport string
import random
import subprocess
import os

def id_generator(size=10, chars=string.ascii_uppercase + string.digits):
    return ''.join(random.choice(chars) for _ in range(size))

def write_meta(meta_path, id, jid):
    with open(meta_path+id+"_jid_"+jid+".info", "w+") as fh:
        print('\t'.join(var_names), file=fh)
        print('\t'.join(values), file=fh)
        print("#################################", file=fh)

def write_slim_sbatch(name, var_names, values, script_path, meta_path, id, time = "4:00:00", mem = "4G"):
    out = "out."+name
    err = "err."+name
    sbatch_name = id+".srun"
    with open(sbatch_name, "w") as fh:
        #SBATCH env variables
        print("#!/bin/bash\n#SBATCH --account=kernlab\n#SBATCH --partition=kern\n#SBATCH --job-name="+name+"\n#SBATCH --time="+time+"\n#SBATCH --mem "+mem+"\n#SBATCH --open-mode=append"+"\n#SBATCH --output="+meta_path+id+".info"+"\n#SBATCH --error="+meta_path+id+".info", file=fh)
        #modules to load on talapas
        print("\nmodule use /projects/apps/shared/modulefiles/\nmodule load python3 tskit SLiM\n", file=fh)
        #slim command
        var_str = ' '.join(["-d "+var_names[i]+"=\\\""+values[i]+"\\\"" for i in range(len(var_names))])
        print("slim -m -t "+var_str+" "+script_path, file=fh)

meta_path = "/home/murillor/projects/greatapes_sims/meta/sims/"
out_path = "/home/murillor/projects/greatapes_sims/output/"

script_path = "/home/murillor/projects/greatapes_sims/scripts/slim/neutral_burnin_greatapes.slim"
name = "neut_burn"
id = id_generator()

var_names = ["path", "prefix", "ancN", "mu", "r", "L", "RAND"]
values = [out_path,"neutral_unmut","10000","1e-9","1.5e-8","50818468",id]

write_slim_sbatch(name, var_names, values, script_path, meta_path, id)
cmd = "sbatch "+id+".srun"
out = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE)

if (out.returncode != 0):
    print("sbatch failed")
    print(out)
else:
    jid = out.stdout.decode().strip().split(" ")[-1]
    write_meta(meta_path, id, jid)

os.remove(id+".srun")
