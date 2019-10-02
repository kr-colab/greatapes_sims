import sys
import string
import random
import subprocess
import os

def id_generator(size=10, chars=string.ascii_uppercase + string.digits):
    return ''.join(random.choice(chars) for _ in range(size))

def write_meta(meta_path, var_names, values, rand, jid):
    with open(meta_path+rand+".meta", "w+") as fh:
        print('\t'.join(var_names), file=fh)
        print('\t'.join(values), file=fh)
    with open(meta_path+"rand_jid.txt", "a") as fh:
        both = [str(n)+"="+str(v) for n,v in zip(var_names, values)]
        print(rand+"\t"+jid+"\t"+"\t".join(both), file=fh)

def write_sim_sh(var_names, values, prefix, out_path ,script_path, meta_path, rand, mutate=0, prefix_mut="", slurm=True, time = "48:00:00", mem = "16G"):
    outfile = out_path+prefix+"_RAND_"+values[var_names.index("RAND")]+".trees"
    out_mut = out_path+prefix_mut+"_RAND_"+values[var_names.index("RAND")]+".trees"
    var_names.append("outfile")
    values.append(outfile)
    sh_name = rand+".sh"
    with open(sh_name, "w") as fh:
        print("#!/bin/bash", file=fh)
        if (slurm):
            #SBATCH env variables
            print("#SBATCH --account=kernlab\n#SBATCH --partition=kern\n#SBATCH --job-name="+prefix+"\n#SBATCH --time="+time+"\n#SBATCH --mem "+mem+"\n#SBATCH --open-mode=append"+"\n#SBATCH --output="+meta_path+rand+".info"+"\n#SBATCH --error="+meta_path+rand+".info", file=fh)
            #modules to load on talapas
            print("\nmodule use /projects/apps/shared/modulefiles/\nmodule load python3 tskit SLiM\n", file=fh)
        #slim command
        var_str = ' '.join(["-d "+var_names[i]+"=\\\""+values[i]+"\\\"" for i in range(len(var_names))])
        print("slim -m -t "+var_str+" "+script_path, file=fh)
        if(mutate):
            print("python overlay.py "+outfile+" "+out_mut+" "+str(mutate), file=fh)


