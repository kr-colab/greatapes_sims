import sys
import string
import random
import subprocess
import os

def id_generator(size=10, chars=string.ascii_uppercase + string.digits):
    return ''.join(random.choice(chars) for _ in range(size))

def write_meta(meta_path, var_names, values, rand, jid):
    with open(meta_path+rand+"_jid_"+jid+".info", "w+") as fh:
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

def remove_mutations(ts, start, end, proportion):
    '''
    This function will return a new tree sequence the same as the input,
    but after removing each non-SLiM mutation within regions specified in lists
    start and end with probability `proportion`, independently. So then, if we
    want to add neutral mutations with rate 1.0e-8 within the regions and 0.7e-8
    outside the regions, we could do
      ts = pyslim.load("my.trees")
      first_mut_ts = msprime.mutate(ts, rate=1e-8)
      mut_ts = remove_mutations(first_mut_ts, start, end, 0.3)
    :param float proportion: The proportion of mutations to remove.
    '''
    new_tables = ts.dump_tables()
    new_tables.mutations.clear()
    mutation_map = [-1 for _ in range(ts.num_mutations)]
    for j, mut in enumerate(ts.mutations()):
        keep_mutation = True
        for i in range(len(start)):
            left = start[i]
            right = end[i]
            assert(left < right)
            if i > 0:
                assert(end[i - 1] <= left)
            if mut.position >= left and mut.position < right and len(mut.metadata) == 0:
                keep_mutation = (random.uniform(0, 1) > proportion)
        if keep_mutation:
            mutation_map[j] = new_tables.mutations.num_rows
            if mut.parent < 0:
                new_parent = -1
            else:
                new_parent = mutation_map[mut.parent]
            new_tables.mutations.add_row(site = mut.site, node = mut.node,
                    derived_state = mut.derived_state,
                    parent = new_parent,
                    metadata = mut.metadata)
    return new_tables.tree_sequence()

def overlay_varmut(fpath, neut_mut, gene_info, type="gene"):
    flist = fpath.split('/')
    foname, fname = flist[-2], flist[-1]
    ts_slim = pyslim.load(fpath).simplify()
    ts_mut = msprime.mutate(ts_slim, neut_mut, keep=True)
    print("Mutated", fpath, "in msprime...", flush=True)
    start, end, del_mut = extract_meta(fpath, gene_info, type)
    if del_mut > 0:
        s1=timer()
        ts = remove_mutations(ts_mut, start, end, del_mut/neut_mut)
        s2 = timer()
        print(("Removed extra mutations in genic regions from", fpath, "... Time elapsed (min):"+str(round((s2-s1)/60,3))), flush=True)
    else:
        ts = ts_mut
    s1=timer()
    ts.dump(fpath[:-6]+"_overlaid.trees")
    s2 = timer()
    print(("Dumped overlaid trees to file", fpath, "... Time elapsed (min):"+str(round((s2-s1)/60,3))), flush=True)
