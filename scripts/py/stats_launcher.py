out_path = "/home/murillor/projects/greatapes_sims/output/"
meta_path = "/home/murillor/projects/greatapes_sims/meta/sims/"
script_path = "/home/murillor/projects/greatapes_sims/scripts/py/win_stats_ts.py"
win_size = 1000000
n_pops = 2
prefix = "neut_split_mut"

tree_files = glob(out_path+prefix+"_RAND_*.trees")
for ts_path in tree_files:
    matches = re.match( r'.+RAND_(.+).trees', ts_path )
    rand_id = matches.groups()[0]
    write_sh(out_path, meta_path, script_path, ts_path, win_size, n_pops, prefix, time = "8:00:00", mem = "4G")
    cmd = "sbatch "+rand_id+".sh"
    out = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE)
    if (out.returncode != 0):
        print("sh failed")
        print(out)
    else:
        jid = out.stdout.decode().strip().split(" ")[-1]
        write_meta(meta_path, var_names, values, rand, jid)

    os.remove(rand+".sh")
