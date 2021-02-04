import argparse

parser = argparse.ArgumentParser(description="Syncs the metadata table, the tree sequences and (optionally) the log files. By default, it removes tree sequences and log files whose random identifier is not in the metadata table. You may use the tree sequences as a basis to sync instead by toggling --sync-with-trees")
parser.add_argument("trees_path", help="Path to tree sequences", required=True)
parser.add_argument("meta_table_path", help="Path to metadata table", required=True)
parser.add_argument("--logs_path", help="Path to logs")
parser.add_argument("--trash_path", help="If provided, files are moved to this"
                    "path instead of being deleted.")
parser.add_argument("-t","--sync-with-trees", help="Deletes rows in the metadata table or files in the log path whose tree sequences are not in the trees_path.", action="store_true")
parser.add_argument("-a", "--sync-nontree-files", help="Also removes files ending in .npz or .pops in the trees_path",action="store_true")
parser.parse_args()

meta = pd.read_csv(meta_table_path, sep="\t", header=None, names=["rand", "jid", "info"])

out_files = np.asarray(os.listdir(out_path))
# for each outfile, check if rand is present in meta
out_remove = np.logical_not([any(r in f for r in meta.rand) for f in out_files])
all_rands = np.array([re.findall(r'RAND_(.+).trees', o)[0] for o in out_files])
out_rands = np.unique(all_rands)
if ("o" in in_command):
    for f in out_files[out_remove]:
        shutil.move(out_path+f, trash_path)
        print("moving "+f+" to "+trash_path)

if("m" in in_command):
    meta[np.logical_not(meta.rand.isin(out_rands))].to_csv(trash_path+"removed_rand_jid.txt", sep="\t", header=False, index=False)
    meta[meta.rand.isin(out_rands)].to_csv(meta_table_path, sep="\t", header=False, index=False)
    print("removed",len(meta[np.logical_not(meta.rand.isin(out_rands))]),"rows not in output")

