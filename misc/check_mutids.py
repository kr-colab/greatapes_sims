import sys
import pyslim, tskit
print("1:tspath")

fh=open("bad_trees.txt","a")
ts_path = sys.argv[1]

ts = pyslim.load(ts_path)
dic={}
for mut in ts.mutations():
    mutids = mut.derived_state.split(',')
    for mid in mutids:
        if mid not in dic:
            dic[mid] = mut.position
        elif dic[mid] != mut.position:
            print(ts_path, file=fh)

fh.close()
