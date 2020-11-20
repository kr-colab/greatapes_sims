import sys
import tskit
import os
import glob

print("arg1: path to folder with trees, arg2: filepath to save output")
path = sys.argv[1]
outf = open(sys.argv[2], 'a')

print("Tree filepath\tMax root time", file=outf)
for f in glob.iglob(path+"*[0-9].trees"):
    ts = tskit.load(f)
    print(f"{f}\t{ts.max_root_time}", file=outf)
