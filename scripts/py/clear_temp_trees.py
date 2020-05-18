import sys
import glob
import os

# the only argument is the final output file path
outfile = sys.argv[1]

# globbing temp files
temp_files = glob.glob(outfile+"[0-9]*")
for f in temp_files:
    os.remove(f)
