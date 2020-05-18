import sys
import glob
import numpy as np

# the only argument is the final output file path
outfile = sys.argv[1]

# globbing temp files
temp_files = glob.glob(outfile+"[0-9]*")
if len(temp_files) > 0:
    gens = np.array([int(t[len(outfile):]) for t in temp_files])
    # outputting the temp file with max gens
    print(temp_files[np.argmax(gens)])
    print(np.max(gens))
