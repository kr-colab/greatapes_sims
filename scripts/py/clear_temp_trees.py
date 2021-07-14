"""
Clears temporary tree sequence files related to the specified final tree sequence path.
It assumes the temporary files to be removed have the same file name as the final ts, but with a numeric suffix.
e.g.,
Final tree sequence path: /home/result_sim.trees
Temporary files to be removed: /home/result_sim.trees10, /home/result_sim.trees100
"""

import sys
import glob
import os

# the only argument is the final output file path
outfile = sys.argv[1]

# globbing temp files
temp_files = glob.glob(outfile+"[0-9]*")
for f in temp_files:
    os.remove(f)
