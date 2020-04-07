import allel
import subprocess, msprime, pyslim
import matplotlib.pyplot as plt
import numpy as np
import os
import random
import re
from glob import glob
import pickle
import sys
import itertools
from timeit import default_timer as timer
import pandas as pd
from stats_funcs import *

out_path, meta_path, ts_path, win_size, n_pops, prefix = sys.argv[1:]


rand_id, L = get_meta(ts_path, meta_path)

win_size, n_pops, L = int(win_size), int(n_pops), int(L)
tmp = win_stats_from_ts(ts_path, rand_id, n_pops, L, win_size)
tmp['rand_id'] = rand_id
tmp.to_csv(out_path+prefix+"_"+str(win_size)+"_stats.csv", mode = 'a', header=False, index=False)
sys.stdout.flush()


