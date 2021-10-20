#!/`usr/bin/env python
# coding: utf-8

import allel
import os
import sys
import zarr
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import subprocess
import itertools
from greatapes_func import vcf_to_zarr

meta = pd.read_csv(sys.argv[1], sep="\t", encoding="ISO-8859-1")
print(meta)
print(meta.from_vcf.unique())
species = list(meta.from_vcf.unique())
pairs = list(itertools.combinations(species, 2))

chromosomes = ["chr" + str(i) for i in range(1, 23)]
f=open(sys.argv[4], "a")
for pair in pairs:
    if pair[0] == "other" or pair[1] == "other":
        continue
    vcf_path1 = sys.argv[2] + pair[0] + ".vcf.gz"
    vcf_path2 = sys.argv[2] + pair[1] + ".vcf.gz"
    zarr_path1 = sys.argv[3] + pair[0] + ".zarr"
    zarr_path2 = sys.argv[3] + pair[1] + ".zarr"

    vcf_to_zarr(vcf_path1, zarr_path1)
    vcf_to_zarr(vcf_path2, zarr_path2)

    #print("loading zarr")
    print(pair)
    callset1 = zarr.open_group(zarr_path1, mode="r")
    callset2 = zarr.open_group(zarr_path2, mode="r")
    for c in chromosomes:
        #print(c)
        #print("passou1")
        is_equal = np.equal(callset1[c+"/variants/ALT"][:][np.isin(callset1[c+"/variants/POS"],callset2[c+"/variants/POS"]),0],callset2[c+"/variants/ALT"][:][np.isin(callset2[c+"/variants/POS"],callset1[c+"/variants/POS"]),0])
        pos_triallelic = callset1[c+"/variants/POS"][:][np.isin(callset1[c+"/variants/POS"],callset2[c+"/variants/POS"])][np.logical_not(is_equal)].astype(np.int)
        #assert (pos_triallelic == callset2[c+"/variants/POS"][np.logical_not(is_equal)].astype(np.int)) and ((sum(is_equal)/len(is_equal))>0.9):
        #    print("Something went wrong. Perc of biallelic is:", sum(is_equal)/len(is_equal))
        for loc in pos_triallelic:
            #remember bed files are 0-indexed
            print(c,loc-1,loc, sep='\t', file=f)
        print(pair,c,"perc of biallelic:",sum(is_equal)/len(is_equal))

f.close()
