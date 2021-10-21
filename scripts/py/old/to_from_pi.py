#!/usr/bin/env python
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
from greatapes_func import *

if (len(sys.argv) < 7):
    print("arg1:path to vcf; arg2: path to masking bed file; arg3: window size; arg4: path to output with stats; arg5: path to ancestral allele file; arg6: path to ref fastas")
    sys.exit()

chromosomes = ["chr" + str(i) for i in range(1, 23)]
contig_len = [247249719, 242951149, 199501827, 191273063,
              180857866, 170899992, 158821424, 146274826,
              140273252, 135374737, 134452384, 132349534,
              114142980, 106368585, 100338915, 88827254,
              78774742, 76117153, 63811651, 62435964, 46944323, 49691432]
# win_size
win_size = sys.argv[3]
path_anc=sys.argv[5]
path_fas=sys.argv[6]
print("loading accessibility mask")
dic_acc = get_acc(sys.argv[2], chromosomes, contig_len)

vcf_path = sys.argv[1]
zarr_path = sys.argv[1][:-6]  + "zarr"

vcf_to_zarr(vcf_path, zarr_path, chromosomes)

print("loading zarr")
callset = zarr.open_group(zarr_path, mode="r")

# print(dic_acc)
# callset.tree(expand=True)
tmp = pd.DataFrame()
for c in chromosomes:
    print(c)
    anc_states, seq=anc_states_from_file(c, path_anc, path_fas, node="inner.20")
    positions = callset[c]["variants/POS"][:]-1
    seq[positions] = callset[c]["variants/ALT"][:,0]
    #seq is gonna contain the derived alleles, and anc_states the ancestral
    site_types={"WW":[["A","T"],["A","T"]], "WS":[["A","T"],["G","C"]], "SW":[["G","C"],["A","T"]], "SS":[["G","C"],["G","C"]]}
    for t in site_types:
        nuc = site_types[t]
        print(anc_states, seq)
        print(anc_states.dtype,seq.dtype)
        is_diff = np.char.not_equal(anc_states,seq)
        is_anc = np.isin(anc_states,nuc[0])
        is_der = np.isin(seq,nuc[1])
        is_type = np.logical_and(np.logical_and(is_anc,is_der),is_diff)
        mask = np.logical_and(dic_acc[c],is_type)
        genotypes, pos, n_gt = callset_to_masked_gt(callset, c, mask)
        ac = genotypes.count_alleles()
        print("calculating stats...")
        pi, windows, n_bases, counts = allel.windowed_diversity(pos, ac,
            size=int(win_size), start=1, stop=contig_len[chromosomes.index(c)],
            is_accessible=mask)
        #print("number of accesible bases is:", sum(mask))
        print("pi computed for state:", t)
    #some sanity checks
        tmp = pd.DataFrame()
        tmp["chr"] = np.full(pi.shape,c)
        tmp["start"] = windows[:,0]
        tmp["end"] = windows[:,1]
        tmp["n_acc"] = n_bases
        tmp["pi"] = pi
        tmp["spp"] = np.full(pi.shape,sys.argv[1][:-7])
        tmp["state"] = np.full(pi.shape,t)
        tmp.to_csv(path_or_buf=sys.argv[4]+"pi_site_types_win_" + win_size + ".csv", mode='a', header=False, index=False)

#parallel python ../../scripts/anc_state_pi.py {} ../masks/All.tri.callability_mask.bed 100000 ../../output/ ::: bonobo_anc.vcf.gz bornean_orangutan_anc.vcf.gz central_chimp_anc.vcf.gz eastern_chimp_anc.vcf.gz eastern_gorilla_anc.vcf.gz homo_anc.vcf.gz nigerian_chimp_anc.vcf.gz sumatran_orangutan_anc.vcf.gz western_chimp_anc.vcf.gz western_gorilla_anc.vcf.gz

