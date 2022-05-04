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
from Bio import SeqIO
import glob

def vcf_to_zarr(vcfs, base_names, zarr_path, chromosomes,force_zarr=False, vcf_suffix=".vcf.gz"):
    for vcf, bname in zip(vcfs, base_names):
        bname = vcf.split("/")[-1][0:-len(vcf_suffix)]
        zarr_out = f"{zarr_path}{bname}"
        is_zarr = os.path.exists(zarr_out)
        if (not is_zarr) or force_zarr:
            print(f"Converting {vcf} to zarr.")
            for chrom in chromosomes:
                print(f"\t{chrom}")
                allel.vcf_to_zarr(vcf, zarr_out, group=chrom, region=chrom, fields="*", log=sys.stdout, overwrite=True)
        else:
            print(f"{vcf} had already been converted to zarr.")

def add_acc_group(zarr_path, bed_path, chroms, spp_to_sbp):
    for spp, sbps in spp_to_sbp.items():
        print(f"Adding accessibility mask {spp} to {sbps}")
        zarrs = [zarr.open_group(f'{zarr_path}{sbp}', mode='r+') for sbp in sbps]
        acc = None
        c = None
        cols = [None]
        if all([f'{c}/accessibility' in z for z in zarrs for c in chroms.chr]):
            print(f"Accessibility mask already there for {sbps}")
            continue
        with open(f'{bed_path}{spp}.callability_mask.bed') as fp:
            for line in fp:
                cols = line.strip().split("\t")
                # started a new chrom
                if c != cols[0]:
                    print(f"\t{cols[0]}")
                    # if this is one of the chroms i nthe zarr
                    if c in zarrs[0]:
                        # if not the first line, time to save to file
                        if c is not None:
                            print(f'Savint {c} to zarr.')
                            for z in zarrs:
                                z[c].create_dataset(name="accessibility", shape=acc.shape, data=acc)
                    if cols[0] in zarrs[0]:
                        acc = np.full(shape=(chroms[chroms.chr==cols[0]].clen.iloc[0]), fill_value=True, dtype='bool')
                    c = cols[0]
                if cols[0] in zarrs[0]:
                    acc[int(cols[1]):int(cols[2])] = np.full(int(cols[2]) - int(cols[1]), False)
        if (c is not None ) and (c in zarrs[0]):
            print(f'Savint {c} to zarr.')
            for z in zarrs:
                z[c].create_dataset(name="accessibility", shape=acc.shape, data=acc)

def callset_to_masked_gt(callset, c, acc):
    '''This function takes a callset, a chromosome and an accessibility mask and
    returns a tuple with the GenotypeArray, an array with the positions, the alternate alleles (assuming sites are biallelic within the VCF) and the
    number of genotyped individuals'''
    vcf_features = {}
    call_shape = callset[c + '/calldata/GT'].shape
    pos = allel.SortedIndex(callset[c + '/variants/POS'])
    # oindex is a zarr function that indexes along a specified dimension
    # removing inacessible rows
    gt = callset[c + '/calldata/GT'].oindex[acc[pos-1], :, :]
    vcf_features["mean_dp"] = np.nanmean(callset[c + '/calldata/DP'].oindex[acc[pos-1], :], axis=1)
    vcf_features["mean_gq"] = np.nanmean(callset[c + '/calldata/GQ'].oindex[acc[pos-1], :], axis=1)
    vcf_features["mq"] = callset[c + '/variants/MQ'].oindex[acc[pos-1]]
    vcf_features["qual"] = callset[c + '/variants/QUAL'].oindex[acc[pos-1]]
    vcf_features["qd"] = callset[c + '/variants/QD'].oindex[acc[pos-1]]
    alt = callset[c + '/variants/ALT'].oindex[acc[pos-1], 0]
    ref = callset[c + '/variants/REF'].oindex[acc[pos-1]]
    pos = pos[acc[pos-1]]
    return(allel.GenotypeArray(gt), pos, ref, alt, call_shape[1], vcf_features)

def merge_allel_arrays(pos1, pos2, ac1, ac2, n1, n2, ref1, ref2, alt1=None, alt2=None):
    """
    Merges the positions arrays 1 and 2, and remaps the allele counts arrays to the merged positions. If arrays with alternative alleles are provided, then the merge checks for and removes triallelic sites, otherwise it assumes the alternate allele in ac1 is the same as in ac2.
    """
    p = np.unique(np.concatenate((pos1,pos2)))
    is_pos1 = np.isin(p, pos1)
    is_pos2 = np.isin(p, pos2)
    # removing any potential triallelic sites
    merged_alt1 = np.full(len(p), '', dtype=object)
    merged_alt2 = merged_alt1.copy()
    if (alt1 is not None) and (alt2 is not None):
        merged_alt1[is_pos1] = alt1
        merged_alt2[is_pos2] = alt2
    is_nottri = np.logical_or(merged_alt1 == merged_alt2, np.logical_or(merged_alt1 == '', merged_alt2 == ''))
    merged_ref1 = np.full(len(p), '', dtype=object)
    merged_ref2 = merged_ref1.copy()
    merged_ref1[is_pos2] = ref2
    merged_ref1[is_pos1] = ref1
    merged_ref2[is_pos1] = ref1
    merged_ref2[is_pos2] = ref2
    assert np.all(merged_ref1==merged_ref2)
    assert np.all(merged_ref1!='')
    assert np.all(merged_ref2!='')
    del merged_ref2
    #now we create new count arrays
    merged_ac1 = np.zeros((len(p),2), dtype=np.int32) #all sites are biallelic
    merged_ac2 = merged_ac1.copy()
    #assume everyone is mono at first
    merged_ac1[:,0] = int(2*n1)
    merged_ac2[:,0] = int(2*n2)
    #for the sites originally present, replace with the observed allele counts
    merged_ac1[is_pos1,:] = ac1
    merged_ac2[is_pos2,:] = ac2
    # filtering out triallelics
    p = p[is_nottri]
    merged_ac1 = merged_ac1[is_nottri]
    merged_alt1 = merged_alt1[is_nottri]
    merged_ac2 = merged_ac2[is_nottri]
    merged_alt2 = merged_alt2[is_nottri]
    merged_ref1 = merged_ref1[is_nottri]
    # preparing array with alleles
    alleles1 = np.full(p.shape, '', dtype=object)
    alleles2 = alleles1.copy()
    # fixed for the ref
    is_fixed_ref1 = merged_ac1[:,1] == 0
    is_fixed_ref2 = merged_ac2[:,1] == 0
    alleles1[is_fixed_ref1] = merged_ref1[is_fixed_ref1]
    alleles2[is_fixed_ref2] = merged_ref1[is_fixed_ref2]
    # fixed for the alt
    is_fixed_alt1 = merged_ac1[:,0] == 0
    is_fixed_alt2 = merged_ac2[:,0] == 0
    alleles1[is_fixed_alt1] = merged_alt1[is_fixed_alt1]
    alleles2[is_fixed_alt2] = merged_alt2[is_fixed_alt2]
    # else it's polymorphic
    is_poly1 = ~np.logical_or(is_fixed_ref1, is_fixed_alt1)
    is_poly2 = ~np.logical_or(is_fixed_ref2, is_fixed_alt2)
    alleles1[is_poly1] = merged_ref1[is_poly1]+'/'+merged_alt1[is_poly1]
    alleles2[is_poly2] = merged_ref1[is_poly2]+'/'+merged_alt2[is_poly2]
    assert np.all(alleles1!='/')
    assert np.all(alleles2!='/')
    #import pdb; pdb.set_trace()
    return (p, merged_ac1, merged_ac2, alleles1, alleles2)

def anc_dict_from_file(c, path_anc, path_fas, node="inner.20"):
    '''This function takes a chrommosome c, the path to the ancestral calls,  and the
    path to the reference fastas, and returns a dictionary with states in the ancestor
    as keys, and boolean np arrays indicating which positions were that state'''
    anc_state = dict()
    states=["A", "T", "C", "G"]
    anc = pd.read_csv(path_anc+c+".calls.txt.gz", header=None, skiprows=1, sep="\t",
            names = ["chromosome","position","leaf.1","leaf.2","leaf.3","leaf.4","leaf.5","leaf.6","leaf.7","leaf.8","leaf.9","leaf.10","leaf.11","inner.12","inner.13","inner.14","inner.15","inner.16","inner.17","inner.18","inner.19","inner.20","hg18","rheMac2"])
    anc["position"] = anc["position"]-1
    anc_alleles = anc[node]
    #anc_alleles = callset[c+'/variants/ANC_ALLELE'][:]
    #anc_pos = callset[c+ '/variants/POS'][:]
    anc_pos = anc["position"]
    record = SeqIO.read(path_fas+c+".fa.masked", "fasta")
    seq = np.array(list(str(record.seq)))
    for state in states:
        anc_state[state] = (seq == state)
        is_state = (anc_alleles==state)
        anc_state[state][anc_pos[is_state]] = True
        anc_state[state][anc_pos[np.invert(is_state)]] = False
    anc_state['N'] = (anc_state['A']+anc_state['T']+anc_state['C']+anc_state['G'] == 0)
    anc_state['all'] = np.full(anc_state['A'].shape, True)
    return(anc_state)

def anc_states_from_file(c, path_anc, path_fas, node="inner.20"):
    '''This function takes a chrommosome c, the path to the ancestral calls,  and the
     path to the reference fastas, and returns a dictionary with states in the ancestor
     as keys, and boolean np arrays indicating which positions were that state'''
    states=["A", "T", "C", "G"]
    anc = pd.read_csv(path_anc+c+".calls.txt", header=None, skiprows=1, sep="\t",
             names = ["chromosome","position","leaf.1","leaf.2","leaf.3","leaf.4","leaf.5","leaf.6","leaf.7","leaf.8","leaf.9","leaf.10","leaf.11","inner.12","inner.13","inner.14","inner.15", "inner.16","inner.17","inner.18","inner.19","inner.20","hg18","rheMac2"])
    anc["position"] = anc["position"]-1
    anc_alleles = anc[node].to_numpy()
    #replace g/t G/C etc with N
    np.place(anc_alleles, np.invert(np.isin(anc_alleles,states)), "N")
    anc_pos = anc["position"]
    record = SeqIO.read(path_fas+c+".fa.masked", "fasta")
    seq = np.array(list(str(record.seq)))
    anc_states = np.copy(seq)
    anc_states[anc_pos] = anc_alleles
    return(anc_states,seq)

