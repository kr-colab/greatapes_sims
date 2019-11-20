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



def windowed_pi_from_ts(ts,win_size,L):
    '''
    This function takes a tree sequence, window size and total lenght of chromosome
    and returns an array of windowed pi'''
    hap = allel.HaplotypeArray(ts.genotype_matrix())
    #print((hap.shape), flush=True)
    pos=np.array([s.position for s in ts.sites()])
    ac = hap.count_alleles()
    pi, windows, n_bases, counts = allel.windowed_diversity(pos, ac, size=win_size, start=1, stop=L)
    return(pi)

def ac_from_ts(ts):
    '''
    This function takes a tree sequence,  and returns tuple with allele counts and positions'''
    hap = allel.HaplotypeArray(ts.genotype_matrix())
    #print((hap.shape), flush=True)
    pos=np.array([s.position for s in ts.sites()])
    ac = hap.count_alleles()
    return(ac, pos)

def remove_mutations(ts, start, end, proportion):
    '''
    This function will return a new tree sequence the same as the input,
    but after removing each non-SLiM mutation within regions specified in lists
    start and end with probability `proportion`, independently. So then, if we
    want to add neutral mutations with rate 1.0e-8 within the regions and 0.7e-8
    outside the regions, we could do
      ts = pyslim.load("my.trees")
      first_mut_ts = msprime.mutate(ts, rate=1e-8)
      mut_ts = remove_mutations(first_mut_ts, start, end, 0.3)
    :param float proportion: The proportion of mutations to remove.
    '''
    new_tables = ts.dump_tables()
    new_tables.mutations.clear()
    mutation_map = [-1 for _ in range(ts.num_mutations)]
    for j, mut in enumerate(ts.mutations()):
        keep_mutation = True
        for i in range(len(start)):
            left = start[i]
            right = end[i]
            assert(left < right)
            if i > 0:
                assert(end[i - 1] <= left)
            if mut.position >= left and mut.position < right and len(mut.metadata) == 0:
                keep_mutation = (random.uniform(0, 1) > proportion)
        if keep_mutation:
            mutation_map[j] = new_tables.mutations.num_rows
            if mut.parent < 0:
                new_parent = -1
            else:
                new_parent = mutation_map[mut.parent]
            new_tables.mutations.add_row(site = mut.site, node = mut.node,
                    derived_state = mut.derived_state,
                    parent = new_parent,
                    metadata = mut.metadata)
    return new_tables.tree_sequence()

def win_pi_sims(path, neut_mut, n_pops, n_sims, T, win_size):
    foname = os.path.basename(path[:-1])
    print(("Base filename:"+foname), flush=True)
    s1=timer()
    if (sys.argv[4] == "gene"):
        print(("Mode: gene locations from file"), flush=True)
        matches = re.match( r'N_(\d+)_mu_(.*)_r_(.*)_deff_(\d+)_L_(\d+)', foname)
        N, del_mut, r, deff, L = matches.groups()
        del_mut = float(del_mut)
        L = int(L)
        if del_mut > 0 :
            start_end = np.loadtxt("/home/murillor/projects/slim_sims/start_end_pos_cons.txt", dtype="int")
            start = start_end[:,0]
            end = start_end[:,1]
    else:
        print(("Mode: uniform gene distribution"), flush=True)
        matches = re.match( r'N_(\d+)_mu_(.*)_r_(.*)_deff_(\d+)_L_(\d+)_L0_(.+)_L1_(\d+)_m_.*', foname)
        N, del_mut, r, deff, L, L0, L1 = matches.groups()
        L, L0, L1 = int(L), int(L0),int(L1)
        del_mut = float(del_mut)
        if del_mut > 0:
            start = []
            end = []
            for i in range(0, L, (L0+L1)):
                start.append(i)
                end.append((i+L1)-1)
            #print((start), flush=True)
            #print((end), flush=True)
    s2 = timer()
    print(("Preparing gene locations... Time elapsed (min):"+str(round((s2-s1)/60,3))), flush=True)
    x = np.arange(n_pops)
    combs = list(itertools.combinations(x, 2))
    pis=np.zeros((len(T),n_sims,n_pops,int(L/win_size)))
    div=np.zeros((len(T), n_sims,len(combs),int(L/win_size)))
    for t in range(len(T)):
        for i in range(n_sims):
            filename = glob(path+str(T[t])+"N_sim_"+str(i)+"_RAND_*[0-9].trees")[0]
            if not os.path.isfile(filename[:-6]+"_overlaid.trees"):
                s1 = timer()
                ts_slim = pyslim.load(filename).simplify()
                ts_mut = msprime.mutate(ts_slim,neut_mut, keep=True)
                s2 = timer()
                print(("Mutated in msprime... Time elapsed (min):"+str(round((s2-s1)/60,3))), flush=True)
                if del_mut > 0:
                    #print(("entrei"), flush=True)
                    s1=timer()
                    ts = remove_mutations(ts_mut, start, end, del_mut/neut_mut)
                    s2 = timer()
                    print(("Removed extra mutations in genic regions... Time elapsed (min):"+str(round((s2-s1)/60,3))), flush=True)
                else:
                    ts = ts_mut
                ts.dump(p)
            else:
                print(("Loading overlaid", filename[:-6], "tree sequence..."), flush=True)
                ts = pyslim.load(filename[:-6]+"_overlaid.trees").simplify()
            #print(("Pi0: ", ts.pairwise_diversity(samples=ts.samples(population=0)),"Pi1: ", ts.pairwise_diversity(samples=ts.samples(population=1))), flush=True)
            s1 = timer()
            for j in range(n_pops):
                pop_ac=ac_from_ts(ts.simplify(ts.samples(population=j)))
                #pis[t,i,j,:] = windowed_pi_from_ts(ts.simplify(ts.samples(population=j)),win_size,L)
                pi, windows, n_bases, counts = allel.windowed_diversity(pop_ac[1], pop_ac[0], size=win_size, start=1, stop=L)
                pis[t,i,j,:] = pi
            s2 = timer()
            print(("Calculating windowed Pi... Time elapsed (min):"+str(round((s2-s1)/60,3))), flush=True)
            hap = allel.HaplotypeArray(ts.genotype_matrix())
            pos=np.array([s.position for s in ts.sites()])

            #ach = hap.count_alleles(subpop=ts.samples(population=0))
            #tmp = ac_from_ts(ts.simplify(ts.samples(population=0)))
            #act = tmp[0]
            #post= tmp[1]
            #print(("ac from ts shape:", act.shape), flush=True)
            #print(("pos from ts shape:", post.shape), flush=True)
            #print(("mean pi from ts shape:", np.mean(allel.windowed_diversity(post, act, size=win_size, start=1, stop=L)[0])), flush=True)
            #print(("ac from subpop:", ach.shape), flush=True)
            #print(("pos from subpop:", pos.shape), flush=True)
            #print(("mean pi from subpop:", np.mean(allel.windowed_diversity(pos, ach, size=win_size, start=1, stop=L)[0])), flush=True)
            s1=timer()
            for k in range(len(combs)):
                ac1 = hap.count_alleles(subpop=ts.samples(population=combs[k][0]))
                ac2 = hap.count_alleles(subpop=ts.samples(population=combs[k][1]))
                dxy, windows, n_bases, counts = allel.windowed_divergence(pos, ac1, ac2, size=win_size, start=1, stop=L)
                div[t,i,k,:] = dxy
            s2 = timer()
            print(("Calculating windowed Dxy... Time elapsed (min):"+str(round((s2-s1)/60,3))), flush=True)

    s1 = timer()
    print((pis.shape), flush=True)
    print((div.shape), flush=True)
    output = open(path+foname+'_pis.pkl', 'wb')
    pickle.dump(pis, output)
    output.close()
    output = open(path+foname+'_div.pkl', 'wb')
    pickle.dump(div, output)
    output.close()

    plt.subplot(2, 1, 1)
    plt.plot(np.transpose(pis[0,0,:]), "-")
    plt.title('0N after split')
    plt.ylabel('Pi')
    plt.subplot(2, 1, 2)
    plt.plot(np.transpose(pis[10,0,:]), "-")
    plt.title('10N after split')
    plt.xlabel('Window')
    plt.ylabel('Pi')
    plt.tight_layout()
    plt.savefig(path+foname+'_landscape.pdf')
    plt.close()

    s2 = timer()
    print(("Saving stats and plots to file... Time elapsed (min):"+str(round((s2-s1)/60,3))), flush=True)

#sys.stdout.flush()

s1 = timer()
path = sys.argv[1]
win_size=int(sys.argv[2])
n_sims=int(sys.argv[3])

neut_mut = 1e-8
n_pops = 2

#T=np.arange(0.1,1.1,step=0.1)
#T=np.concatenate([T,np.array([2.0])])
#T=np.around(T, 1)
#T=np.arange(0,11,step=1)
T=np.concatenate([np.arange(0,2.2,0.2), np.array([10.])])
T = [float("%.1f"%t) for t in T]

s2 = timer()
print(("Initializing... Time elapsed (min):"+str(round((s2-s1)/60,3))), flush=True)

win_pi_sims(path, neut_mut, n_pops, n_sims, T, win_size)

"""

path = "/Users/murillo/Drive/phd/w19/rotation/trees/"
folders = glob("N_*_mu_*_L0_*")
paths = [path+folders[i]+"/" for i in range(len(folders))]

total_mut = 1e-8
n_pops = 2
n_sims=1
win_size = 500000
T=np.arange(1,11,step=1)

for i in range(len(paths)):
    path = paths[i]
    foname = os.path.basename(path[:-1])
    print((foname), flush=True)
    matches = re.match( r'N_(\d+)_mu_(.*)_r_(.*)_deff_(\d+)_L_(\d+)_L0_(.+)_L1_(\d+)', foname)
    N, mu, r, deff, L, L0, L1 = matches.groups()
    del_mut = mu
    win_pi_sims(path, neut_mut, del_mut, n_pops, n_sims, T, win_size)"""

#python3 mimulus_win_pi.py /Users/murillo/Drive/phd/w19/rotation/trees/N_10000_mu_0_r_2e-08_deff_0_L_20000000_L0_0_L1_0_m_0/ 500000 1

