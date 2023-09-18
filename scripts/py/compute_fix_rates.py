#!/usr/bin/env python
# coding: utf-8

# In[2]:


from datetime import datetime
import numpy as np
def allele_frequencies(ts, sample_sets=None, mode="site"):
    if sample_sets is None:
        sample_sets = [ts.samples()] 
    n = np.array([len(x) for x in sample_sets])
    def f(x):
        return x/n
    return ts.sample_count_stat(sample_sets, f, len(sample_sets), windows='sites', polarised=True, mode=mode, strict=False, span_normalise=False)

def get_sites_coeff(tseq):
    """
    Gets the selection coefficient for the mutations within each site.
    If more than one mutation happened at a site, this will return np.nan.
    Output is a numpy array with shape (ts.num_sites,).
    """
    site_coeff = np.full(tseq.num_sites, 0.0, dtype=float)
    for j, s in enumerate(tseq.sites()):
        coeffs = []
        for m in s.mutations:
            coeffs.append(stdpopsim.ext.selection_coeff_from_mutation(tseq, m))
        if len(coeffs) >= 1:
            site_coeff[j] = coeffs[0]
            if len(coeffs) > 1:
                site_coeff[j] = np.nan
    return site_coeff


# In[3]:


import tskit
import numpy
from helper_functions import *
import stdpopsim
import pyslim


# In[4]:


rand_id="T7X58A009J7ZPPZ"
rand_id="R9N6TUIB2YIKHVX"
rand_id="R54VZWUI2F4SIK6"
rand_id="8CDA8OUW1O0B2BN"
#rand_id="370AE9Z4DTPI14R"


# In[5]:


popspath = f"../../output/{rand_id}/{rand_id}_rep0.pops"


# In[6]:


tspath = f"../../output/{rand_id}/{rand_id}_rep0.union.recap.mut.trees"


# In[7]:


ex_path = f"../../output/maps/{rand_id}_exons.tsv"


# In[8]:


with open(popspath, "r") as f:
    pops = eval(f.readline())


# In[9]:


exons = pd.read_csv(ex_path,sep="\t", header=None)
# removing extraneous columns?
exons = exons.iloc[:,:3]


# In[10]:


exon_bpoints = exons.iloc[:,1:3].to_numpy().flatten()


# In[11]:


total_bp_exons = np.sum((exons.iloc[:,2]-exons.iloc[:,1]).to_numpy())


# In[12]:


print(f"{datetime.now()} | Total bp exons: {total_bp_exons}", flush=True)


# In[13]:


ts = tskit.load(tspath)


# In[14]:


print(f"{datetime.now()} | ts loaded", flush=True)


# In[15]:


ts = pyslim.update(ts)


# In[16]:


print(f"{datetime.now()} | pyslim update completed", flush=True)



# The highest time for nodes that are samples (these would be at the great apes branch)
timemax = np.max(ts.tables.nodes[ts.tables.nodes.flags==1].time)

#timemax = 1000000 # just to test!!
# In[ ]:


tsd = ts.decapitate(time = timemax)

#tsd = ts
print(f"{datetime.now()} | Decapitation with time {timemax} completed -- prev max root time {ts.max_root_time}, now {tsd.max_root_time}", flush=True)

# 4min
samples = sample_from_ts(tsd, sample_size=10)
print(samples)
print(f"{datetime.now()} | Got samples from decapitated ts", flush=True)



# In[ ]:


is_within_exon = (np.searchsorted(exon_bpoints,tsd.tables.sites.position) % 2) == 1


# In[ ]:


print(f"{datetime.now()} | Built boolean array is_within_exon", flush=True)


# In[ ]:


pos_sites = np.array([
        sum(md["selection_coeff"]>0 for m in s.mutations for md in m.metadata["mutation_list"])
        for s in tsd.sites()
]) > 0
neg_sites = np.array([
        sum(md["selection_coeff"]<0 for m in s.mutations for md in m.metadata["mutation_list"])
        for s in tsd.sites()
]) > 0


# In[ ]:


print(f"{datetime.now()} | Built arrays with pos_sites and neg_sites", flush=True)


# In[ ]:


# 15min
freqs = allele_frequencies(tsd, sample_sets=[samples[5]])
print(f"{datetime.now()} | Computed site allele frequencies", flush=True)
#bfreqs = allele_frequencies(tsd, sample_sets=[samples[5]], mode="branch")


# In[ ]:


#print(f"{datetime.now()} | Computed branch allele frequencies", flush=True)


# In[ ]:


total_fix_rate = np.sum(freqs)/(timemax*tsd.sequence_length)
total_fix_rate
# Rate of fixations per generation per bp across the SLiM history of humans
# Mutation rate ~2e-8


# In[ ]:


print(f"{datetime.now()} | Site total fixation rate: {total_fix_rate}, sum of freqs {np.sum(freqs)}", flush=True)

#print(f"{datetime.now()} | Branch total fixation rate: {(np.sum(bfreqs)/(timemax*tsd.sequence_length))*2e-8}, sum of freqs {np.sum(bfreqs)}", flush=True)

# In[ ]:


pos_fix_rate_in_ex = np.sum(freqs[pos_sites])/(timemax*total_bp_exons)
print(pos_fix_rate_in_ex, 2e-8*0.01)# change the 0.01 to whatever the prop of beneficial was for this sim


# In[ ]:


print(f"{datetime.now()} | Pos fix rate in exons: {pos_fix_rate_in_ex}", flush=True)


# In[ ]:


(pos_fix_rate_in_ex)*total_bp_exons
# number of positive fixations per generation


# In[ ]:


print(f"{datetime.now()} | Number of positive fixations per generation: {(pos_fix_rate_in_ex)*total_bp_exons}", flush=True)


# In[ ]:


fix_rate_in_ex = np.sum(freqs[is_within_exon])/(timemax*total_bp_exons)
fix_rate_in_ex


# In[ ]:


print(f"{datetime.now()} | Fix rate within exons: {fix_rate_in_ex}", 
      flush=True)


# In[ ]:


(pos_fix_rate_in_ex)/fix_rate_in_ex
# what percentage of fixations were driven by pos sel


# In[ ]:


print(f"{datetime.now()} | Percentage of fixations within exon in pos sel sites: {(pos_fix_rate_in_ex)/fix_rate_in_ex}", 
      flush=True)

