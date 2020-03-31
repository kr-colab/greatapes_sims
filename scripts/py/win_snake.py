import itertools
import pandas as pd
import numpy as np
import string
import random
import datetime as dt
import math

def expand_grid(data_dict):
    """Create a dataframe from every combination of given values."""
    rows = itertools.product(*data_dict.values())
    return pd.DataFrame.from_records(rows, columns=data_dict.keys())

def id_generator(size=15, chars=string.ascii_uppercase + string.digits):
    return ''.join(random.choice(chars) for _ in range(size))

def get_par_string(row, col_names=["siminterval", "L", "recfile", "exonfile", "mu","delprop", "delcoef","posprop", "poscoef", "N", "gens"]):
    row_values = row.values.astype('str').tolist()
    return(' '.join(["-d "+col_names[i]+"=\\\""+row_values[i]+"\\\"" for i in                  range(len(col_names))]))

seed=1384393
random.seed(seed)
np.random.seed(seed)
#these are the paths to all files/dir we will need
#this table contains info for all edges we will simulate
edges_path ="/home/murillor/projects/greatapes_sims/meta/edges_meta.tsv"
meta_path = "/home/murillor/projects/greatapes_sims/meta/sims/"
out_path = "../../output/"
#path to tsv file with recombination rates
rec_file = "/home/murillor/projects/greatapes_sims/meta/rec_rate_hg18.txt"
#path to tsv with exon annotations
ex_file = "/home/murillor/projects/greatapes_sims/meta/exons_hg18.bed"
chr_file = "/home/murillor/projects/greatapes_sims/meta/hg18.chrom.sizes"
overl_path = "/home/murillor/projects/greatapes_sims/scripts/py/overlay.py"
#path to slim recipe
recipe_path = "/home/murillor/projects/greatapes_sims/scripts/slim/recipe_test.slim"
nreps=1
rescale_factor = 1
burn_gen = 2 # burn_gen N gens of burnin

#setting some parameters
chrom=["chr12"]
padding=np.array([0])
sample_size=10
win_size=1000
L=list(win_size+2*padding)
total_mu = 1.66e-8/rescale_factor
# proportions should be fractions
delprops = [0.6]
posprops = [0]
delcoefs = [-0.03]
poscoefs = [0]
sel_params = {"L":L,"delprop":delprops, "delcoef":delcoefs, "posprop":posprops, "poscoef":poscoefs}
# cols with params
pars = ["siminterval", "L", "recfile", "exonfile", "mu","delprop", "delcoef", "posprop","poscoef", "N", "gens"]
# suck up params and hash them
edges_meta = pd.read_csv(edges_path,sep="\t")
edges_meta["edge"] = edges_meta["edge"].str.replace('_','-')
edges_meta["parent"] = edges_meta["parent"].str.replace('_','-')
#print(edges_meta)
edges_info = edges_meta[["edge","parent","N","gens"]].copy()

#edges_info.gens=1000 #this is here just for testing purposes
edges_info.N=np.ceil(edges_info.N*rescale_factor).astype(int)
rescale_factor = str(rescale_factor)
edges_info.iloc[0,3] = edges_info.iloc[0,2]*burn_gen

#getting chromosome sizes
chr_sizes = pd.read_csv(chr_file,sep="\t", header=None)
chr_sizes.columns=["chr","length"]

#making a data frame that is going to hold all combinations of parameters
tmp=pd.DataFrame()
for i, row in edges_info.iterrows():
    row=row.to_dict()
    for key in row:
        row[key] = [row[key]]
    row.update(sel_params)
    tmp = pd.concat([tmp,expand_grid(row)])

#these boolean masks are here bc some parameter combs are nonsensical (e.g., no sim should have 0 proportion of positive mutations and pos coeff non-zero)
both_p = np.logical_and(tmp.posprop>0,tmp.poscoef > 0)
both_p_zero = np.logical_and(tmp.posprop==0,tmp.poscoef == 0)
both_d = np.logical_and(tmp.delprop>0,tmp.delcoef < 0)
both_d_zero = np.logical_and(tmp.delprop==0,tmp.delcoef == 0)
use_param_comb = np.logical_and(np.logical_or(both_p,both_p_zero),np.logical_or(both_d, both_d_zero))
tmp = tmp[use_param_comb]

windowed=pd.DataFrame()
for c in chrom:
    clen=chr_sizes[chr_sizes.chr==c].length.item()
    windows = []
    start=0
    stop=clen
    last=False
    for window_start in range(start, stop, win_size):
        # determine window stop
        window_stop = window_start + win_size
        if window_stop >= stop:
            # last window
            window_stop = stop
            last = True
        else:
            window_stop -= 1
        windows.append([window_start, window_stop])
        if last:
            break
    windows = pd.DataFrame(windows)
    windows.columns=["start","end"]
    windows['chr']=c
    windowed=pd.concat([windowed,windows])
windowed_tmp = (
      tmp.assign(key=1)
      .merge(windowed.assign(key=1), on="key")
      .drop("key", axis=1)

)
windowed_tmp=windowed_tmp.loc[windowed_tmp.start<2.1*win_size]
tmp = windowed_tmp
#putting everything on our master data.frame (all params)
#tmp["L"] = L
tmp['padding'] = (tmp.L-win_size)/2
tmp["mu"] = (tmp.posprop+tmp.delprop)*total_mu
tmp["siminterval"] = np.where(tmp.N>500000,"500","")
tmp["exonfile"] = ex_file
tmp["recfile"] = rec_file
tmp["numid"] = tmp.groupby(["chr","start","end","L","delprop","delcoef","posprop","poscoef"]).grouper.label_info
rands=np.array([id_generator() for i in range((tmp.numid.max()+1))])
tmp["rand_id"] = rands[tmp["numid"].tolist()]
tmp["outfile_pre"] = tmp.edge+"_"+tmp.rand_id+"_"+tmp.chr+"_"+tmp.start.astype('str')+"_"+tmp.end.astype('str')
assert not tmp.outfile_pre.duplicated().any(), "one of the outfile names is duplicated"
tmp.outfile_pre=out_path+tmp.outfile_pre
tmp.reset_index(drop=True, inplace=True)
#finally assembling our parameter string which will be passed to the slim recipe
tmp["par_string"]=tmp[pars].apply(get_par_string, axis=1)

# writing tmp to file
today = dt.datetime.today().strftime('%m%d%Y')
tmp["date"] = today
tmp.to_csv("../../output/sims_info.tsv", sep="\t", header=False, index=False, mode="a")

# I'm getting the terminals bc not all trees are supposed to be overlayed with mutations and recapped
terminals=list(set(edges_meta.edge.tolist()) - set(edges_meta.parent.tolist()))
terminals.sort()
#terminals=["eastern-chimp"]
terminal_prefixes = np.array(tmp[tmp.edge.isin(terminals)].outfile_pre)

#dealing with replicates
replicated = np.repeat(terminal_prefixes, nreps)
pad=len(str(nreps))
uniqs=np.unique(replicated, return_counts=True)[1]
assert len(set(uniqs)) == 1, "Not all sims are being replicated equally."
rep_suf = np.tile(["rep"+str(i).zfill(pad) for i in range(uniqs[0])], len(uniqs))
all_replicates = [f"{x1}_{x2}.trees" for x1,x2 in zip(replicated,rep_suf)]

# overlay filenames
overl = [f"{x1}_{x2}_overl.trees" for x1,x2 in zip(replicated,rep_suf)]

# single stats filenames
single = [f"{x1}_{x2}_singlestats.tsv" for x1,x2 in zip(replicated, rep_suf)]

# multi pop stats filenames
print(terminals)
pairs = list(itertools.combinations(terminals,2))
pair_strings = ["sp1_"+spp[0]+"_sp2_"+spp[1] for spp in pairs]
print(pairs)
multi = [out_path+s+"_"+rid+"_rep"+str(rep).zfill(pad)+"_multistats.tsv" for s in pair_strings for rid in rands for rep in range(nreps)]

# this is here to make sure only treeseqs that are not inputs to the branch rule make it to be overlayed
ruleorder: root > branch > overlay > single_pop_stats > multi_pop_stats

rule all:
    input:  overl#single+multi #+ tmp[tmp.edge.isin(terminals)].overl_out.tolist()

rule multi_pop_stats:
    input: spp1="../../output/{edge1}_{rand_id}_rep{rep}_overl.trees", spp2="../../output/{edge2}_{rand_id}_rep{rep}_overl.trees"
    output: "../../output/sp1_{edge1}_sp2_{edge2}_{rand_id}_rep{rep}_multistats.tsv"
    params: win_size=win_size, n=sample_size, L = lambda wildcards: tmp[(tmp.edge==wildcards.edge1) & (tmp.rand_id==wildcards.rand_id)].L.item()
    resources: mem_mb=32000, runtime=10*24*60
    threads: 2
    shell: "python multipop_stats_from_trees.py {input.spp1},{input.spp2} {output} {wildcards.edge1},{wildcards.edge2} {wildcards.rand_id} {wildcards.rep} {params.win_size} {params.L} {params.n}"

rule single_pop_stats:
    input: "../../output/{edge}_{rand_id}_rep{rep}_overl.trees"
    output: "../../output/{edge}_{rand_id}_rep{rep}_singlestats.tsv"
    params: win_size=win_size,  n=sample_size, L = lambda wildcards: tmp[(tmp.edge==wildcards.edge) & (tmp.rand_id==wildcards.rand_id)].L.item()
    benchmark: "../../benchmarks/{edge}_{rand_id}_rep{rep}.singlestats.benchmark.txt"
    resources: mem_mb=24000, runtime=10*24*60
    threads: 2
    shell: "python stats_from_tree.py {input} {output} {wildcards.edge} {wildcards.rand_id} {wildcards.rep} {params.win_size} {params.L} {params.n}"

rule overlay:
    input: "../../output/{edge}_{rand_id}_rep{rep}.trees"
    params: mut_rate=total_mu, recapN=edges_meta[edges_meta.edge=="great-apes"].N.item(),rec_hap_path=rec_hap_file, ex_file_path=ex_file, sel_mut_rate=lambda wildcards: tmp[(tmp.edge==wildcards.edge) & (tmp.rand_id==wildcards.rand_id)].mu.item()
    output: "../../output/{edge}_{rand_id}_rep{rep}_overl.trees"
    resources: mem_mb=128000, runtime=10*24*60
    threads: 2
    benchmark: "../../benchmarks/{edge}_{rand_id}_rep{rep}.overl.benchmark.txt"
    shell: "python overlay.py {input} {output} {params.mut_rate} {params.recapN} {params.rec_hap_path} {params.ex_file_path} {params.sel_mut_rate}"

rule root:
    input: lambda wildcards: "../../meta/"+tmp[(tmp.edge=="great-apes") & (tmp.rand_id==wildcards.rand_id)].chr+
#### parei aquiii!
    params: recipe = recipe_path, s = lambda wildcards: tmp[(tmp.edge=="great-apes") & (tmp.rand_id==wildcards.rand_id)].par_string.item(), rescale=rescale_factor
    output: "../../output/great-apes_{rand_id}_rep{rep}.trees"
    output: "../../output/great-apes_{rand_id}_rep{rep}.trees"
    benchmark: "../../benchmarks/great-apes_{rand_id}_rep{rep}.slim.benchmark.txt"
    log: "great-apes_{rand_id}_rep{rep}.log"
    resources: mem_mb=56000, runtime=10*24*60
    threads: 2
    shell:  "slim -m -t {params.s} -d outfile='\"{output}\"' -d rescf='{params.rescale}' {params.recipe}"

rule branch:
    input: lambda wildcards: "../../output/"+edges_info[edges_info.edge==wildcards.edge].parent.item()+"_"+wildcards.rand_id+"_rep"+wildcards.rep+".trees"
    params: recipe=recipe_path,s=lambda wildcards: tmp[(tmp.edge==wildcards.edge) & (tmp.  rand_id==wildcards.rand_id)].par_string.item(), rescale=rescale_factor
    wildcard_constraints: edge="(?!great-apes).*", rep="[0-9]+"
    output: "../../output/{edge}_{rand_id}_rep{rep}.trees"
    benchmark: "../../benchmarks/{edge}_{rand_id}_rep{rep}.slim.benchmark.txt"
    resources: mem_mb=64000, runtime=10*24*60
    threads: 2
    shell: "slim -m -t {params.s} -d path_population_tree='\"{input}\"' -d outfile='\"{output}\"' -d rescf='{params.rescale}' {params.recipe}"
