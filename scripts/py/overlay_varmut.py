import tskit
import msprime
import numpy as np
import argparse
import pandas as pd
from helper_functions import *

def remove_msp_mutations(ts):
    ts = pyslim.update(ts)
    tables = ts.dump_tables()
    tables.mutations.clear()
    max_id = -1
    for mut in ts.mutations():
        for d in mut.derived_state.split(","):
            max_id = max(max_id, int(d))
        if mut.metadata['mutation_list'][-1]["mutation_type"] <= 2:
            assert mut.metadata['mutation_list'][-1]["selection_coeff"] != 0
            tables.mutations.append(mut)
        else:
            assert mut.metadata['mutation_list'][-1]["selection_coeff"] == 0
    tables.compute_mutation_parents()
    tables.compute_mutation_times()
    tables.sort()
    new_ts = tables.tree_sequence()
    assert ts.num_mutations-new_ts.num_mutations, "No mutations were removed!"
    print(f"# of mutations before: {ts.num_mutations}, # of muts after: {new_ts.num_mutations}, # of removed: {ts.num_mutations-new_ts.num_mutations}", flush=True)
    return new_ts, max_id

## stat win size must agree with mut win size

parser = argparse.ArgumentParser(description='Gets stats from unioned tree sequence. Assumes ')
parser.add_argument('infilepath', type=str)
parser.add_argument('outfilepath', type=str)
parser.add_argument('rand_id', type=str)
parser.add_argument('win_size', type=lambda x: int(float(x)))
parser.add_argument('total_mut_rate', type=float)
parser.add_argument('region_rate', type=float)
parser.add_argument('sigma', type=float, help="rates are sampled from (1+N(0,sigma))*total_mut_rate")
parser.add_argument('coords_dict', type=str, help="String of a dictionary with padded and non-padded start and ends of the chromosomic region. Assumes one chromosome only!")
parser.add_argument('--seed', type=int, default=2735, required=False)

#args = vars(parser.parse_args(["../../output/4JU4ZW0RTTFYVAS/4JU4ZW0RTTFYVAS_rep0.union.recap.mut.trees", 
#                               "../../output/varmut/4JU4ZW0RTTFYVAS/rand-id_4JU4ZW0RTTFYVAS_rep_0_mut-win-size_1000000_sigma_0.016_total-mut-rate_2e-08.union.recap.mut.trees", 
#                               "4JU4ZW0RTTFYVAS",
#                               "1000000",
#                               "2e-08",
#                               "1.2001e-08",
#                               "0.016",
#                               "{'chr': 'chr12', 'start': 0, 'end': 132000000, 'padded_start': 0.0, 'padded_end': 132000000}"]))


#args

args = vars(parser.parse_args())

coords= eval(args['coords_dict'])
region_rate=float(args["region_rate"])

start = int(coords["start"] - coords["padded_start"])
end = int(coords["end"]-coords["start"] + start)
clen = int(coords["padded_end"]-coords["padded_start"])
print(end, start, flush=True)

ex_path = f"../../output/maps/{args['rand_id']}_exons.tsv"
# exons file
exons = pd.read_csv(ex_path,sep="\t",header=None)
# removing extraneous columns?
exons = exons.iloc[:,:3]
breaks = exons.iloc[:,1:].to_numpy()

#getting win breaks
breakpoints = []
if start > 0:
    breakpoints.append(0)
breakpoints += list(range(start,clen+1, args["win_size"]))
assert breakpoints[-1] <= clen
if breakpoints[-1] < clen:
    breakpoints.append(clen)
win_breaks = np.array(breakpoints)

# drawing mutation rates
rng = np.random.default_rng(args["seed"])
coefs = rng.normal(0,args["sigma"], len(win_breaks)-1)
rates = (1+coefs)*args["total_mut_rate"]

all_breaks = np.unique(np.sort(np.concatenate([win_breaks, breaks.flatten()])))
all_rates = np.full(all_breaks.shape[0]-1, 0, dtype=np.float32)
# filling out all_rates with the rate from the corresponding window
all_rates = rates[np.searchsorted(win_breaks, all_breaks)[1:]-1]
# subtracting the non-neutral rate from the exons
for i,j in np.searchsorted(all_breaks, breaks):
    all_rates[i:j] -= region_rate
#np.savez("test_rates.npz", win_breaks, breaks, all_breaks, rates, all_rates, region_rate)

assert len(all_rates) == len(all_breaks)-1
# just making sure they all got filled
assert np.all(all_rates>-0.5)
# making sure that we don't get a lot of rates zeroe'd out because this would interfere with the variance in neut mut rate
assert np.sum(rates-region_rate<0) < 0.05*len(rates)
all_rates[all_rates<0] = 0.0
if region_rate > 1e-12:
    assert np.sum(rates*np.diff(win_breaks)) > np.sum(all_rates*np.diff(all_breaks))
print(all_rates, all_breaks, flush=True)
rmap = msprime.RateMap(position=all_breaks, rate=all_rates)

#loading ts and remutating it
ts = tskit.load(args["infilepath"])
slim_gen = int(float(ts.provenance(ts.num_provenances-2).record.split("start_time")[-1].split(" ")[1][:-1]))
ts, _ = remove_msp_mutations(ts)

print("Mutating!", flush=True)
model = msprime.SLiMMutationModel(type=3, next_id=1)
ts = msprime.sim_mutations(ts, start_time=slim_gen, model=model, rate=args["total_mut_rate"], keep=True)
print(f"# of muts after first remutating: {ts.num_mutations}", flush=True)
ts = msprime.sim_mutations(ts, end_time=slim_gen, model=model, rate=rmap, keep=True)
print(f"# of muts after secondremutating: {ts.num_mutations}", flush=True)


ts.dump(args["outfilepath"])
