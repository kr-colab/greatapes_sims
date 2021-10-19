import argparse
import numpy as np

def npzfiles_to_dict(npzs):
    return {key: [npz[key] for npz in npzs] for key in npzs[0].keys()}

parser = argparse.ArgumentParser(description='Join stats from split up simulations')
parser.add_argument('outpath', type=str, help="out file name (.npz)")
parser.add_argument('paths', type=str, help="paths to the npz files to be joined", nargs="+")

args = vars(parser.parse_args())

np_objs = [np.load(p) for p in args['paths']]

# asserting npz objects have the same contents
# and labels are the same
for i in range(len(np_objs)-1):
    assert np_objs[i].keys() == np_objs[i+1].keys()
    assert (np_objs[i]['labels'] == np_objs[i+1]['labels']).all()

np_dict = npzfiles_to_dict(np_objs)

merged_dict = {k: np.concatenate(v) if k != "labels" else v[0] for k,v in np_dict.items() }

# saving to output
np.savez(args['outpath'], chroms=merged_dict["chrom"], windows=merged_dict["windows"], coord_windows=merged_dict["coord_windows"], dxy=merged_dict["dxy"], labels=merged_dict["labels"])
