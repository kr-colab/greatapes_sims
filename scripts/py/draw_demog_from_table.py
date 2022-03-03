tree_table_path = "../../data/meta/edges_meta.tsv"
burnin = 2 # 2*N burn in
out_nwk_fname = "../../output/great_apes.nwk"

import pandas as pd
import msprime
import demes
import demesdraw

df = pd.read_csv(tree_table_path, sep="\t")
df.loc[0, "gens"] = burnin * df.loc[0].N

def df_to_nwk(focal, df, parent_col="parent", edge_col="edge", time_col="edge_age_kya"):
    focal_row = df[df[edge_col] == focal]
    focal_time = focal_row[time_col].item()
    tree_str = f'{focal}:{focal_time}'
    is_root = focal_row["parent"].isnull().item()
    if is_root:
        tree_str += ";"
    childs = df[df[parent_col] == focal]
    if len(childs) == 0:
        return tree_str
    else:
        assert len(childs) == 2
        c1, c2 = childs[edge_col]
        c1_str = df_to_nwk(c1, df)
        c2_str = df_to_nwk(c2, df)
        tree_str = f'({c1_str},{c2_str}){tree_str}'
    return tree_str

nwk_tree = df_to_nwk("great_apes", df)
with open(out_nwk_fname, "w") as nwk_fh:
    print(nwk_tree, file=nwk_fh)

pop_sizes = {row["edge"]: row["N"] for i, row in df.iterrows()}

positions = {
                "great_apes": 6.15625,
                "orangutans": 8.5,
                "bornean_orangutan": 8,
                "sumatran_orangutan": 9,
                "african_apes": 3.8125,
                "gorilla": 6.5,
                "eastern_gorilla": 6,
                "western_gorila": 7,
                "human_pan": 1.125,
                "humans": 0,
                "pan": 2.25,
                "bonobo": 1,
                "chimps": 3.5,
                "eastern_central": 2.5,
                "nigerian_western": 4.5,
                "central_chimp": 2,
                "eastern_chimp": 3,
                "nigerian_chimp": 4,
                "western_chimp": 5
                
            }
demog = msprime.Demography.from_species_tree(nwk_tree, initial_size = pop_sizes)
dgraph = demog.to_demes()
dgraph.time_units = "Kya"
w = 1.2 * demesdraw.utils.size_max(dgraph)
for k, v in positions.items():
    positions[k] = v * w
fig, ax = demesdraw.utils.get_fig_axes(scale=1.85)
ax = demesdraw.tubes(dgraph, ax, positions=positions) 
ax.figure.savefig("../../output/tubes_demog_great_apes.pdf")
