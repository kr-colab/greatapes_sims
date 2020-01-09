import pandas as pd

def ruleMaker(name):
    print(f'rule {name}_edge:\n'
    f'\tinput: tmp[tmp.edge=="{name}"].path_population_tree.tolist()\n'
    f'\tparams: recipe=recipe_path,'
    f'\tout_p = out_path,s=tmp[tmp.edge=="{name}"].par_string.tolist()\n'
    f'\toutput: tmp[tmp.edge=="{name}"].outfile.tolist()\nshell: "slim'
    f'-m -t {{params.s}} -d '
    f'outfile=\'\\"{{params.out_p}}{{output}}\\"\''
    f'{{params.recipe}}"\n'
         )

edges_path ="/home/murillor/projects/greatapes_sims/meta/edges_meta.tsv"

edges_meta = pd.read_csv(edges_path,sep="\t")

[ruleMaker(n) for n in edges_meta[~edges_meta.parent.isna()].edge.values.tolist()]
