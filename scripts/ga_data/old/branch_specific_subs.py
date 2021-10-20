import pandas as pd
import numpy as np

path_anc="/home/murillor/projects/greatapes/data/anc_allele/"
path_out="/home/murillor/projects/greatapes/data/subs/"
nodes = pd.read_csv("/home/murillor/projects/greatapes/data/meta/node_parent.tsv", sep="\t")

# looping through internal nodes and getting the branch subs
for i, row in nodes.loc[nodes.parent.notnull()].iterrows():
    node_col = row['node']
    node_name = row['name']
    parent_col = row['parent']
    parent_name = nodes.loc[nodes.node==parent_col].name.item()
    for chrom in ["chr"+str(i+1) for i in range(22)]:
        ancestral = pd.read_csv(f"{path_anc}{chrom}.calls.txt", header=None, skiprows=1, sep="\t",
             names = ["chromosome","position","leaf.1","leaf.2","leaf.3","leaf.4","leaf.5",
                      "leaf.6","leaf.7","leaf.8","leaf.9","leaf.10","leaf.11","inner.12",
                      "inner.13","inner.14","inner.15", "inner.16","inner.17","inner.18",
                      "inner.19","inner.20","hg18","rheMac2"])
        # test if state at focal node is different then state at ancestral
        # ignoring case (lower vs upper reflect confidence in the call)
        ancestral[parent_col]=ancestral[parent_col].str.upper()
        ancestral[node_col]=ancestral[node_col].str.upper()
        if parent_col=="inner.20":
	    #where it is ambiguous, change it to rheMac2 state
            #but only if that state is within the possible inferred states for the great apes node
            rhe_in_parent = ancestral.apply(lambda x: x['rheMac2'] in x[parent_col], axis=1)
            ancestral.loc[rhe_in_parent,parent_col]=ancestral.loc[rhe_in_parent,"rheMac2"]
        is_diff = ancestral[node_col]!=ancestral[parent_col]
        # making sure there are no ambiguities (A/T)
        is_node_letter = ancestral[node_col].str.contains(r'^[a-zA-Z]$')
        is_parent_letter = ancestral[parent_col].str.contains(r'^[a-zA-Z]$')
        is_all=np.logical_and(np.logical_and(is_diff,is_node_letter),is_parent_letter)
        substitutions = ancestral.loc[is_all][["chromosome","position",parent_col,node_col]]
        substitutions.columns = ["chromosome","position",parent_name, node_name]
        substitutions.to_csv(f"{path_out}{chrom}_substitutions_{parent_name}_{node_name}.tsv", sep="\t",index=False)
        # getting a vcf formatted with ref as ancestral and alt as derived
        vcf = substitutions.copy()
        vcf.columns = ["#CHROM", "POS", "REF", "ALT"]
        vcf['ID'] = '.'
        vcf = vcf[["#CHROM", "POS", "ID", "REF", "ALT"]]
        vcf['QUAL'] = 40
        vcf['FILTER'] = "PASS"
        vcf['INFO'] = ""
        vcf.to_csv(f"{path_out}{chrom}_substitutions_{parent_name}_{node_name}.vcf", sep="\t",index=False)
