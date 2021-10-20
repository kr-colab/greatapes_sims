#!/usr/bin/env python3
import pandas as pd
from functools import reduce
from Bio import SeqIO
import numpy as np
import sys

chr = "chr"+sys.argv[1]

Homo = pd.read_csv("Homo_alleles.txt."+chr, sep="\t", header=None, names=["chr", "pos", "alleles"])
Pan_troglodytes = pd.read_csv("Pan_troglodytes_alleles.txt."+chr, sep="\t", header=None, names=["chr", "pos", "alleles"])
Pan_paniscus = pd.read_csv("Pan_paniscus_alleles.txt."+chr, sep="\t", header=None, names=["chr", "pos", "alleles"])
Gorilla = pd.read_csv("Gorilla_alleles.txt."+chr, sep="\t", header=None, names=["chr", "pos", "alleles"])
Pongo_abelii = pd.read_csv("Pongo_abelii_alleles.txt."+chr, sep="\t", header=None, names=["chr", "pos", "alleles"])
Pongo_pygmaeus = pd.read_csv("Pongo_pygmaeus_alleles.txt."+chr, sep="\t", header=None, names=["chr", "pos", "alleles"])
Macaque = pd.read_csv(chr+"_macaque.txt", sep="\t", header=None, names=["chr", "pos", "alleles"])

dfs = [Homo, Pan_troglodytes, Pan_paniscus, Gorilla, Pongo_abelii, Pongo_pygmaeus, Macaque]

joined = reduce(lambda left,right: pd.merge(left,right,on=["chr", "pos"], how='outer'), dfs)

joined.columns = ["chr", "pos", "Homo", "Pan_troglodytes", "Pan_paniscus", "Gorilla", "Pongo_abelii", "Pongo_pygmaeus","Macaque"]

record = SeqIO.read(chr+".fa.masked", "fasta")

seq = np.array(list(str(record.seq)))

joined.Homo[joined.Homo.isna()]= seq[joined[joined.Homo.isna()].pos-1]
joined.Pan_troglodytes[joined.Pan_troglodytes.isna()]= seq[joined[joined.Pan_troglodytes.isna()].pos-1]
joined.Pan_paniscus[joined.Pan_paniscus.isna()]= seq[joined[joined.Pan_paniscus.isna()].pos-1]
joined.Gorilla[joined.Gorilla.isna()]= seq[joined[joined.Gorilla.isna()].pos-1]
joined.Pongo_abelii[joined.Pongo_abelii.isna()]= seq[joined[joined.Pongo_abelii.isna()].pos-1]
joined.Pongo_pygmaeus[joined.Pongo_pygmaeus.isna()]= seq[joined[joined.Pongo_pygmaeus.isna()].pos-1]

joined.to_csv(chr+"_all_states.tsv", sep="\t", index=False, header=True)

