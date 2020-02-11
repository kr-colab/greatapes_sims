import pyslim
from stats_funcs import *

paths=["../../output/eastern-gorilla_CF96YXF0YAQNNT0_rep04_overl.trees","../../output/sumatran-orangutan_CF96YXF0YAQNNT0_rep04_overl.trees"]

all_ac=[]
all_pos=[]
win_size = 1000
L=1000
for p in paths:
	ts = pyslim.load(p)
	acs, pos = acs_from_ts(ts, 10)
	all_ac.append(acs[0])
	all_pos.append(pos)

x = np.arange(len(all_ac))
combs = list(itertools.combinations(x, 2))

for k in range(len(combs)):
	ac1=all_ac[combs[k][0]]
	pos1=all_pos[combs[k][0]]
	ac2=all_ac[combs[k][1]]
	pos2=all_pos[combs[k][1]]
	pos = np.concatenate((pos1,pos2))
	pos.sort()
	new_ac1 = np.full((len(pos),2),-1)
	new_ac2 = np.full((len(pos),2),-1)
	new_ac1[np.where(np.isin(pos,pos1)==True)] = ac1
	new_ac2[np.where(np.isin(pos,pos2)==True)] = ac2
	n_sampled1 = np.sum(ac1,axis=1)[0]
	n_sampled2 = np.sum(ac2,axis=1)[0]
	# filling out the empty genotypes, which should be fixed to the reference alleles
	# if a different allele fixed, then the genotype matrix return by tskit would be e.g. [0,20]
	# confirm with Peter
	new_ac1[np.where(new_ac1==np.array([-1,-1]))[0]] = np.array([n_sampled1,0])
	new_ac2[np.where(new_ac2==np.array([-1,-1]))[0]] = np.array([n_sampled2,0])
	dxy, windows, n_bases, counts = allel.windowed_divergence(pos, new_ac1,  new_ac2, size=win_size, start=1, stop=L)

