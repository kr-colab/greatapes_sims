library(reticulate)

library(data.table)

library(argparse)


use_condaenv("apes")
use_python("/opt/anaconda3/envs/apes/bin/python")

np <- reticulate::import("numpy", convert=T)

neutsup="O814WK8MN3UOF27JIR"

filepath = paste0("../../output/joined_stats/sup-rand-id_", neutsup,"_rep_0_win-size_1000000_sample-size_10.npz")


parser <- ArgumentParser(description='Get plots from npz file with dxy matrix, labels and windows')
parser$add_argument('npzfilepath',  type="character", help="Path to .npz file")
parser$add_argument('outpath',  type="character", help="Output path")


filepath = parser$parse_args()$npzfilepath
outpath = parser$parse_args()$outpath


colorder=c("value", "stat", "spp1", "spp2", "chr", "start", "end", "n_acc_bases")

sup_rand_id = strsplit(filepath, "_")[[1]][3]
win_size = as.integer(strsplit(filepath, "_")[[1]][7])

npz = np$load(filepath)
npz$files
dxy_mat = npz$f[["dxy"]]
coord_starts = npz$f[["coord_windows"]]
chroms = npz$f[["chroms"]]
labels = matrix(npz$f[["labels"]], nrow=ncol(dxy_mat), byrow=T) # fixing row to column major
labels = t(apply(labels, 1, sort)) # alphabetical sorting species names
rownames(dxy_mat) = paste0("row", 1:nrow(dxy_mat))
colnames(dxy_mat) = paste0("col", 1:ncol(dxy_mat))
rownames(labels) = paste0("col", 1:ncol(dxy_mat)) #row in labels is col in dxy_mat
colnames(labels) = c("spp1", "spp2")
names(coord_starts) = paste0("row", 1:nrow(dxy_mat))
coords = data.table(rn = names(coord_starts), chr=chroms, start=coord_starts, end = coord_starts+win_size)
labels = data.table(labels, variable=rownames(labels))
dxy = melt(data.table(dxy_mat, keep.rownames = TRUE), id.vars = c("rn"))

dxy = dxy[labels, on=.(variable)]
dxy = dxy[coords, on=.(rn)]
dxy$stat = fifelse(dxy$spp1==dxy$spp2, "pi", "dxy")
dxy[,c("rn", "variable") := NULL]
dxy$start = as.integer(dxy$start)
dxy$end = as.integer(dxy$end)
dxy$n_acc_bases = dxy$end - dxy$start

write_table(dxy[,colorder], paste0(outpath,sup_rand_id,"-pidxy_win-size",win_size,".tsv"), sep="\t", quote=FALSE, row.names=FALSE)