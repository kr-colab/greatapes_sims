library(reticulate)
library(GenomicRanges)
library(data.table)
library(ggplot2)
library(fuzzyjoin)
library(ape)
library(tidyverse)
library(broom)
library(patchwork)
library(argparse)

get_mrca_tmrca = function(pair) {
    
    get_ancestors = function(e, vec, edges) {
        #print(e)
        #print(vec)
        par = edges[edges$edge==e,]$parent
        if (! par %in% edges$edge) {
            return(vec)
        }
        vec = c(vec, par)
        get_ancestors(par, vec, edges)
    }
    
    get_time = function(e, edges) {
        return(edges[edges$edge==e,]$gens)
    }
    
    ancs1 = get_ancestors(pair[1], pair[1], edges)
    ancs2 = get_ancestors(pair[2], pair[2], edges)
    mrca = intersect(ancs1, ancs2)[1]
    tmrca = (sum(sapply(ancs1, get_time, edges)[1:(which(ancs1==mrca)-1)]) + sum(sapply(ancs2, get_time, edges)[1:(which(ancs2==mrca)-1)]))
    return(c(mrca, tmrca))
}

rid_to_desc = function(rid, rep, win_size, sep="_") {
    row = sims_info[sims_info$sup_rand_id==rid,][1,]
    return(paste("rep", rep, "stat-win-size", win_size, "resc", row$rescf, "sim-win-len", row$win_len, "pad", row$padding, "del-rate", row$delrate, "del-s", row$delcoef, "pos-rate", row$posrate, "pos-s", row$poscoef, sep=sep))
}


use_condaenv("greatapes")
use_python("/opt/anaconda3/envs/greatapes/bin/python")

np <- reticulate::import("numpy", convert=T)
print("igothere")

neutsup="O814WK8MN3UOF27JIR"
wbgsonlysup="RD24T8SVUEW8FTANCW"
sbgsonlysup="4GROMS6FHMNW5X9L72"
filepath = "../../output/joined_stats/sup-rand-id_O814WK8MN3UOF27JIR_rep_0_win-size_1000000_sample-size_10.npz"


parser <- ArgumentParser(description='Get plots from npz file with dxy matrix, labels and windows')
parser$add_argument('npzfilepath',  type="character")

filepath = parser$parse_args()$npzfilepath

outpath = "../../output/figs/"
edges = read.table("../../meta/edges_meta.tsv", sep="\t", header=T, fill=T)
edges$edge = str_replace(edges$edge, "_", "-")
edges$parent = str_replace(edges$parent, "_", "-")
sims_info = read.table("../../output/sims_info.tsv", sep="\t")
header = read.table("../../output/header_sims_info.tsv", header=T)
colnames(sims_info) = colnames(header)

sup_rand_id = neutsup
sup_rand_id = strsplit(filepath, "_")[[1]][3]
rep = strsplit(filepath, "_")[[1]][5]
win_size = as.integer(strsplit(filepath, "_")[[1]][7])
desc = rid_to_desc(sup_rand_id, rep, win_size)
spaced_desc = str_replace_all(desc, "_", " ")

rand_id = sims_info[sims_info$sup_rand_id == sup_rand_id,]$rand_id[1]
rec_map_path = paste0("../../meta/maps/",rand_id, "_recrate.tsv")
exon_map_path = paste0("../../meta/maps/",rand_id, "_exons.tsv")
rec = fread(rec_map_path)
colnames(rec) = c("chr", "start", "end", "rate")
exon = fread(exon_map_path, drop = c(4))
colnames(exon) = c("chr", "start", "end")

npz = np$load(filepath)
npz$files
dxy_mat = npz$f[["dxy"]]
coord_starts = npz$f[["coord_windows"]]
labels = matrix(npz$f[["labels"]], nrow=ncol(dxy_mat), byrow=T) # fixing row to column major
labels = t(apply(labels, 1, sort)) # alphabetical sorting species names
rownames(dxy_mat) = paste0("row", 1:nrow(dxy_mat))
colnames(dxy_mat) = paste0("col", 1:ncol(dxy_mat))
rownames(labels) = paste0("col", 1:ncol(dxy_mat)) #row in labels is col in dxy_mat
colnames(labels) = c("spp1", "spp2")
names(coord_starts) = paste0("row", 1:nrow(dxy_mat))
coords = data.table(rn = names(coord_starts), start=coord_starts, end = coord_starts+win_size)
labels = data.table(labels, variable=rownames(labels))
dxy = melt(data.table(dxy_mat, keep.rownames = TRUE), id.vars = c("rn"))

dxy = dxy[labels, on=.(variable)]
dxy = dxy[coords, on=.(rn)]
dxy$stat = fifelse(dxy$spp1==dxy$spp2, "pi", "dxy")
dxy$combo = paste(dxy$spp1, dxy$spp2, sep="_")
dxy[stat == "pi"]$combo = dxy[stat == "pi"]$spp1
dxy[,c("rn", "variable") := NULL]
dxy$start = as.integer(dxy$start)
dxy$end = as.integer(dxy$end)
x=sapply(unique(dxy[dxy$stat=="dxy",]$combo), function(x) get_mrca_tmrca(strsplit(x, split="_")[[1]]))
mrcas =data.table(t(x), keep.rownames = TRUE)
colnames(mrcas) = c("combo", "mrca", "tmrca")
mrcas[, tmrca := as.numeric(tmrca)]
dxy=mrcas[dxy, on=.(combo)]
# joining with rec rate and getting mean rate by window
dxy = as.data.table(interval_left_join(dxy, rec, by = c("start", "end"), minoverlap=2)) # minoverlap bc end is not inclusive
dxy[, percent:= (end.y-start.y)/(end.x-start.x)]
dxy=dxy[, .(mean_rec = weighted.mean(rate, percent)), by=setdiff(colnames(dxy), c("rate", "start.y", "end.y", "chr", "percent"))]
setnames(dxy, c("start.x", "end.x"), c("start","end"))

# joining with exons and getting percent overlap
dxy = as.data.table(interval_left_join(dxy, exon, by = c("start", "end"), minoverlap=2)) # minoverlap bc end is not inclusive
dxy[, percent:= (end.y-start.y)/(end.x-start.x)]
dxy[is.na(dxy$percent)]$percent=0
dxy=dxy[, .(ex_overlap = sum(percent)), by=setdiff(colnames(dxy), c("start.y", "chr", "end.y", "percent"))]
setnames(dxy, c("start.x", "end.x"), c("start","end"))

dxy[,.(cor = cor(value, mean_rec, use="complete.obs")), by=.(spp1, spp2)]
dxy[,.(cor = cor(value, ex_overlap, use="complete.obs")), by=.(spp1, spp2)]

cor_rec = dxy %>%
    group_by(spp1,spp2) %>% 
    do(tidy(cor.test(.$value, .$mean_rec)))

cor_ex = dxy %>%
    group_by(spp1,spp2) %>% 
    do(tidy(cor.test(.$value, .$ex_overlap)))

ggplot(data=dxy, aes(y=value, x=ex_overlap, col=tmrca)) + geom_point() + facet_wrap(vars(spp1, spp2), scales="free_y") + theme_bw() + labs(y="Pi/Dxy", x="% exon overlap", col="TMRCA") + scale_colour_viridis_c(direction=-1)
ggsave(filename=paste0(outpath,"ex-overlap-scatter_",desc,"_",sup_rand_id,".pdf"), width = 500, height = 400, units = "mm")

ggplot(data=dxy, aes(y=value, x=mean_rec, col=tmrca)) + geom_point() + facet_wrap(vars(spp1, spp2), scales="free_y") + theme_bw() + labs(y="Pi/Dxy", x="Rec rate", col="TMRCA") + scale_colour_viridis_c(direction=-1)
ggsave(filename=paste0(outpath,"rec-scatter_",desc,"_",sup_rand_id,".pdf"), width = 500, height = 400, units = "mm")

p_ex = ggplot(cor_ex, aes(x=p.value)) + geom_histogram() + theme_bw() + labs(y="Count", x="P-value", title="Correlation between pi/dxy and % overlap with exons")

p_rec = ggplot(cor_rec, aes(x=p.value)) + geom_histogram() + theme_bw() + labs(y="Count", x="P-value", title="Correlation between pi/dxy and recombination rate")

(p_ex / p_rec) + plot_annotation(title=spaced_desc)

ggsave(filename=paste0(outpath,"pval-dist_",desc,"_",sup_rand_id,".pdf"), width = 210, height = 297, units = "mm")

col10palette = c("#cb5e95", "#9bd345", "#9242c5", "#7ecf93", "#665ea9", "#cba657", "#523240", "#c3533b", "#99afc0", "#53673c")

p_pi = ggplot(data = dxy[dxy$stat=="pi"], aes(x=start, y=value, group=combo)) + geom_line(aes(col=combo)) + scale_colour_manual(values=col10palette) + theme_bw(base_size=12) + labs(y="Diversity", x="Window", col="Species")

p_dxy_mrca = ggplot(data = dxy[dxy$stat=="dxy"], aes(x=start, y=value, group=combo)) + geom_line(aes(col=mrca)) + scale_colour_manual(values=col10palette) + theme_bw(base_size=12) + labs(y="Divergence", x="Window", col="MRCA")


p_dxy_tmrca = ggplot(data = dxy[dxy$stat=="dxy"], aes(x=start, y=value, group=combo)) + geom_line(aes(col=tmrca)) + scale_colour_viridis_c(direction=-1) + theme_bw(base_size=12) + labs(y="Divergence", x="Window", col="TMRCA")


(p_pi / p_dxy_mrca / p_dxy_tmrca) + plot_annotation(title=spaced_desc)

ggsave(filename=paste0(outpath,"pidxy-landscape_",desc,"_",sup_rand_id,".pdf"), width = 210, height = 297, units = "mm")


with(dxy[dxy$combo=="humans",],
     plot(start,mean_rec))
with(dxy[dxy$combo=="humans",],
     plot(start,ex_overlap))
with(dxy[dxy$combo=="humans",],
     plot(start,mean_rec))
with(dxy[dxy$combo=="humans",],
     plot(mean_rec,ex_overlap))
