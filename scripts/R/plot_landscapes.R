library(renv)
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
library(GGally)

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
    if(!is.na(pair[2])) {
        ancs2 = get_ancestors(pair[2], pair[2], edges)
        mrca = intersect(ancs1, ancs2)[1]
        tmrca = (sum(sapply(ancs1, get_time, edges)[1:(which(ancs1==mrca)-1)]) + sum(sapply(ancs2, get_time, edges)[1:(which(ancs2==mrca)-1)]))
    } else {
        mrca = pair[1]
        tmrca = sapply(ancs1, get_time, edges)[1]
    }
    return(c(mrca, tmrca))
}

rid_to_desc = function(rid, rep, win_size, sep="_") {
    row = sims_info[sims_info$sup_rand_id==rid,][1,]
    return(paste("rep", rep, "stat-win-size", win_size, "resc", row$rescf, "sim-win-len", row$win_len, "pad", row$padding, "del-rate", row$delrate, "del-s", row$delcoef, "pos-rate", row$posrate, "pos-s", row$poscoef, sep=sep))
}


use_condaenv("apes")
use_python("/opt/anaconda3/envs/apes/bin/python")

np <- reticulate::import("numpy", convert=T)
print("igothere")

neutsup="O814WK8MN3UOF27JIR"
sbgssup="4GROMS6FHMNW5X9L72"
sbgsspossup="GGKXIAL2IXDKAJYXAN"
filepath = paste0("../../output/joined_stats/sup-rand-id_", sbgsspossup,"_rep_0_win-size_1000000_sample-size_10.npz")


parser <- ArgumentParser(description='Get plots from npz file with dxy matrix, labels and windows')
parser$add_argument('npzfilepath',  type="character")

filepath = parser$parse_args()$npzfilepath
filepath=filepath



outpath = "../../output/figs/"
edges = fread("../../meta/edges_meta.tsv", sep="\t", header=T, fill=T)
edges$edge = str_replace(edges$edge, "_", "-")
edges$parent = str_replace(edges$parent, "_", "-")
sims_info = read.table("../../output/sims_info.tsv", sep="\t")
header = read.table("../../output/header_sims_info.tsv", header=T)
colnames(sims_info) = colnames(header)
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
x=sapply(unique(dxy$combo), function(x) get_mrca_tmrca(strsplit(x, split="_")[[1]]))
mrcas =data.table(t(x), keep.rownames = TRUE)
colnames(mrcas) = c("combo", "mrca", "tmrca")
mrcas[, tmrca := as.numeric(tmrca)]
setkey(mrcas, mrca)
setkey(edges, edge)
mrcas = mrcas[edges[,c("edge", "ancestral_coding")],]
dxy=mrcas[dxy, on=.(combo)]
dxy$combo = factor(dxy$combo, levels = mrcas$combo[order(ifelse(mrcas$combo==mrcas$mrca, 0.1, 1)*mrcas$ancestral_coding, mrcas$tmrca)], labels=str_replace(mrcas$combo[order(ifelse(mrcas$combo==mrcas$mrca, 0.1, 1)*mrcas$ancestral_coding, mrcas$tmrca)], "_", " ")) # levels ordered by within vs between spp comparison, tmrca, and ancestral coding (tree order)
u=unique(mrcas[,c("mrca", "ancestral_coding", "tmrca")])
dxy$mrca = factor(dxy$mrca, levels=unique(u$mrca[order(u$tmrca, u$ancestral_coding)]))
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
dxy[, ex_overlap := 100*ex_overlap]

## REMOVING 0 rec rate window
dxy=dxy[dxy$mean_rec!=0,]

dxy[,.(cor = cor(value, mean_rec, use="complete.obs", method="spearman")), by=.(spp1, spp2)]
dxy[,.(cor = cor(value, ex_overlap, use="complete.obs", method="spearman")), by=.(spp1, spp2)]

cor_rec = dxy %>%
    group_by(spp1,spp2) %>% 
    do(tidy(cor.test(.$value, .$mean_rec, method="spearman")))

cor_ex = dxy %>%
    group_by(spp1,spp2) %>% 
    do(tidy(cor.test(.$value, .$ex_overlap, method="spearman")))


#### HUMANS ONLY###

humans = dxy[dxy$spp1=="humans"|dxy$spp2=="humans",]
humans$combo=droplevels(humans$combo)
humans$combo = factor(humans$combo, levels=c("humans", "central-chimp humans", "eastern-chimp humans", "humans nigerian-chimp", "humans western-chimp", "bonobo humans", "eastern-gorilla humans", "humans western-gorila", "bornean-orangutan humans", "humans sumatran-orangutan"), labels=c("humans", "humans central-chimp", "humans eastern-chimp", "humans nigerian-chimp", "humans western-chimp", "humans bonobo", "humans eastern-gorilla", "humans western-gorila", "humans bornean-orangutan", "humans sumatran-orangutan"))
drop_subspp = c("humans eastern-chimp", "humans nigerian-chimp", "humans western-chimp", "humans bonobo", "humans eastern-gorilla","humans bornean-orangutan")
humans=humans[!humans$combo %in% drop_subspp,]
humans$combo=droplevels(humans$combo)
ggplot(data=humans, aes(y=value, x=ex_overlap, col=mean_rec)) + geom_point(size=1) + facet_wrap(vars(combo), scales="free_y", labeller = label_wrap_gen(width=18)) + theme_bw(base_size=24) + labs(y="Pi and Dxy", x="% Exon", col="Rec rate", subtitle=desc) + scale_colour_viridis_c(option="magma", direction=-1)
ggsave(filename=paste0(outpath,"exon-scatter-humans-subset_",desc,"_",sup_rand_id,".pdf"), width = 20, height = 4, units = "in")

ggplot(data=humans, aes(y=value, x=mean_rec, col=ex_overlap)) + geom_point(size=1) + facet_wrap(vars(combo), scales="free_y", labeller = label_wrap_gen(width=18), nrow=1) + theme_bw(base_size=24) + labs(y="Pi and Dxy", x="Rec rate", col="% Exon", subtitle=desc) + scale_colour_viridis_c(direction=-1)
ggsave(filename=paste0(outpath,"rec-scatter-humans-subset_",desc,"_",sup_rand_id,".pdf"), width = 20, height = 4, units = "in")

col5palette = c("#b96440", "#864bc6", "#749c51", "#b24c82","#6d81b0")
droplevs = c("western-chimp", "nigerian-chimp", "eastern-chimp", "eastern-gorilla", "bornean-orangutan")
sub=dxy[(!dxy$combo %in% droplevs) & (dxy$stat=="pi"),]
sub$combo = droplevels(sub$combo)
p_pi = ggplot(data = sub, aes(x=start, y=value, group=combo)) + geom_line(aes(col=combo)) + scale_colour_manual(values=col5palette) + theme_bw(base_size=24) + labs(y="Diversity", x="Window", col="Species") + guides(colour = guide_legend(override.aes = list(size=5))) + theme(legend.position = "bottom")
p_pi
ggsave(filename=paste0(outpath,"pi-landscape_subset_chr12_",desc,"_",sup_rand_id,".pdf"), width = 14, height = 7, units = "in")

t=dcast(dxy[dxy$stat=="pi" & !(dxy$spp1 %in% droplevs) & !(dxy$spp2 %in% droplevs) ,], formula=start+end~combo, value.var='value')
cormat=cor(t[,-c("start", "end")], use="complete.obs", method="spearman")
library(ggtree)
library(colorspace)
tree <- read.tree(text = "(((humans:3765,(central-chimp:871,bonobo:871):2894):1854,western-gorila:5669):5519,sumatran-orangutan:11138);")
tree$edge.length=(tree$edge.length*1000)/20
ptips=ggtree(tree, ladderize=F) + geom_tiplab(size=4) + geom_treescale(x=0,y=0,width=100000)
gheatmap(ptips, cormat, offset=340000, low="white", high="black", colnames_position = "bottom", font.size=4, colnames_angle=45, colnames_offset_y=-.2, width=1.2) + scale_fill_continuous_diverging(palette="Blue-Red 3",limits=c(-1,1)) + labs(fill="Corr")

ggplot(data = reshape2::melt(cormat), aes(x=Var1, y=Var2, fill=value)) + geom_tile() + scale_fill_continuous_diverging(palette="Blue-Red 3",limits=c(-1,1))+ theme_minimal(base_size=14) + labs(title="Corr between Pi/Dxy",subtitle=spaced_desc, fill="Corr") + xlab("") + ylab("")  + theme(axis.text.x = element_text(angle = 45, vjust = 1.1, hjust=1))
ggsave(filename=paste0(outpath,"pi-cormat_genome-wide_subset_",desc,"_",sup_rand_id,".pdf"), width = 5.8, height = 5.2, units = "in")

########


ggplot(data=dxy, aes(y=value, x=ex_overlap, col=mean_rec)) + geom_point() + facet_wrap(vars(combo), scales="free_y", labeller = label_wrap_gen(width=18)) + theme_bw() + labs(y="Pi/Dxy", x="% Exon", col="Rec rate", title=spaced_desc) + scale_colour_viridis_c(option="magma", direction=-1)
ggsave(filename=paste0(outpath,"exon-scatter_",desc,"_",sup_rand_id,".pdf"), width = 500, height = 400, units = "mm")

ggplot(data=dxy, aes(y=value, x=mean_rec, col=ex_overlap)) + geom_point() + facet_wrap(vars(combo), scales="free_y", labeller = label_wrap_gen(width=18)) + theme_bw() + labs(y="Pi/Dxy", x="Rec rate", col="% Exon", title=spaced_desc) + scale_colour_viridis_c(direction=-1)
ggsave(filename=paste0(outpath,"rec-scatter_",desc,"_",sup_rand_id,".pdf"), width = 500, height = 400, units = "mm")

p_ex = ggplot(cor_ex, aes(x=p.value)) + geom_histogram() + theme_bw() + labs(y="Count", x="P-value", title="Correlation between pi/dxy and % overlap with exons")

p_rec = ggplot(cor_rec, aes(x=p.value)) + geom_histogram() + theme_bw() + labs(y="Count", x="P-value", title="Correlation between pi/dxy and recombination rate")

(p_ex / p_rec) + plot_annotation(title=spaced_desc)


ggsave(filename=paste0(outpath,"pval-dist_",desc,"_",sup_rand_id,".pdf"), width = 210, height = 297, units = "mm")

col10palette = c("#cb5e95", "#9bd345", "#9242c5", "#7ecf93", "#665ea9", "#cba657", "#523240", "#c3533b", "#99afc0", "#53673c")
col9palette = c("#7dcd5b","#904cc2","#cbb354","#4b2f51","#8bc6af","#c65381","#575c3a","#c2593b","#9898c4")
p_pi = ggplot(data = dxy[dxy$stat=="pi"], aes(x=start, y=value, group=combo)) + geom_line(aes(col=combo)) + scale_colour_manual(values=col10palette) + theme_bw(base_size=12) + labs(y="Diversity", x="Window", col="Species") + theme(legend.margin=margin(t=0, r=0, b=0, l=-1, unit="cm"))
p_pi

p_dxy_mrca = ggplot(data = dxy[dxy$stat=="dxy"], aes(x=start, y=value, group=combo)) + geom_line(aes(col=mrca)) + scale_colour_manual(values=col9palette) + theme_bw(base_size=12) + labs(y="Divergence", x="Window", col="MRCA") + theme(legend.margin=margin(t=0, r=0, b=0, l=-1, unit="cm"))
p_dxy_mrca

p_dxy_tmrca = ggplot(data = dxy[dxy$stat=="dxy"], aes(x=start, y=value, group=combo)) + geom_line(aes(col=tmrca, linetype=mrca)) + scale_colour_viridis_c(direction=-1) + theme_bw(base_size=12) + labs(y="Divergence", x="Window", col="TMRCA") + theme(legend.margin=margin(t=0, r=0, b=0, l=-1, unit="cm"))
p_dxy_tmrca

(p_pi / p_dxy_mrca / p_dxy_tmrca) + plot_annotation(title=spaced_desc)

cor_rec_ex = with(dxy[dxy$combo=="humans",],cor.test(mean_rec,ex_overlap, use="complete.obs"))

coly1 = "#8a501a"
coly2 = "#7a387c"
p_rec_ex = ggplot(data=dxy[dxy$combo=="humans",], aes(x=start)) + geom_line(aes(y=mean_rec), col=coly1) + geom_line( aes(y=ex_overlap / 2), col=coly2) + scale_y_continuous(name = "Rec rate", sec.axis = sec_axis(~.*2, name="% Exon")) + theme_bw(base_size = 12) + theme(axis.title.y = element_text(color = coly1), axis.title.y.right = element_text(color = coly2))  + annotate("text", label=paste0("Corr = ", round(cor_rec_ex$estimate, 2)), x=1e8, y = 5.5) + xlab("Window") + theme(legend.margin=margin(t=0, r=0, b=0, l=-1, unit="cm"))

(p_pi / p_dxy_mrca / p_rec_ex) + plot_annotation(title=spaced_desc, theme=theme(title = element_text(size=8))) + plot_layout(heights=c(3,3,1))

(p_pi / p_dxy_mrca) + plot_annotation(title=spaced_desc)

ggsave(filename=paste0(outpath,"pidxy-landscape_",desc,"_",sup_rand_id,".pdf"), width = 210, height = 297, units = "mm")

with(dxy[dxy$combo=="humans",],
     plot(start,mean_rec))
with(dxy[dxy$combo=="humans",],
     plot(start,ex_overlap))
with(dxy[dxy$combo=="humans",],
     plot(start,mean_rec))
with(dxy[dxy$combo=="humans",],
     plot(mean_rec,ex_overlap))

with(dxy[dxy$combo=="humans",],cor.test(mean_rec,ex_overlap, use="complete.obs"))$estimate


library("colorspace")
dxy$combo = factor(dxy$combo, levels = str_replace(mrcas$combo[order(mrcas$combo!=mrcas$mrca, mrcas$ancestral_coding)], "_", " "), ) # levels ordered by within vs between spp comparison, tmrca, and ancestral coding (tree order)
widedxy = dcast(dxy, start + end ~ combo, value.var="value")
cormat = cor(widedxy[,-c("start", "end")])
ggplot(data = reshape2::melt(cormat), aes(x=Var1, y=Var2, fill=value)) + geom_tile() + scale_fill_continuous_diverging(limits=c(-1,1))+ theme_bw(base_size=10) + labs(title="Corr between Pi/Dxy",subtitle=spaced_desc, fill="Corr") + xlab("") + ylab("")  + theme(title=element_text(size=8), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(filename=paste0(outpath,"cor-mat_",desc,"_",sup_rand_id,".pdf"), width = 210, height = 190, units = "mm")
