### Bulk level analysis ###
### DE analysis for pseudobulk RA 3dst ### 

# removes all glob variables from env
rm(list = ls())

# Loading libraries
suppressMessages(suppressWarnings(library(scran,warn.conflicts = F, quietly = T)))
suppressMessages(suppressWarnings(library(ggplot2,warn.conflicts = F, quietly = T)))
suppressMessages(suppressWarnings(library(cowplot,warn.conflicts = F, quietly = T)))
suppressMessages(suppressWarnings(library(pheatmap,warn.conflicts = F, quietly = T)))
suppressMessages(suppressWarnings(library(stringr,warn.conflicts = F, quietly = T)))
suppressMessages(suppressWarnings(library(grid,warn.conflicts = F, quietly = T)))

# Read the fuctions file 
setwd("/Users/svickovi/Library/Mobile Documents/com~apple~CloudDocs/Desktop/morphoSPOT/3dst_repo/3dst/data")
source('../functions/Read_Functions.R')

# Where are you raw csv expression files located? 
path_samples = list.files("../data/", pattern = glob2rx("RA*norm.exp*"))

# Set output directory that will contain all plots and output gene files
path_output = "/Users/svickovi/Library/Mobile Documents/com~apple~CloudDocs/Desktop/morphoSPOT/plot_outputs/"

stats = matrix(ncol = 6, nrow = 1)
for (ra in unique(sapply(str_split(path_samples, "_"),"[[",1))){
  
  # List avaiable raw matrices for normalization 
  files_norm = path_samples[grep(path_samples, pattern = ra)]
  
  # Load all raw counts matrices in your env
  sample_counter = 1
  for (i in files_norm){
    print(i)
    assign(paste0(ra, "_bulk_mat", sample_counter), as.matrix(rowMeans(get(load(i)))))
    sample_counter = sample_counter + 1
  }
  raw_mat_loaded_in_env <- grep(paste0(ra, "_bulk_"),names(.GlobalEnv),value=TRUE)
  ravg_tmp = do.call(mbind,  mget(raw_mat_loaded_in_env))
  ravg = as.matrix(rowMeans(ravg_tmp))
  assign(paste0(ra, "_mat_avg"), ravg)
  
}  
# Combine all the norm data in one pseudobulk matrix
pseudo_mat_loaded_in_env <- grep("bulk_mat", names(.GlobalEnv),value=TRUE)
RA.norm = do.call(mbind, mget(pseudo_mat_loaded_in_env))
RA.norm[is.na(RA.norm)] <- 0
colnames(RA.norm) = pseudo_mat_loaded_in_env

# Get varible genes 
all.exp.values.norm = RA.norm

# this takes out per section variance
# make a sce object
section_labels = sapply(str_split(pseudo_mat_loaded_in_env, "_"),"[[",1)
sce <- SingleCellExperiment(list(logcounts=as.matrix(all.exp.values.norm)), # this is based on norm counts so object names is logcounts
                            colData=DataFrame(section=section_labels),
                            rowData=DataFrame(label=row.names(all.exp.values.norm)))

# get variable genes
print("Getting variable genes ...")
dec3 <- modelGeneVar(sce, block=sce$section)

# some QC plots to evaluate variance per section
per.block <- dec3$per.block
par(mfrow=c(3, 2))
for (i in seq_along(per.block)) { # check that variacle is correctly computed when taking into account sections/experiments
  decX <- per.block[[i]]
  plot(decX$mean, decX$total, xlab="Mean log-expression",
       ylab="Variance", main=names(per.block)[i])
  curve(metadata(decX)$trend(x), col="blue", add=TRUE)
}
dev.off()

# Get the top genes
top.hvgs <- getTopHVGs(dec3, n=2000)

# Subset on only variable genes
all.exp.values.norm = all.exp.values.norm[rownames(all.exp.values.norm) %in% top.hvgs,]
all.exp.values.norm = genes_cleanup(all.exp.values.norm)

# Hierachical clustering on tsne matrix
hc = hclust(dist(t(all.exp.values.norm)), method = "ward.D2")
f = 5
clusters = as.matrix(cutree(hc, f))
col.panel = c("#FED8B1", "gray85", "gray76",  "gray61", "gray90", "gray99")
row.names(clusters) = colnames(all.exp.values.norm)

print("Plotting dendogram ... Please adjust k and rerun if needed")
hc_plot = plot(hc, labels = F, sub="", xlab = "", main = "Pseudobulk Cluster Dendogram")
rect.hclust(hc, f, border = col.panel[1:f])

m3_ref = all.exp.values.norm # precaution
colnames(m3_ref)

# Split into pseudobulk cluster groups
for (i in 1:f){
  assign(paste0("mg", i), m3_ref[,clusters == i])
  secs = sapply(strsplit(colnames(get(paste0("mg", i))), "_"), "[[",1)
  
  for (j in 1:length(secs)){
    mtmp =  get(paste0("mg", i))[,grep(j, colnames(get(paste0("mg", i))), value = T)]
    write.table(mtmp, file = paste0(path_output, "PseudoBulk_Cluster",i,".", j,".csv"), sep =",", quote =F, row.names = F, col.names = T)
  }
  
}

# Do bimodial test where degg[,1] is the p-value and degg[,2] is the relative difference between the groups
genes = ""
for (i in 1:f){
  print(paste0("Peforming DE analysis for cluster: ", i))
  degg <- bimod.diffExp.test(get(paste0("mg", i)), m3_ref[,!colnames(m3_ref) %in% colnames(get(paste0("mg", i)))], row.names(get(paste0("mg", i))))
  assign(paste0("clus", i, ".genes"), rownames(subset(degg, degg[,1] < 0.001 & degg[,2] > 0.5)))
  assign(paste0("clus", i, ".genes.values"), degg[get(paste0("clus", i, ".genes")),])
  genes = c(genes, get(paste0("clus", i, ".genes")))
  write.table(get(paste0("clus", i, ".genes.values")), file = paste0(path_output, "DEGs.", "PseudoBulk_Cluster",i, ".csv"), quote = F, row.names = T, col.names = T, sep =",")
}

# Clean up gene names (remove house keepers)
genes = genes[-1]
house.keeping.genes = row.names(read.delim("../data/house_keeping_final.txt", row.names = 1)) # found in ./data
genes <- intersect(genes, outersect(genes, row.names(house.keeping.genes)))

# PCA
pc <- prcomp(m3_ref[genes,])
rot <- pc[2]
rot = as.data.frame(rot)
rot = rot[,1:2]
colnames(rot) = c("PCA1","PCA2")
rot$RA = sapply(strsplit(row.names(rot), "_"),"[[",1)

ggplot(rot, aes(x=PCA1, y=PCA2, color=RA)) + 
  geom_point() +
  scale_color_manual(name = "RA patient biopsy", values = c("firebrick4","firebrick","firebrick3","firebrick2","firebrick1"))

# How much variation is explained with each PC? 
summary(pc)
eigs <- pc$sdev^2
eigs[1] / sum(eigs)
eigs[2] / sum(eigs)
eigs[3] / sum(eigs)

# Plot Heatmap
quantile.range <- quantile(m3_ref, probs = seq(0, 1, 0.05))
palette.breaks <- seq(quantile.range["0%"], quantile.range["100%"], 0.1)

# use http://colorbrewer2.org/ to find optimal divergent color palette (or set own)
col.pal = colorRampPalette(c("black","purple","darkorchid1","gold1"))
color.palette  <- col.pal(length(palette.breaks) - 1)

annotdf <- data.frame(row.names = rownames(rot), RA = rot$RA)
cols = c("firebrick4","firebrick","firebrick3","firebrick2","firebrick1")
names(cols) = paste0("RA", 1:5)
mycolors <- list(RA = cols)

heat <- pheatmap(m3_ref[genes,], show_colnames = F , cluster_rows = T, cluster_cols = T, 
                 col=color.palette, breaks = palette.breaks,clustering_distance_cols = "euclidean",
                 scale="none", trace = "none",density.info = "none", cexRow = 0.05,
                 annotation_col = annotdf, annotation_colors = mycolors, fontsize_row = 10,
                 cellheight = 2, cellwidth = 5)

add.flag(heat,
         kept.labels = c("IGLL5", "FN1", "CD52", "MS4A1", "CXCL12", "MMP9", "CLU", "VIM", "COL6A2", "CD4", "XBP1", "JUN", "LTB", "TYROBP"),
         repel.degree = 0)







