# Installing DPT
### DE analysis for RA1 3dst ### 

# Installing DPT if needed
library(devtools)
# devtools::install_url('https://www.helmholtz-muenchen.de/fileadmin/ICB/software/destiny/destiny_1.3.4.tar.gz')
# devtools::install_url('https://www.helmholtz-muenchen.de/fileadmin/ICB/software/dpt_0.6.0.tar.gz')

# Loading libraries
library(dpt)
library(gplots)
library(destiny)
library(scran)
library(plot3D)
library(Rtsne)
library(jackstraw)
library(gdata)
library(ggplot2)
library(akima)
library(matrixStats)

# Source f(x) file
source('../functions/Read_Functions.R') # in ./functions

# Load data as R objects for RA1
load('RA1_norm.exp.values.1') # in ./data
load('RA1_norm.exp.values.2') # in ./data
load('RA1_norm.exp.values.3') # in ./data
load('RA1_norm.exp.values.4') # in ./data
load('RA1_1_selected_adjusted_spots_3D_manual_app') # in ./data
load('RA1_2_selected_adjusted_spots_3D_manual_app') # in ./data
load('RA1_3_selected_adjusted_spots_3D_manual_app') # in ./data
load('RA1_4_selected_adjusted_spots_3D_manual_app') # in ./data

# Combine in one file
RA.norm = cbind(m1, m2, m3, m4) # RA1
m = RA.norm # this is a combined matrix for all sections from RA1 patient with 4 norm expression values (individual files avaiable in ./data)

# Load all infilatrates for RA1
inf_s1 = read.csv("RA1_1_all_inf.csv", header = T, row.names = 1) # files found on SCP
inf_s2 = read.csv("RA1_2_all_inf.csv", header = T, row.names = 1) # files found on SCP
inf_s3 = read.csv("RA1_3_all_inf.csv", header = T, row.names = 1) # files found on SCP
inf_s4 = read.csv("RA1_4_all_inf.csv", header = T, row.names = 1) # files found on SCP

### Take only infiltrates that are present in all sections according to annotations
# eg. inf1_all have values tracking "inf1" accross all 4 sections
# naming in "all_inf.csv" files was arbitrary so new objects are made here
inf1_all = c(row.names(subset(inf_s1, infiltrate == "Inf3")), row.names(subset(inf_s2, infiltrate == "Inf1")),
             row.names(subset(inf_s3, infiltrate == "Inf14")), row.names(subset(inf_s4, infiltrate == "Inf4")))
inf2_all = c(row.names(subset(inf_s1, infiltrate == "Inf4")), row.names(subset(inf_s2, infiltrate == "Inf3")),
             row.names(subset(inf_s3, infiltrate == "Inf15")), row.names(subset(inf_s4, infiltrate == "Inf5")))
inf3_all = c(row.names(subset(inf_s1, infiltrate == "Inf11")), row.names(subset(inf_s2, infiltrate == "Inf7")),
             row.names(subset(inf_s3, infiltrate == "Inf8")), row.names(subset(inf_s4, infiltrate == "Inf11")))
inf4_all = c(row.names(subset(inf_s1, infiltrate == "Inf16" | infiltrate == "Inf17")), row.names(subset(inf_s2, infiltrate == "Inf11")),
             row.names(subset(inf_s3, infiltrate == "Inf4")), row.names(subset(inf_s4, infiltrate == "Inf6")))
inf5_all = c(row.names(subset(inf_s1, infiltrate == "Inf19")), row.names(subset(inf_s2, infiltrate == "Inf13")),
             row.names(subset(inf_s3, infiltrate == "Inf1")), row.names(subset(inf_s4, infiltrate == "Inf1")))

# Subset per infiltrate as the naming ie. Inf1 presumes we are following the same infiltrate (given 3D alignment) in all present sections
RA.norm.inf1 = RA.norm[,colnames(RA.norm) %in% inf1_all]
RA.norm.inf2 = RA.norm[,colnames(RA.norm) %in% inf2_all]
RA.norm.inf3 = RA.norm[,colnames(RA.norm) %in% inf3_all]
RA.norm.inf4 = RA.norm[,colnames(RA.norm) %in% inf4_all]
RA.norm.inf5 = RA.norm[,colnames(RA.norm) %in% inf5_all]

# Make an object with all annotated infiltrates
all_inf = c(row.names(inf_s1),row.names(inf_s2), row.names(inf_s3), row.names(inf_s4))

### Make avg gene exp barplot per inf and section
gen_names = c("CD52","MS4A1","FN1") # list of marker genes

# Calculate avg gene expression per inf and section
inf1_rm = rm_section(RA.norm.inf1[gen_names,]) + 1
inf2_rm = rm_section(RA.norm.inf2[gen_names,]) + 1
inf3_rm = rm_section(RA.norm.inf3[gen_names,]) + 1
inf4_rm = rm_section(RA.norm.inf4[gen_names,]) + 1
inf5_rm = rm_section(RA.norm.inf5[gen_names,]) + 1

# Calculate sd of gene expression per inf and sections
sd1_rm = sd_section(RA.norm.inf1[gen_names,])
sd2_rm = sd_section(RA.norm.inf2[gen_names,])
sd3_rm = sd_section(RA.norm.inf3[gen_names,])
sd4_rm = sd_section(RA.norm.inf4[gen_names,])
sd5_rm = sd_section(RA.norm.inf5[gen_names,])

# Calculate size of gene expression per inf and sections
s1_rm = size_section(RA.norm.inf1[gen_names,])
s2_rm = size_section(RA.norm.inf2[gen_names,])
s3_rm = size_section(RA.norm.inf3[gen_names,])
s4_rm = size_section(RA.norm.inf4[gen_names,])
s5_rm = size_section(RA.norm.inf5[gen_names,])

# Plot for each Inf separately if needed
myData = data.frame(as.numeric(inf1_rm), as.numeric(sd1_rm/sqrt(s1_rm)))
colnames(myData) = c("Avg", "SEM")
row.names(myData) = c(paste("Section1", gen_names, sep ="_"),paste("Section2", gen_names, sep ="_"),
                      paste("Section3", gen_names, sep ="_"),paste("Section4", gen_names,  sep ="_"))
limits <- aes(ymax = myData$Avg + myData$SEM, ymin = myData$Avg - myData$SEM)
p <- ggplot(data = myData, aes(x = factor(sapply(strsplit(row.names(myData), split = "_"),"[[",1)), y = Avg, fill = rep(c("1_CD52","2_MS4A1","3_FN1"), 4)))
p + geom_bar(stat = "identity", position = position_dodge(0.9)) +
  geom_errorbar(limits, position = position_dodge(0.9), width = 0.25) +
  labs(x = "", y = "Avg expression") +
  scale_fill_manual(name = "Gene", values=c("#066799",  "#7392CB","#CDCC63"))

### Combine all the inf norm data in one matrix
inf_all = mbind(RA.norm.inf1, RA.norm.inf2, RA.norm.inf3, RA.norm.inf4, RA.norm.inf5)

# Run DPT on the whole inf matrix
dm = run_dpt(inf_all)

# Do hierarchical clustering
hc = hclust(dist(as.matrix(dm$RA1, dm$RA2, dm$RA3)),method = "ward.D")
clusters = cutree(hc, k = 2) # 2 clusters
names(clusters) = rownames(dm) 
col.clusters = c(rep("#F16049",length(clusters[clusters == 1])), rep("#50BCEC",length(clusters[clusters == 2])))
scatter3D(dm$RA1, dm$RA2, dm$RA3, pch=16, col = col.clusters)

### Do DEG analysis based on Nature Biotechnology 33, 495–502 (2015) doi:10.1038/nbt.3192
inf_all = genes_cleanup(inf_all)
m3_ref = inf_all # Precaution

# Split inro cluster groups
mg1 <- m3_ref[,clusters == 1]
mg2 <- m3_ref[,clusters == 2]

# Do bimodial test where degg[,1] is the p-value and degg[,2] is the relative difference between the groups
degg <- bimod.diffExp.test(mg1, m3_ref[,!colnames(m3_ref) %in% colnames(mg1)], row.names(mg1))
clus1.genes <- rownames(subset(degg, degg[,1] < 0.001 & (degg[,2] > 1 | degg[,2] < -1)))
clus1.genes.values = degg[clus1.genes,]

degg <- bimod.diffExp.test(mg2,m3_ref[,!colnames(m3_ref) %in% colnames(mg2)], row.names(mg2))
clus2.genes <- rownames(subset(degg, degg[,1] < 0.001 & (degg[,2] > 1 | degg[,2] < -1)))
clus2.genes.values = degg[clus2.genes,]

# Clean up gene names (remove house keepers)
genes = unique(c(clus1.genes,clus2.genes))
house.keeping.genes = row.names(read.delim("house_keeping_final.txt", row.names = 1)) # found in ./data
genes <- intersect(genes, outersect(genes,row.names(house.keeping.genes)))

### Make distance heatmaps of clones parallel to tSNE
ra.1 = sapply(strsplit(rownames(dm),"_"),"[[",1)
colfunc <- colorRampPalette(c("blue","yellow" ,"red"))
color.ra.1 = c(rep(colfunc(4)[1], length(ra.1[ra.1 == "X1"])), rep(colfunc(4)[2], length(ra.1[ra.1 == "X2"])),rep(colfunc(4)[3], length(ra.1[ra.1 == "X3"])),rep(colfunc(4)[4], length(ra.1[ra.1 == "X4"])))
names(color.ra.1) = rownames(data) 

annotation_col = as.matrix(c(rep("#f15f48", ncol(mg1)),rep("#4ebceb", ncol(mg2))))
rownames(annotation_col) = colnames(m3_ref)
annotation_col = c(annotation_col[annotation_col == "#f15f48",],annotation_col[annotation_col == "#4ebceb",],annotation_col[annotation_col == "#8a3795",])
double.col = cbind(annotation_col, color.ra.1)

cbind(names(annotation_col) == names(color.ra.1))
ann=matrix(nrow=length(annotation_col), ncol = 2)
rownames(ann) = names(annotation_col)
for (name in names(annotation_col)){
  ann[name,] = cbind(annotation_col[name], color.ra.1[name])
}

#Assign color annotations for plotting in 3D for inf
ann.col.inf = matrix(nrow=ncol(RA.norm), ncol = 1)
rownames(ann.col.inf) = colnames(RA.norm)
for (sam in colnames(RA.norm)){
  if ((!sam %in% row.names(ann)) == TRUE) ann.col.inf[sam,1] = 4
  else {
    if (ann[sam,1] == "#f15f48") { ann.col.inf[sam,1] = 3}
    if (ann[sam,1] == "#4ebceb") { ann.col.inf[sam,1] = 2}
    if (ann[sam,1] == "#8a3795") { ann.col.inf[sam,1] = 1}
  }
}

pdf("Heatmap_Inf.pdf")

source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")

quantile.range <- quantile(m3_ref, probs = seq(0, 1, 0.01))
palette.breaks <- seq(quantile.range["60%"], quantile.range["100%"], 0.1)

# use http://colorbrewer2.org/ to find optimal divergent color palette (or set own)
col.pal = colorRampPalette(c("black","purple","gold1"))
color.palette  <- col.pal(length(palette.breaks) - 1)

heatmap.3(m3_ref[genes, ], Rowv = T, Colv = F, col=color.palette, breaks = palette.breaks,
          scale="none",trace = "none",density.info = "none", ColSideColors=ann)

dev.off()

ann.col = matrix(nrow=nrow(ann), ncol = 1)
rownames(ann.col) = rownames(ann)
for (sam in rownames(ann)){
  if ((ann[sam,1] == "#FED8B1") & (ann[sam,2] == "#FF7F00")){ ann.col[sam,1] = "red"}
  else {ann.col[sam,1] = ann[sam,1] }
}

# Make tSNE clustering on all spots from the RA1 biopsy 
all.exp.values.norm = RA.norm
fit = trendVar(all.exp.values.norm)
decomp = decomposeVar(all.exp.values.norm, fit)
top.hvgs = order(decomp$bio, decreasing=TRUE)
decomp_bio = decomp[top.hvgs,][1:1500,]

# Subset on only variable genes
all.exp.values.norm = all.exp.values.norm[rownames(all.exp.values.norm) %in% rownames(decomp_bio),]
all.exp.values.norm = genes_cleanup(all.exp.values.norm)

# Color points based on only subsetting inflitrates
tsne_inf = all.exp.values.norm[,colnames(all.exp.values.norm) %in% all_inf]
tsne_rest = all.exp.values.norm[,!colnames(all.exp.values.norm) %in% all_inf]
all_matrix = cbind(tsne_inf, tsne_rest)

# How many components are important? 
per_out = permutationPA(all_matrix, B = 10, threshold = 0.1, verbose = TRUE, seed = NULL)

# Run tSNE again
tsne_out <- Rtsne(t(all_matrix), dims = 3, initial_dims = 13, theta = 0.5, check_duplicates = FALSE, pca = TRUE, perplexity = 60, max_iter = 1000, verbose = TRUE) # Run TSNE; initial_dims = 10 in the manuscript

# Make matrix of tsne data output
tsne <- as.data.frame(tsne_out$Y)
t1 <- as.matrix(tsne[1])
t2 <- as.matrix(tsne[2])
t3 <- as.matrix(tsne[3])
tsne_m <- cbind (t1,t2,t3)
rownames(tsne_m) = colnames(all_matrix)

# Hierachical clustering on tsne matrix
hc = hclust(dist(as.matrix(tsne_m[,1], tsne_m[,2])),method = "ward.D2")
clusters = as.matrix(cutree(hc, k = 4))
row.names(clusters) = rownames(tsne_m)

# Assign colors to clusters
ct = 1
col.clusters = matrix (ncol = 1, nrow = nrow(clusters))
for (cl in clusters[,1]){
  if (cl == 1) col.clusters[ct,] = "#FED8B1"
  if (cl == 2) col.clusters[ct,] = "gray61"
  if (cl == 3) col.clusters[ct,] = "gray76"
  if (cl == 4) col.clusters[ct,] = "gray90"
  ct = ct + 1
}
row.names(col.clusters) = row.names(clusters)

# Make tSNE 3D plots for output
scatter3D(tsne_m[,1],tsne_m[,2],tsne_m[,3], phi = 0,theta = -180, axes=FALSE, ann=FALSE,pch=16, col=col.clusters)
col.clusters[row.names(col.clusters) %in% colnames(tsne_inf),] <- "red"
col.clusters[!row.names(col.clusters) %in% colnames(tsne_inf),] <- "grey51"
scatter3D(tsne_m[,1],tsne_m[,2],tsne_m[,3], phi = 0,theta = -180, axes=FALSE, ann=FALSE,pch=16, col=col.clusters)

#Do DEG analysis
m3_ref = all_matrix
mg1 <- m3_ref[,clusters == 1]
mg2 <- m3_ref[,clusters == 2]
mg3 <- m3_ref[,clusters == 3]
mg4 <- m3_ref[,clusters == 4]

degg <- bimod.diffExp.test(mg1, m3_ref[,!colnames(m3_ref) %in% colnames(mg1)], row.names(mg1))
clus1.genes <- rownames(subset(degg, degg[,1] < 0.001 & (degg[,2] > 0.5 | degg[,2] < -0.5)))
clus1.genes.values = degg[clus1.genes,]

degg <- bimod.diffExp.test(mg2,m3_ref[,!colnames(m3_ref) %in% colnames(mg2)], row.names(mg2))
clus2.genes <- rownames(subset(degg, degg[,1] < 0.001 & (degg[,2] > 0.5 | degg[,2] < -0.5)))
clus2.genes.values = degg[clus2.genes,]

degg <- bimod.diffExp.test(mg3,m3_ref[,!colnames(m3_ref) %in% colnames(mg3)], row.names(mg3))
clus3.genes <- rownames(subset(degg, degg[,1] < 0.001 & (degg[,2] > 0.5 | degg[,2] < -0.5)))
clus3.genes.values = degg[clus3.genes,]

degg <- bimod.diffExp.test(mg4,m3_ref[,!colnames(m3_ref) %in% colnames(mg4)], row.names(mg4))
clus4.genes <- rownames(subset(degg, degg[,1] < 0.001 & (degg[,2] > 0.5 | degg[,2] < -0.5)))
clus4.genes.values = degg[clus4.genes,]

### Assign colors and clean up genes
annotation_col = as.matrix(c(rep("#FED8B1", ncol(mg1)),rep("gray61", ncol(mg2)),rep("gray76", ncol(mg3)),rep("gray85", ncol(mg4))))
genes = unique(c(clus1.genes,clus2.genes,clus3.genes,clus4.genes))
house.keeping.genes = row.names(read.delim("house_keeping_final.txt", row.names = 1))
genes <- genes[!genes %in% row.names(house.keeping.genes)]

annotation_col = c(annotation_col[annotation_col == "#FED8B1",],annotation_col[annotation_col == "gray61",],annotation_col[annotation_col == "gray76",],annotation_col[annotation_col == "gray85",])
inf_colors = c(rep("#FF7F00", ncol(tsne_inf)), rep("grey51", ncol(tsne_rest)))
names(annotation_col) = c(colnames(mg1), colnames(mg2), colnames(mg3),colnames(mg4))
names(inf_colors) = colnames(all_matrix)

ann=matrix(nrow=length(annotation_col), ncol = 2) # this ann goes as input to plot.clusters
rownames(ann) = names(annotation_col)
for (name in names(annotation_col)){
  ann[name,] = cbind(annotation_col[name],inf_colors[name])
}
ann_to_plot = ann # this ann goes as input to plot.clusters

pdf("Heatmap_tSNE.pdf")

source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")

quantile.range <- quantile(m3_ref, probs = seq(0, 1, 0.05))
palette.breaks <- seq(quantile.range["60%"], quantile.range["100%"], 0.1)

# use http://colorbrewer2.org/ to find optimal divergent color palette (or set own)
col.pal = colorRampPalette(c("black","purple","darkorchid1","gold1"))
color.palette  <- col.pal(length(palette.breaks) - 1)

heatmap.3(m3_ref[genes,], Rowv = T, Colv = F, col=color.palette, breaks = palette.breaks,
          scale="none",trace = "none",density.info = "none", ColSideColors = ann, labCol = 0, cexRow = 0.2)

dev.off()

# Which barcodes aren't part of the annotated infiltrates? Why? 
bar_clus1 = mg1[genes, colnames(mg1) %in% colnames(tsne_inf)]
ncol(bar_clus1)/length(all_inf) # percentage of total annotated inflitrates
bar_clus2 = mg2[genes, colnames(mg2) %in% colnames(tsne_inf)]
ncol(bar_clus2)/length(all_inf) # percentage of total annotated inflitrates
bar_clus3 = mg3[genes, colnames(mg3) %in% colnames(tsne_inf)]
ncol(bar_clus3)/length(all_inf) # percentage of total annotated inflitrates
bar_clus4 = mg4[genes, colnames(mg4) %in% colnames(tsne_inf)]
1/length(all_inf)# percentage of total annotated inflitrates

# All clusters together == m3_ref
all_bar = cbind(mg1, mg2, mg3, mg4)
colnames(all_bar) = c(rep("B1", ncol(mg1)),rep("B2", ncol(mg2)),
                      rep("B3", ncol(mg3)),rep("B4", ncol(mg4)))

#### Make barplot of avg expression per interseting marker gene
gen_names = c("CCL19","CXCL13","LTB","PRG4","MMP3","CD52","MS4A1","FN1")

# Avg expression and std error
bar_rm = rm_section(all_bar[gen_names,])
bar_sd = sd_section(all_bar[gen_names,])/sqrt(size_section(all_bar[gen_names,]))

# Plot barplot
m_both = bar_rm # Precaution
sd_both = bar_sd # Precaution
myData = data.frame(as.numeric(m_both), as.numeric(sd_both))
rownames(myData) = c(paste(gen_names, "Cluster1", sep="_"), paste(gen_names, "Cluster2", sep="_"),
                     paste(gen_names, "Cluster3", sep="_"),paste(gen_names, "Cluster4", sep="_"))
limits <- aes(ymax = myData$as.numeric.m_both. + myData$as.numeric.sd_both.,
              ymin = myData$as.numeric.m_both. - myData$as.numeric.sd_both.)
colvar = rep(c("#FED8B1","gray61","gray76","gray85"), length(gen_names))
names(colvar) = rownames(myData)
p <- ggplot(data = myData, aes(x = rownames(myData), y = myData$as.numeric.m_both.), fill = colvar)
p + ylab(c("Avg Log10 expression"))+ ylim(-1, 5)+ theme_bw() + geom_bar(stat = "identity", position = position_dodge(width = 0.9), fill = colvar) + geom_errorbar(mapping = limits, stat = "identity", position = position_dodge(width = 0.9), width = 0.3)+ theme(axis.line.y = element_line(colour = "black"), axis.title.x = element_blank(), axis.ticks.x= element_blank(), axis.title.y = element_text(), axis.text.x = element_text(angle = 30),panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())

#Plot all barcode clusters over the tissue for Fig3a
ann.col = matrix(nrow=nrow(ann_to_plot), ncol = 1)
rownames(ann.col) = rownames(ann_to_plot)
for (sam in rownames(ann)){
  if ((ann[sam,1] == "#FED8B1") & (ann[sam,2] == "#FF7F00")){ ann.col[sam,1] = "red"}
  else {ann.col[sam,1] = ann[sam,1] }
}

ann.cluster = matrix(nrow=nrow(ann_to_plot), ncol = 1)
rownames(ann.cluster) = rownames(ann_to_plot)
for (sam in rownames(ann)){
  if ((ann[sam,1] == "#FED8B1") & (ann[sam,2] == "#FF7F00") | (ann[sam,1] == "gray61") & (ann[sam,2] == "#FF7F00") | (ann[sam,1] == "gray76") & (ann[sam,2] == "#FF7F00") | (ann[sam,1] == "gray85") & (ann[sam,2] == "#FF7F00")){ ann.cluster[sam,1] = 5}
  else {
    if (ann[sam,1] == "#FED8B1") { ann.cluster[sam,1] = 4}
    if (ann[sam,1] == "gray61") { ann.cluster[sam,1] = 3}
    if (ann[sam,1] == "gray76") { ann.cluster[sam,1] = 2}
    if (ann[sam,1] == "gray85") { ann.cluster[sam,1] = 1}}
}

#Plot 3D genes as a heatmap (not morphological)
#Plot 3D heatmaps
#RA1
plot.gene.2d.inf.4("RA1", ann.col.inf, m1, m2, m3, m4, s1, s2, s3, s4, x=40, y=20, transparency=1, min=1, max=4)
plot.gene.2d.cluster.4("RA1", ann.cluster, m1, m2, m3, m4, s1, s2, s3, s4, x=40, y=20, transparency=1, min=1, max=5)
plot.gene.2d.4("RA1_Con", "ACTB", m1, m2, m3, m4, s1, s2, s3, s4, x=40, y=20, transparency=1, min=-2, max=6, con = T)
