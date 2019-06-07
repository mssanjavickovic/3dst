### DE analysis for RA2 3dst ### 

# Installing DPT
library(devtools)
# devtools::install_url('https://www.helmholtz-muenchen.de/fileadmin/ICB/software/destiny/destiny_1.3.4.tar.gz')
# devtools::install_url('https://www.helmholtz-muenchen.de/fileadmin/ICB/software/dpt_0.6.0.tar.gz')

# Loading libraries
library(dpt)
library(destiny)
library(scran)
library(plot3D)
library(Rtsne)
library(jackstraw)
library(gdata)
library(ggplot2)
library(akima)

# Source f(x) file
source('Read_Functions.R') # in ./functions

# Load data as R objects for RA2
load('RA2_norm.exp.values.1') # in ./data
load('RA2_norm.exp.values.2') # in ./data
load('RA2_norm.exp.values.3') # in ./data
load('RA2_norm.exp.values.4') # in ./data
load('RA2_norm.exp.values.5') # in ./data
load('RA2_norm.exp.values.6') # in ./data
load('RA2_norm.exp.values.7') # in ./data
load('RA2_1_selected_adjusted_spots_3D_manual_app') # in ./data
load('RA2_2_selected_adjusted_spots_3D_manual_app') # in ./data
load('RA2_3_selected_adjusted_spots_3D_manual_app') # in ./data
load('RA2_4_selected_adjusted_spots_3D_manual_app') # in ./data
load('RA2_5_selected_adjusted_spots_3D_manual_app') # in ./data
load('RA2_6_selected_adjusted_spots_3D_manual_app') # in ./data
load('RA2_7_selected_adjusted_spots_3D_manual_app') # in ./data

# Combine in one file
RA.norm = cbind(m1, m2, m3, m4, m5, m6, m7) # RA2
m = RA.norm # this is a combined matrix for all sections from RA1 patient with 4 norm expression values and 7 RA2 sections (individual files avaiable in ./data)

# Load all infilatrates for RA2
inf_s1 = read.csv("RA2_1_all_inf.csv", header = T, row.names = 1) # files found on SCP
inf_s2 = read.csv("RA2_2_all_inf.csv", header = T, row.names = 1) # files found on SCP
inf_s3 = read.csv("RA2_3_all_inf.csv", header = T, row.names = 1) # files found on SCP
inf_s4 = read.csv("RA2_4_all_inf.csv", header = T, row.names = 1) # files found on SCP
inf_s5 = read.csv("RA2_5_all_inf.csv", header = T, row.names = 1) # files found on SCP
inf_s6 = read.csv("RA2_6_all_inf.csv", header = T, row.names = 1) # files found on SCP
inf_s7 = read.csv("RA2_7_all_inf.csv", header = T, row.names = 1) # files found on SCP

# Combine infiltrates from all sections into one object
all.inf = rbind(inf_s1, inf_s2, inf_s3, inf_s4, inf_s5, inf_s6, inf_s7)

# Subset per infiltrate as the naming ie. Inf1 presumes we are following the same infiltrate (given 3D alignment) in all present sections
inf1 = subset(all.inf, infiltrate == "Inf1")
RA.norm.inf1 = RA.norm[,colnames(RA.norm) %in% row.names(inf1)]
inf1.spots = matrix(row.names(inf1))

inf2 = subset(all.inf, infiltrate == "Inf2")
RA.norm.inf2 = RA.norm[,colnames(RA.norm) %in% row.names(inf2)]
inf2.spots = matrix(row.names(inf2))

inf3 = subset(all.inf, infiltrate == "Inf3")
RA.norm.inf3 = RA.norm[,colnames(RA.norm) %in% row.names(inf3)]
inf3.spots = matrix(row.names(inf3))

inf4 = subset(all.inf, infiltrate == "Inf4")
RA.norm.inf4 = RA.norm[,colnames(RA.norm) %in% row.names(inf4)]
inf4.spots = matrix(row.names(inf4))

inf5 = subset(all.inf, infiltrate == "Inf5")
RA.norm.inf5 = RA.norm[,colnames(RA.norm) %in% row.names(inf5)]
inf5.spots = matrix(row.names(inf5))

inf6 = subset(all.inf, infiltrate == "Inf6")
RA.norm.inf6 = RA.norm[,colnames(RA.norm) %in% row.names(inf6)]
inf6.spots = matrix(row.names(inf6))

### Combine all the inf norm data in one matrix
inf_all = mbind(RA.norm.inf1, RA.norm.inf2, RA.norm.inf3, RA.norm.inf4, RA.norm.inf5)

# Run DPT on the whole inf matrix
dm = run_dpt(inf_all)

# Assign colors
ra = sapply(strsplit(rownames(dm),"_"),"[[",1)
colfunc <- colorRampPalette(c("blue", "yellow" ,"red"))
color.ra = c(rep(colfunc(7)[1], length(ra[ra == "X1"])), rep(colfunc(7)[2], length(ra[ra == "X2"])),rep(colfunc(7)[3], length(ra[ra == "X3"])),rep(colfunc(7)[4], length(ra[ra == "X4"])),rep(colfunc(7)[5], length(ra[ra =="X5"])),rep(colfunc(7)[6], length(ra[ra == "X6"])),
                   rep(colfunc(7)[7], length(ra[ra == "X7"])))

# Plot PCA pt results
plot(dm$RA1, dm$RA2, pch = 16, col=c(color.ra))
library(plot3D)
scatter3D(dm$RA1, dm$RA2, dm$RA3, pch=16, col = color.ra) # Color based on sections

# Do hierarchical clustering and plot dendogram
hc = hclust(dist(as.matrix(dm$RA1, dm$RA2, dm$RA3)),method = "ward.D2")
clusters = cutree(hc, k = 3)
col.clusters = c(rep("black",length(clusters[clusters == 1])), rep("blue",length(clusters[clusters == 2])), rep("red",length(clusters[clusters == 3])))
scatter3D(dm$RA1, dm$RA2, dm$RA3, pch=16, col = col.clusters)

# Do DE on the 3 clusters
pdf("HC_Inf_Clusters.pdf", useDingbats = F)
clusters = as.matrix(cutree(hc, k = 3))
row.names(clusters) = rownames(data)
col.clusters = c(rep("#f15f48",length(clusters[clusters == 1])), rep("#4ebceb",length(clusters[clusters == 2])), rep("#8a3795",length(clusters[clusters == 3])))
plot(hc, labels = clusters, cex = 0.5)
dev.off()

# Plot PCA based on HClust
pdf("Scatter_3D_clusters.pdf", useDingbats = F)
scatter3D(dm$RA1, dm$RA2, dm$RA3, pch=16, col = col.clusters)
dev.off()

#Do DE analysis
m3_ref = t(inf_all)

# Remove ribosomal, MT, ambiguous (from HTseq-count) and ENSG (no corresponding Refseq ID)
m3_ref = genes_cleanup(m3_ref)

# Subset data per cluster
mg1 <- m3_ref[,clusters == 1]
mg2 <- m3_ref[,clusters == 2]
mg3 <- m3_ref[,clusters == 3]

# Run DE between clusters in one-vs-all
degg <- bimod.diffExp.test(mg1, m3_ref[,!colnames(m3_ref) %in% colnames(mg1)], row.names(mg1))
clus1.genes <- rownames(subset(degg, degg[,1] < 0.001 & (degg[,2] > 0.5 | degg[,2] < -0.5)))
clus1.genes.values = degg[clus1.genes,]

degg <- bimod.diffExp.test(mg2,m3_ref[,!colnames(m3_ref) %in% colnames(mg2)], row.names(mg2))
clus2.genes <- rownames(subset(degg, degg[,1] < 0.001 & (degg[,2] > 0.5 | degg[,2] < -0.5)))
clus2.genes.values = degg[clus2.genes,]

degg <- bimod.diffExp.test(mg3,m3_ref[,!colnames(m3_ref) %in% colnames(mg3)], row.names(mg3))
clus3.genes <- rownames(subset(degg, degg[,1] < 0.001 & (degg[,2] > 0.5 | degg[,2] < -0.5)))
clus3.genes.values = degg[clus3.genes,]

pdf("Heatmap_Inf.pdf")
library(gplots)
source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")

# Assign color ranges
quantile.range <- quantile(m3_ref, probs = seq(0, 1, 0.01))
palette.breaks <- seq(quantile.range["60%"], quantile.range["100%"], 0.1)

# Assign colors to clusters
nam = names(clusters[order(clusters[,1]),])
annotation_col = as.matrix(c(rep("#f15f48", ncol(mg1)),rep("#4ebceb", ncol(mg2)),rep("#8a3795", ncol(mg3))))
rownames(annotation_col) = nam
annotation_col = c(annotation_col[annotation_col == "#f15f48",],annotation_col[annotation_col == "#4ebceb",],annotation_col[annotation_col == "#8a3795",])
double.col = cbind(annotation_col, color.ra)

cbind(names(annotation_col) == names(color.ra))
ann=matrix(nrow=length(annotation_col), ncol = 2)
rownames(ann) = names(annotation_col)
for (name in names(annotation_col)){
  ann[name,] = cbind(annotation_col[name],color.ra[name])
}

# use http://colorbrewer2.org/ to find optimal divergent color palette (or set own)
col.pal = colorRampPalette(c("black","purple","gold1"))
color.palette  <- col.pal(length(palette.breaks) - 1)

heatmap.3(m3_ref[genes,nam], Rowv = T, Colv = F, col=color.palette, breaks = palette.breaks,
          scale="none",trace = "none",density.info = "none", ColSideColors=ann)

dev.off()

#Plot all barcode inf clusters over the tissue
ann.col.inf = matrix(nrow=ncol(RA.norm), ncol = 1) # this object goes for plotting tissue overlays in the plot.2d function
rownames(ann.col.inf) = colnames(RA.norm)
for (sam in colnames(RA.norm)){
  if ((!sam %in% row.names(ann)) == TRUE) ann.col.inf[sam,1] = 4
    else {
      if (ann[sam,1] == "#f15f48") { ann.col.inf[sam,1] = 3}
      if (ann[sam,1] == "#4ebceb") { ann.col.inf[sam,1] = 2}
      if (ann[sam,1] == "#8a3795") { ann.col.inf[sam,1] = 1}
    }
}

# Take only top variable genes using the scran package
all.exp.values.norm = RA.norm
fit = trendVar(all.exp.values.norm)
decomp = decomposeVar(all.exp.values.norm, fit)
top.hvgs = order(decomp$bio, decreasing=TRUE)
decomp_bio = decomp[top.hvgs,][1:1500,]

# Subset on only variable genes
all.exp.values.norm = all.exp.values.norm[rownames(all.exp.values.norm) %in% rownames(decomp_bio),]
all.exp.values.norm = genes_cleanup(all.exp.values.norm)

# Color points based on only subsetting inflitrates ### This is a check point
tsne_inf = all.exp.values.norm[,colnames(all.exp.values.norm) %in% all.inf$new.barcode.names]
tsne_rest = all.exp.values.norm[,!colnames(all.exp.values.norm) %in% all.inf$new.barcode.names]
all_matrix = cbind(tsne_inf, tsne_rest) # this matrix is used in DE

# How many components are important? 
per_out = permutationPA(as.matrix(t(all_matrix)), B = 1, threshold = 0.01, verbose = TRUE, seed = NULL)

# Run tSNE again
tsne_out <- Rtsne(as.matrix(t(all_matrix)), dims = 3, pca = TRUE, initial_dims = per_out$r, theta = 0.5, check_duplicates = FALSE, perplexity = 60, max_iter = 10000, verbose = TRUE) # Run TSNE

# Make matrix of tsne data output
tsne <- as.data.frame(tsne_out$Y)
t1 <- as.matrix(tsne[1])
t2 <- as.matrix(tsne[2])
t3 <- as.matrix(tsne[3])
tsne_m <- cbind (t1,t2,t3)
rownames(tsne_m) = colnames(all_matrix)

# Hierachical clustering
hc = hclust(dist(data.frame(tsne_m[,1],tsne_m[,2])),method = "ward.D2")
clusters = as.matrix(cutree(hc, k = 4)) # 4 clusters
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

#Do DE analysis
m3_ref = all_matrix
mg1 <- m3_ref[,clusters == 1]
mg2 <- m3_ref[,clusters == 2]
mg3 <- m3_ref[,clusters == 3]
mg4 <- m3_ref[,clusters == 4]

# Do DE in one-vs-all
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

### Assing colors to clusters and annotation values
annotation_col = c(annotation_col[annotation_col == "#FED8B1",],annotation_col[annotation_col == "gray61",],annotation_col[annotation_col == "gray76",],annotation_col[annotation_col == "gray85",])
inf_colors = c(rep("#F49ACC", ncol(tsne_inf)), rep("grey51", ncol(tsne_rest)))
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

heatmap.3(all_matrix[genes,nam], Rowv = T, Colv = F, col=color.palette, breaks = palette.breaks,
          scale="none",trace = "none",density.info = "none", cexRow = 0.4, ColSideColors = ann, labCol = 0)
dev.off()

# Which barcodes are part of the annotated infiltrates?
bar_clus1 = mg1[genes,colnames(mg1) %in% row.names(all.inf)]
ncol(bar_clus1)/nrow(all.inf) # percentage of total annotated inflitrates
bar_clus2 = mg2[genes,colnames(mg2) %in% row.names(all.inf)]
ncol(bar_clus2)/nrow(all.inf) # percentage of total annotated inflitrates
bar_clus3 = mg3[genes,colnames(mg3) %in% row.names(all.inf)]
ncol(bar_clus3)/nrow(all.inf) # percentage of total annotated inflitrates
bar_clus4 = mg4[genes,colnames(mg4) %in% row.names(all.inf)]
ncol(bar_clus4)/nrow(all.inf) # percentage of total annotated inflitrates
bar_clus5 = mg5[genes,colnames(mg5) %in% row.names(all.inf)]
ncol(bar_clus5)/nrow(all.inf)# percentage of total annotated inflitrates

### Make distance heatmaps of clones parallel to tSNE for Infliltrates only
pdf("Heatmap_tSNE_Inf.pdf")
all_bar = mbind(bar_clus1, bar_clus2, bar_clus3, bar_clus4)

# Assign colors
annotation_col = as.matrix(c(rep("#FED8B1", ncol(bar_clus1)),rep("gray61", ncol(bar_clus2)),rep("gray76", ncol(bar_clus3)),rep("gray85", ncol(bar_clus4))))
names(annotation_col) =c(colnames(bar_clus1),colnames(bar_clus2),colnames(bar_clus3),colnames(bar_clus4))
inf_colors = c(rep("#FF7F00", ncol(tsne_inf)))
names(inf_colors) =c(colnames(tsne_inf))

ann=matrix(nrow=length(annotation_col), ncol = 2)
rownames(ann) = names(annotation_col)
for (name in names(annotation_col)){
  ann[name,] = cbind(annotation_col[name],inf_colors[name])
}

heatmap.3(all_bar, Rowv = T, Colv = T, col=color.palette, breaks = palette.breaks,cexRow = 0.4,
          scale="none",trace = "none",density.info = "none", ColSideColors = ann, labCol = 0)
dev.off()

# Plot average expression of marker genes per cluster
gen_names = c("CCL21","CD6","LTB","MMP3","CD52","MS4A1","CCL19", "FN1")

cs_all_1 = bar_clus1[gen_names,]
rm_cs_all_1 = rowMeans(cs_all_1)
sd_cs_all_1 = apply(cs_all_1, 1, sd)
sd_cs_all_1 = sd_cs_all_1/sqrt(ncol(bar_clus1))

cs_all_2 = bar_clus2[gen_names,]
rm_cs_all_2 = rowMeans(cs_all_2)
sd_cs_all_2 = apply(cs_all_2, 1, sd)
sd_cs_all_2 = sd_cs_all_2/sqrt(ncol(bar_clus2))

cs_all_3 = bar_clus3[gen_names,]
rm_cs_all_3 = rowMeans(cs_all_3)
sd_cs_all_3 = apply(cs_all_3, 1, sd)
sd_cs_all_3 = sd_cs_all_3/sqrt(ncol(bar_clus3))

cs_all_4 = bar_clus4[gen_names,]
rm_cs_all_4 = rowMeans(cs_all_4)
sd_cs_all_4 = apply(cs_all_4, 1, sd)
sd_cs_all_4 = sd_cs_all_4/sqrt(ncol(bar_clus4))

m_both = cbind(rm_cs_all_1, rm_cs_all_2, rm_cs_all_3, rm_cs_all_4)
sd_both = cbind(sd_cs_all_1,sd_cs_all_2, sd_cs_all_3,sd_cs_all_4)

myData = data.frame(as.numeric(m_both), as.numeric(sd_both))
rownames(myData) = c(paste(row.names(cs_all_1), "Cluster1", sep="_"),paste(row.names(cs_all_1), "Cluster2", sep="_"),
                     paste(row.names(cs_all_1), "Cluster3", sep="_"),paste(row.names(cs_all_1), "Cluster4", sep="_"))
limits <- aes(ymax = myData$as.numeric.m_both. + myData$as.numeric.sd_both.,
              ymin = myData$as.numeric.m_both. - myData$as.numeric.sd_both.)
colvar = rep(c("#FED8B1", "gray61", "gray76", "gray85"), length(gen_names))
names(colvar) = rownames(myData)

pdf("Avg_genes.pdf")
p <- ggplot(data = myData, aes(x = rownames(myData), y = myData$as.numeric.m_both.))
p + ylab(c("Avg Log10 expression"))+ ylim(-3, 10)+ theme_bw() + geom_bar(stat = "identity", position = position_dodge(width = 0.9), fill = colvar) + geom_errorbar(mapping = limits, stat = "identity", position = position_dodge(width = 0.9), width = 0.3)+ theme(axis.line.y = element_line(colour = "black"), axis.title.x = element_blank(), axis.ticks.x= element_blank(), axis.title.y = element_text(), axis.text.x = element_text(angle = 30),panel.border = element_blank(), panel.grid.major = element_blank(),                                                                                                         panel.grid.minor = element_blank())
dev.off()

#Assign colors
ann.col = matrix(nrow=nrow(ann_to_plot), ncol = 1)
unique(ann_to_plot[,1])
unique(ann_to_plot[,2])
rownames(ann.col) = rownames(ann_to_plot)
for (sam in rownames(ann)){
  if ((ann[sam,1] == "#FED8B1") & (ann[sam,2] == "#FF7F00")){ ann.col[sam,1] = "red"}
  else {ann.col[sam,1] = ann[sam,1] }
}

ann.cluster = matrix(nrow=nrow(ann_to_plot), ncol = 1)
rownames(ann.cluster) = rownames(ann_to_plot)
for (sam in rownames(ann_to_plot)){
  if ((ann_to_plot[sam,1] == "#FED8B1") & (ann_to_plot[sam,2] == "#F49ACC") | (ann_to_plot[sam,1] == "gray61") & (ann_to_plot[sam,2] == "#F49ACC") | (ann_to_plot[sam,1] == "gray76") & (ann_to_plot[sam,2] == "#F49ACC") | (ann_to_plot[sam,1] == "gray85") & (ann_to_plot[sam,2] == "#F49ACC")){ ann.cluster[sam,1] = 5}
  else {
    if (ann_to_plot[sam,1] == "#FED8B1") { ann.cluster[sam,1] = 4}
    if (ann_to_plot[sam,1] == "gray61") { ann.cluster[sam,1] = 3}
    if (ann_to_plot[sam,1] == "gray76") { ann.cluster[sam,1] = 2}
    if (ann_to_plot[sam,1] == "gray85") { ann.cluster[sam,1] = 1}}
}

#Plot 3D heatmaps
#RA2
plot.gene.2d.inf.9("RA2", ann.col.inf, m1, m2, m3, m4, m5, m6, m7, s4, s5, s6, s7, s8, s10, s11, x=40, y=20, transparency=1, min=1, max=4)
plot.gene.2d.cluster.9("RA2", ann.cluster, m1, m2, m3, m4, m5, m6, m7, s4, s5, s6, s7, s8, s10, s11, x=40, y=20, transparency=1, min=1, max=5)
plot.gene.2d.9("RA2_Con", "FOSB", m1, m2, m3, m4, m5, m6, m7, s4, s5, s6, s7, s8, s10, s11, x=40, y=20, transparency=1, min=-2, max=6, con = T)
