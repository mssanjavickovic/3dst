### Bulk level analysis ###
# Loading libraries
library(devtools)
library(scran)
library(Rtsne)
library(jackstraw)
library(gdata)

# Load data as R objects for RA1
setwd("/Users/svickovi/Library/Mobile Documents/com~apple~CloudDocs/Desktop/morphoSPOT/3dst_repo/3dst/data")
load('RA1_norm.exp.values.1') # in ./data
load('RA1_norm.exp.values.2') # in ./data
load('RA1_norm.exp.values.3') # in ./data
load('RA1_norm.exp.values.4') # in ./data

# Combine into avg expression profile per section
r2b1 <- rowMeans(m1)
r2b2 <- rowMeans(m2)
r2b3 <- rowMeans(m3)
r2b4 <- rowMeans(m4)

# # Bind into one matrix
r2 = cbind(r2b1, r2b2, r2b3, r2b4)

# Do same type of analysis as for spatial data
all.exp.values.norm = r2
fit = trendVar(all.exp.values.norm)
decomp = decomposeVar(all.exp.values.norm, fit)
top.hvgs = order(decomp$bio, decreasing=TRUE)
decomp_bio = decomp[top.hvgs,][1:500,]

# Subset on only variable genes
all.exp.values.norm = all.exp.values.norm[rownames(all.exp.values.norm) %in% rownames(decomp_bio),]
all.exp.values.norm = genes_cleanup(all.exp.values.norm)
genes_r2b = rownames(all.exp.values.norm)

# Hierachical clustering
hc = hclust(dist(t(all.exp.values.norm)), method = "ward.D2")
clusters = as.matrix(cutree(hc, k = 2))
row.names(clusters) = colnames(all.exp.values.norm)
col.clusters = c(rep("#FED8B1",length(clusters[clusters == 1])),rep("gray61",length(clusters[clusters == 2])),
                 rep("gray76",length(clusters[clusters == 3])),rep("gray90",length(clusters[clusters == 4])),rep("black",length(clusters[clusters == 5])))
plot(hc, labels = clusters, cex = 0.5)

#Do DEG analysis
m3_ref = all.exp.values.norm
mg1 <- m3_ref[,clusters == 1]
mg2 <- m3_ref[,clusters == 2]

degg <- bimod.diffExp.test(mg1, as.matrix(m3_ref[,!colnames(m3_ref) %in% colnames(mg1)]), row.names(mg1))
clus1.genes <- rownames(subset(degg, degg[,1] < 0.001 & (degg[,2] > 0.5 | degg[,2] < -0.5)))

### Make distance heatmaps of clones parallel to tSNE
genes = unique(c(clus1.genes))
house.keeping.genes = read.delim("house_keeping_final.txt", row.names = 1)
genes <- intersect(genes, outersect(genes, row.names(house.keeping.genes))) # no unique DE genes found between avg exp per section

# Bulk data analysis
# Load data as R objects for RA2
load('RA2_norm.exp.values.1') # in ./data
load('RA2_norm.exp.values.2') # in ./data
load('RA2_norm.exp.values.3') # in ./data
load('RA2_norm.exp.values.4') # in ./data
load('RA2_norm.exp.values.5') # in ./data
load('RA2_norm.exp.values.6') # in ./data
load('RA2_norm.exp.values.7') # in ./data

# Combine into avg expression profile per section
r1b1 <- rowMeans(m1)
r1b2 <- rowMeans(m2)
r1b3 <- rowMeans(m3)
r1b4 <- rowMeans(m4)
r1b5 <- rowMeans(m5)
r1b6 <- rowMeans(m6)
r1b7 <- rowMeans(m7)

# Bind into one matrix
r1 = cbind(r1b1, r1b2, r1b3, r1b4, r1b5, r1b6, r1b7)

# Do same type of analysis as for spatial data
all.exp.values.norm = r2
fit = trendVar(all.exp.values.norm)
decomp = decomposeVar(all.exp.values.norm, fit)
top.hvgs = order(decomp$bio, decreasing=TRUE)
decomp_bio = decomp[top.hvgs,][1:500,]

# Subset on only variable genes
all.exp.values.norm = all.exp.values.norm[rownames(all.exp.values.norm) %in% rownames(decomp_bio),]
all.exp.values.norm = genes_cleanup(all.exp.values.norm)
genes_r1b = rownames(all.exp.values.norm)

# Hierachical clustering
hc = hclust(dist(t(all.exp.values.norm)), method = "ward.D2")
clusters = as.matrix(cutree(hc, k = 3))
row.names(clusters) = colnames(all.exp.values.norm)
col.clusters = c(rep("#FED8B1",length(clusters[clusters == 1])),rep("gray61",length(clusters[clusters == 2])),
                 rep("gray76",length(clusters[clusters == 3])),rep("gray90",length(clusters[clusters == 4])),rep("black",length(clusters[clusters == 5])))
plot(hc, labels = clusters, cex = 0.5)

#Do DEG analysis
m3_ref = all_matrix
mg1 <- m3_ref[,clusters == 1]
mg2 <- m3_ref[,clusters == 2]
mg3 <- m3_ref[,clusters == 3]

degg <- bimod.diffExp.test(mg1, m3_ref[,!colnames(m3_ref) %in% colnames(mg1)], row.names(mg1))
clus1.genes <- rownames(subset(degg, degg[,1] < 0.001 & (degg[,2] > 0.5 | degg[,2] < -0.5)))

degg <- bimod.diffExp.test(mg2,m3_ref[,!colnames(m3_ref) %in% colnames(mg2)], row.names(mg2))
clus2.genes <- rownames(subset(degg, degg[,1] < 0.001 & (degg[,2] > 0.5 | degg[,2] < -0.5)))

degg <- bimod.diffExp.test(mg3,m3_ref[,!colnames(m3_ref) %in% colnames(mg3)], row.names(mg3))
clus3.genes <- rownames(subset(degg, degg[,1] < 0.001 & (degg[,2] > 0.5 | degg[,2] < -0.5)))

### List DE genes
genes = unique(c(clus1.genes,clus2.genes,clus3.genes))
house.keeping.genes = read.delim("house_keeping_final.txt", row.names = 1)
genes <- intersect(genes, outersect(genes,row.names(house.keeping.genes))) # no unique DE genes found between avg exp per section

# Do same type of analysis as for spatial data
exp.values = mbind(r1, r2)
exp.values = exp.values[complete.cases(exp.values),]

# Do same type of analysis as for spatial data
all.exp.values.norm = exp.values
fit = trendVar(all.exp.values.norm)
decomp = decomposeVar(all.exp.values.norm, fit)
top.hvgs = order(decomp$bio, decreasing=TRUE)
decomp_bio = decomp[top.hvgs,][1:500,]

# Subset on only variable genes
all.exp.values.norm = all.exp.values.norm[rownames(all.exp.values.norm) %in% rownames(decomp_bio),]
all.exp.values.norm = genes_cleanup(all.exp.values.norm)
genes_rb_all = rownames(all.exp.values.norm)

# Hierachical clustering
hc = hclust(dist(t(all.exp.values.norm)), method = "ward.D2")
clusters = as.matrix(cutree(hc, k = 2))
row.names(clusters) = colnames(all.exp.values.norm)
col.clusters = c(rep("#FED8B1",length(clusters[clusters == 1])),rep("gray61",length(clusters[clusters == 2])),
                 rep("gray76",length(clusters[clusters == 3])),rep("gray90",length(clusters[clusters == 4])),rep("black",length(clusters[clusters == 5])))
plot(hc, labels = clusters, cex = 0.5)
dev.off()

# PCA
pdf("PCA_RA1_vs_RA2.pdf", useDingbats = F) # this one went in the paper
pc <- prcomp(all.exp.values.norm)
rot <- pc[2]
rot = as.data.frame(rot)
par(mar=c(5.1, 4.1, 4.1, 9.1), xpd=TRUE)
plot(rot$rotation.PC1, rot$rotation.PC2)
cf = c("#E6E5F3", "#EE99A7", "#F04F55", "#ED2224", "#AED8E6", "#98BDDA", "#81A2CE", "#6D88C2", "#586DB3", "#4354A5", "#323E99")
plot(rot$rotation.PC1, rot$rotation.PC2, pch=16,cex=3,col=cf)
legend("topright",inset=c(0,0),row.names(rot), pch=16, xpd = TRUE, horiz = FALSE, bty = "n", cex = 1,col=cf)
dev.off()

summary(pc)
eigs <- pc$sdev^2
eigs[1] / sum(eigs)
eigs[2] / sum(eigs)
eigs[3] / sum(eigs)

#Do DEG analysis
m3_ref = all.exp.values.norm
mg1 <- m3_ref[,clusters == 1]
mg2 <- m3_ref[,clusters == 2]

degg <- bimod.diffExp.test(mg1, m3_ref[,!colnames(m3_ref) %in% colnames(mg1)], row.names(mg1))
clus1.genes <- rownames(subset(degg, degg[,1] > 0.001 & (degg[,2] > 1 | degg[,2] < -1)))
clus1.genes.values = degg[clus1.genes,]

degg <- bimod.diffExp.test(mg2,m3_ref[,!colnames(m3_ref) %in% colnames(mg2)], row.names(mg2))
clus2.genes <- rownames(subset(degg, degg[,1] < 0.001 & (degg[,2] > 1 | degg[,2] < -1)))
clus2.genes.values = degg[clus2.genes,]

### List of de genes
genes = unique(c(clus1.genes,clus2.genes))
house.keeping.genes = row.names(read.delim("house_keeping_final.txt", row.names = 1))
genes <- intersect(genes, outersect(genes,row.names(house.keeping.genes)))

# Plot heatmap of de genes
pdf("Heatmap2_RA_Bulk.pdf")

source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")

annotation_col=cf
names(annotation_col) = c(colnames(all.exp.values.norm))
annotation_col = as.matrix(annotation_col)
quantile.range <- quantile(m3_ref, probs = seq(0, 1, 0.05))
palette.breaks <- seq(quantile.range["60%"], quantile.range["100%"], 0.1)

# use http://colorbrewer2.org/ to find optimal divergent color palette (or set own)
col.pal = colorRampPalette(c("black","purple","darkorchid1","gold1"))
color.palette  <- col.pal(length(palette.breaks) - 1)

heatmap.3(all.exp.values.norm[genes,], Rowv = T, Colv = F, col=color.palette, breaks = palette.breaks,
          scale="none",trace = "none",density.info = "none", ColSideColors = annotation_col, labCol = 0, cexRow = 0.3)

dev.off()
