### DE analysis for RA 3dst ### 

# removes all glob variables from env
rm(list = ls())

# Loading libraries
suppressMessages(suppressWarnings(library(devtools,warn.conflicts = F, quietly = T)))
suppressMessages(suppressWarnings(library(gplots,warn.conflicts = F, quietly = T)))
suppressMessages(suppressWarnings(library(destiny,warn.conflicts = F, quietly = T)))
suppressMessages(suppressWarnings(library(dpt,warn.conflicts = F, quietly = T)))
suppressMessages(suppressWarnings(library(scran,warn.conflicts = F, quietly = T)))
suppressMessages(suppressWarnings(library(plot3D,warn.conflicts = F, quietly = T)))
suppressMessages(suppressWarnings(library(Rtsne,warn.conflicts = F, quietly = T)))
suppressMessages(suppressWarnings(library(jackstraw,warn.conflicts = F, quietly = T)))
suppressMessages(suppressWarnings(library(rgl,warn.conflicts = F, quietly = T)))
suppressMessages(suppressWarnings(library(gdata,warn.conflicts = F, quietly = T)))
suppressMessages(suppressWarnings(library(ggplot2,warn.conflicts = F, quietly = T)))
suppressMessages(suppressWarnings(library(akima,warn.conflicts = F, quietly = T)))
suppressMessages(suppressWarnings(library(matrixStats,warn.conflicts = F, quietly = T)))
suppressMessages(suppressWarnings(library(cowplot,warn.conflicts = F, quietly = T)))
suppressMessages(suppressWarnings(library(pheatmap,warn.conflicts = F, quietly = T)))
suppressMessages(suppressWarnings(library(stringr,warn.conflicts = F, quietly = T)))
suppressMessages(suppressWarnings(library(grid,warn.conflicts = F, quietly = T)))

setwd('/Users/svickovi/Library/Mobile Documents/com~apple~CloudDocs/Desktop/morphoSPOT/3dst_repo/3dst/de_analysis')
#setwd('/Users/sanjavickovic/Desktop/morphoSPOT/3dst_repo/3dst/de_analysis')

# Which samples do you want to use? 
norm_samples = "RA5"

# Where are you norm expression R files located? 
path_samples = "../data/"

# Set output directory that will contain all cell_typing pdf plots and output gene files
path_output = "../../../plot_outputs/"

# Source f(x) file
source('../functions/Read_Functions.R') # in ./functions

# Load data as R objects for selected RA patient
files_norm = list.files(pattern = glob2rx(paste0(norm_samples, "_norm.exp.values.*")), path = path_samples)
norm_mat_loaded_in_env = ""
section_labels = ""
all_ann_inf = ""
for (i in 1:length(files_norm)){
  load(paste0(path_samples, files_norm[i]))
  load(paste0(path_samples, norm_samples, "_", i, "_selected_adjusted_spots_3D_manual_app"))
  assign(paste0("inf_s", i), read.csv(paste0("../data/", norm_samples, "_", i, "_all_inf.csv"), header = T, row.names = 1)) # files found on SCP
  norm_mat_loaded_in_env = c(norm_mat_loaded_in_env, paste0("m", i) )
  section_labels = c(section_labels, sapply(strsplit(colnames(get(paste0("m", i))), "_"), "[[", 1))
  all_ann_inf = c(all_ann_inf, row.names(get(paste0("inf_s", i))))
}
norm_mat_loaded_in_env = norm_mat_loaded_in_env[-1]
section_labels = section_labels[-1]
all_ann_inf = all_ann_inf[-1]

# Combine all the norm data in one matrix
RA.norm = do.call(mbind, mget(norm_mat_loaded_in_env))
RA.norm[is.na(RA.norm)] <- 0
m = RA.norm # this is a combined matrix for all sections from a single RA patient # precaution

### Take only infiltrates that are present in all sections according to annotations
inf_all = read.csv(paste0("../data/", norm_samples, "_zstack_Infs.csv"), header = T)

# Subset per infiltrate as the naming ie. Inf1 presumes we are following the same infiltrate (given 3D alignment) in all present sections
all_inf = ""
for (i in 1:length(unique(inf_all$Inf.))){
  assign(paste0("RA.norm.inf", i), as.matrix(RA.norm[,colnames(RA.norm) %in% row.names(subset(inf_all, inf_all$Inf. == paste0("Inf", i)))]))
  all_inf = c(all_inf, colnames(get(paste0("RA.norm.inf", i)))) # Make an object with all annotated infiltrates
}
all_inf = all_inf[-1]

### Make avg gene exp barplot per inf and section
gen_names = c("CD52","MS4A1","FN1") # list of marker genes

# Calculate avg gene expression per inf and section
for (i in 1:length(unique(inf_all$Inf.))){
  assign(paste0("inf", i, "_rm"), rm_section(get(paste0("RA.norm.inf", i))[gen_names,]))
  assign(paste0("sd", i, "_rm"), sd_section(get(paste0("RA.norm.inf", i))[gen_names,]))
  assign(paste0("s", i, "_rm"), size_section(get(paste0("RA.norm.inf", i))[gen_names,]))
}

# Plot avg gene expression for each Inf separately
myplots <- vector("list", length(unique(inf_all$Inf.)))
for(i in 1:length(unique(inf_all$Inf.))){
  myData = data.frame(as.numeric(get(paste0("inf", i, "_rm"))),as.numeric(sweep(get(paste0("sd", i, "_rm")), 2, sqrt(get(paste0("s", i, "_rm"))), "/")))
  colnames(myData) = c("Avg", "SEM")
  nms = ""
  for (j in 1:length(unique(section_labels))){
    nms = c(nms, paste0("Section",j, "_", gen_names))
  }
  nms = nms[-1]
  row.names(myData) = nms
  limits <- aes(ymax = Avg + SEM, ymin = Avg - SEM)
  
  myplots[[i]] <- ggplot(data = myData, aes(x = factor(sapply(strsplit(row.names(myData), split = "_"),"[[",1)), y = Avg, fill = rep(c("Gene1: CD52","Gene2: MS4A1","Gene3: FN1"), length(unique(section_labels))))) + geom_bar(stat = "identity", position = position_dodge(0.99)) + 
    geom_errorbar(limits, position = position_dodge(0.99), width = 0.25) + 
    labs(x = "", y = "Avg expression") + 
    scale_fill_manual(name = "Gene", values=c("#066799",  "#7392CB","#CDCC63")) +
    labs(title = paste0("Infiltrate ", i))
}
save_plot(paste0(path_output,"Average_genes_per_infiltrate_", norm_samples, ".pdf"), plot_grid(plotlist=myplots), ncol = 3, nrow = 2)

# Run analysis on infiltrates only
#k = 2 #RA1
#k = 3 #RA2
#k = 3 #RA3
#k = 3 #RA4
#k = 2 #RA5
col.inf.clusters = run_inf_analysis(RA.norm[,colnames(RA.norm) %in% all_inf], k=2, norm_samples, path_output) # number of clusters

# Run spatial clustering analysis on all ST spots per biopsy
#f = 4 #RA1
#f = 3 #RA2
#f = 3 #RA3
#f = 4 #RA4
#f = 4 # RA5
col.spatial.clusters = run_spatial_cluster_analysis(RA.norm, 4, read_tsne_from_memory="yes", norm_samples, path_output)

# Which barcodes and clusters are part of the annotated infiltrates? 
inf_per_cluster(RA.norm, all_ann_inf, col.spatial.clusters)

#### Make barplot of avg expression per interesting marker genes
gen_names = c("CCL19","CXCL13","LTB","PRG4","MMP3","CD52","MS4A1","FN1", "TYROBP")
save_plot(paste0(path_output,"Average_genes_per_cluster_", norm_samples, ".pdf"), avg_genes_barplot(RA.norm, gen_names, col.spatial.clusters), ncol = 1, nrow = 1)

### Make gene-to-gene correlation plots 
# df1 = data.frame(t(RA.norm[c("CD52", "MS4A1"),]))
# df1[,3] = col.spatial.clusters[row.names(df1),]
#df2 = data.frame(t(RA.norm[c("CD52", "MS4A1"),]))
#df2[,3] = col.spatial.clusters[row.names(df2),]
#df3 = data.frame(t(RA.norm[c("CD52", "MS4A1"),]))
#df3[,3] = col.spatial.clusters[row.names(df3),]
# df4 = data.frame(t(RA.norm[c("CD52", "MS4A1"),]))
# df4[,3] = col.spatial.clusters[row.names(df4),]
#df5 = data.frame(t(RA.norm[c("CD52", "MS4A1"),]))
#df5[,3] = col.spatial.clusters[row.names(df5),]
#t1 = c(mean(df1[df1[,1] & df1[,3] == "#FED8B1",][,1]),mean(df2[df2[,1] & df2[,3] == "#FED8B1",][,1]),mean(df3[df3[,1] & df3[,3] == "#FED8B1",][,1]))
#t2 = c(mean(df4[df4[,1] & df4[,3] == "#FED8B1",][,1]),mean(df5[df5[,1] & df5[,3] == "#FED8B1",][,1]))
#t.test(t1, t2, paired = FALSE, alternative = "greater")$p.value

# Change the confidence interval fill color
ggplot(df, aes(x=df[,1], y=df[,2], color = df[,3])) + 
  geom_point(shape=18, color = df[,3])

nrow(df[df[,1] > 0.5 & df[,2] < 0.5,])/nrow(df) # in the whole data
nrow(df[df[,1] < 0.5 & df[,2] > 0.5 & df[,3] == "#FED8B1",])/nrow(df[df[,3] == "#FED8B1",]) 

#Plot all barcode clusters over the tissue for Fig3a
# ann.col = matrix(nrow=nrow(ann_to_plot), ncol = 1)
# rownames(ann.col) = rownames(ann_to_plot)
# for (sam in rownames(ann)){
#   if ((ann[sam,1] == "#FED8B1") & (ann[sam,2] == "#FF7F00")){ ann.col[sam,1] = "red"}
#   else {ann.col[sam,1] = ann[sam,1] }
# }

#Plot 3D genes as a heatmap (not morphological)
#Assign color annotations for plotting in 3D for inf and clusters
ann.col.inf = assign_col_inf_numbers(RA.norm, col.inf.clusters) # presumes max 4 clusters
ann.cluster = assign_col_cluster_numbers(col.spatial.clusters, all_ann_inf) # presumes max 3 clusters

# Plot spatial heatmaps per cell type
setwd(path_output)

if (norm_samples == 'RA1' | norm_samples =='RA3' ) {
  plot.gene.2d.inf.4(norm_samples, ann.col.inf, m1, m2, m3, m4, s1, s2, s3, s4, x=40, y=20, transparency=1, min=1, max=4)
  plot.gene.2d.cluster.4(norm_samples, ann.cluster, m1, m2, m3, m4, s1, s2, s3, s4, x=40, y=20, transparency=1, min=1, max=5)
  plot.gene.2d.4(norm_samples, "CD79A", m1, m2, m3, m4, s1, s2, s3, s4, x=40, y=20, transparency=1, min=0, max=6, con = T)
}

if (norm_samples == 'RA2'){
  plot.gene.2d.inf.7(norm_samples, ann.col.inf, m1, m2, m3, m4, m5, m6, m7, s1, s2, s3, s4, s5, s6, s7, x=40, y=20, transparency=1, min=1, max=4)
  plot.gene.2d.cluster.7(norm_samples, ann.cluster, m1, m2, m3, m4, m5, m6, m7, s1, s2, s3, s4, s5, s6, s7, x=40, y=20, transparency=1, min=1, max=5)
  plot.gene.2d.7(norm_samples, "LTB", m1, m2, m3, m4, m5, m6, m7, s1, s2, s3, s4, s5, s6, s7, x=40, y=40, transparency=1, min=0, max=4, con = T)
} 
if (norm_samples == 'RA4') {
  plot.gene.2d.inf.5(norm_samples, ann.col.inf, m1, m2, m3, m4, m5, s1, s2, s3, s4, s5, x=20, y=40, transparency=1, min=1, max=4)
  plot.gene.2d.cluster.5(norm_samples, ann.cluster, m1, m2, m3, m4, m5, s1, s2, s3, s4, s5, x=30, y=40, transparency=1, min=1, max=5)
  plot.gene.2d.5(norm_samples, "ACTB", m1, m2, m3, m4, m5, s1, s2, s3, s4, s5, x=20, y=40, transparency=1, min=0, max=6, con = T)
}
if (norm_samples == 'RA5') {
  plot.gene.2d.inf.3(norm_samples, ann.col.inf, m1, m2, m3, s1, s2, s3, x=30, y=30, transparency=1, min=1, max=4)
  plot.gene.2d.cluster.3(norm_samples, ann.cluster, m1, m2, m3,  s1, s2, s3, x=30, y=30, transparency=1, min=1, max=5)
  plot.gene.2d.3(norm_samples, "ACTB", m1, m2, m3, s1, s2, s3, x=30, y=30, transparency=1, min=0, max=6, con = T)
}
