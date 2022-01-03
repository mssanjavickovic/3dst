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
suppressMessages(suppressWarnings(library(ggpubr,warn.conflicts = F, quietly = T)))
suppressMessages(suppressWarnings(library(ggsignif,warn.conflicts = F, quietly = T)))


# set wd
setwd('/Users/sanjavickovic/Desktop/morphoSPOT/3dst_repo/3dst/de_analysis')

# Which samples do you want to use? 
norm_samples = "RA2"

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
  assign(paste0("inf_all_", i, "_rm"), rm_section_all(get(paste0("RA.norm.inf", i))[gen_names,]))
  assign(paste0("sd", i, "_rm"), sd_section(get(paste0("RA.norm.inf", i))[gen_names,]))
  assign(paste0("s", i, "_rm"), size_section(get(paste0("RA.norm.inf", i))[gen_names,]))
  
}

# Plot avg gene expression for each Inf separately (barplots)
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

# Plot avg gene expression for each Inf separately (barplots)
myplots <- vector("list", length(unique(inf_all$Inf.))*length(unique(sapply(strsplit(row.names(inf_all), "_"), "[[", 1))))
counter = 1
myData_collect = matrix(nrow = 1, ncol = 4)
colnames(myData_collect) = c("Avg", "Section", "Gene", "inf")
for(i in 1:length(unique(inf_all$Inf.))){
  myData = data.frame(get(paste0("inf_all_", i, "_rm")))
  myData = t(myData)
  row.names(myData) = str_replace(sapply(strsplit(row.names(myData), "\\."), "[[",1), "X", "Section")
  myData_0 = matrix(ncol = 3, nrow = 1)
  for (c in 1:ncol(myData)){
    myData_0 = rbind(myData_0, cbind(as.matrix(myData[,c]), as.matrix(row.names(myData)), as.matrix(rep(colnames(myData)[c], nrow(myData)))))
  }
  myData_0 = myData_0[-1,]
  colnames(myData_0) = c("Avg", "Section", "Gene")
  myData_0 = as.data.frame(myData_0)
  myData_0$Avg = as.numeric(myData_0$Avg)
  row.names(myData_0) = seq(1, nrow(myData_0))
  myData_0$inf = rep(paste0("inf", i), nrow(myData_0))
  myData_collect = rbind(myData_collect, myData_0)
  
  for (sec in 1:length(unique(myData_0$Section))){
    myData_sec = myData_0[myData_0$Section == unique(myData_0$Section)[sec],]
    myplots[[counter]] <- ggplot(data = myData_sec, aes(x = Gene, y = Avg, fill = Gene)) + 
      geom_boxplot() + 
      geom_jitter(shape=16, position=position_jitter(0.01)) +
      labs(x = "", y = "Avg expression") + 
      scale_fill_manual(name = "Gene", values=c("#066799",  "#7392CB","#CDCC63")) +
      labs(title = paste0("Infiltrate#", i, " Section#", sec))
    counter = counter + 1

  }
}
myData_collect = myData_collect[-1,]
save_plot(paste0(path_output,"Average_genes_per_infiltrate_boxplots_", norm_samples, ".pdf"), plot_grid(plotlist=myplots, ncol = length(unique(myData_0$Section)), nrow = length(unique(inf_all$Inf.))), base_width = 15, base_height =10)

# export data for manuscript figshare
row.names(myData_collect) = seq(1, nrow(myData_collect))
#write.table(myData_collect, file = "/Users/sanjavickovic/Desktop/morphoSPOT/Manuscript_Nov2021/figshare/SuppFig5b.csv", sep = ",", quote = F)

# Run analysis on infiltrates only
#k = 2 #RA1
k = 3 #RA2
#k = 3 #RA3
#k = 3 #RA4
#k = 2 #RA5
#k = 3 #RA6
col.inf.clusters = run_inf_analysis(RA.norm[,colnames(RA.norm) %in% all_inf], k=k, norm_samples, path_output) # number of clusters

# Run spatial clustering analysis on all ST spots per biopsy
#f = 4 #RA1
f = 3 #RA2
#f = 3 #RA3
#f = 4 #RA4
#f = 4 # RA5
#f = 4 # RA6
col.spatial.clusters = run_spatial_cluster_analysis(RA.norm, f, read_tsne_from_memory="yes", norm_samples, path_output)

# Which barcodes and clusters are part of the annotated infiltrates? 
inf_per_cluster(RA.norm, all_ann_inf, col.spatial.clusters)

#### Make barplot of avg expression per interesting marker genes
gen_names = c("CCL19","CXCL13","LTB","PRG4","MMP3", "CD52","MS4A1","FN1","TYROBP")
save_plot(paste0(path_output,"Average_genes_per_cluster_", norm_samples, ".pdf"), avg_genes_barplot(RA.norm, gen_names, col.spatial.clusters), ncol = 1, nrow = 1, base_width = 10, base_height = 5)
save_plot(paste0(path_output,"Average_genes_per_cluster_boxplot_", norm_samples, ".pdf"), avg_genes_box_old(RA.norm, gen_names, col.spatial.clusters), ncol = 1, nrow = 1, base_width = 10, base_height = 10)
write.table(avg_genes_box_old_data(RA.norm, gen_names, col.spatial.clusters), file = paste0("/Users/sanjavickovic/Desktop/morphoSPOT/Manuscript_Nov2021/figshare/Average_genes_per_cluster_boxplot_", norm_samples, ".csv"), sep = ",", quote = F)

#Plot 3D genes as a heatmap (not morphological)
#Assign color annotations for plotting in 3D for inf and clusters
ann.col.inf = assign_col_inf_numbers(RA.norm, col.inf.clusters) # presumes max 4 clusters
ann.cluster = assign_col_cluster_numbers(col.spatial.clusters, all_ann_inf) # presumes max 3 clusters

# Plot spatial heatmaps per cell type
setwd(path_output)

if (norm_samples == 'RA1' | norm_samples =='RA3' | norm_samples =='RA6') {
  plot.gene.2d.inf.4(norm_samples, ann.col.inf, m1, m2, m3, m4, s1, s2, s3, s4, x=40, y=20, transparency=1, min=1, max=3)
  plot.gene.2d.cluster.4(norm_samples, ann.cluster, m1, m2, m3, m4, s1, s2, s3, s4, x=40, y=20, transparency=1, min=1, max=3)
  plot.gene.2d.4(norm_samples, "CCR7", m1, m2, m3, m4, s1, s2, s3, s4, x=100, y=50, transparency=1, min=0, max=2, con = T)
}

if (norm_samples == 'RA2'){
  plot.gene.2d.inf.7(norm_samples, ann.col.inf, m1, m2, m3, m4, m5, m6, m7, s1, s2, s3, s4, s5, s6, s7, x=40, y=20, transparency=1, min=1, max=4)
  plot.gene.2d.cluster.7(norm_samples, ann.cluster, m1, m2, m3, m4, m5, m6, m7, s1, s2, s3, s4, s5, s6, s7, x=40, y=20, transparency=1, min=1, max=5)
  plot.gene.2d.7(norm_samples, "CD74", m1, m2, m3, m4, m5, m6, m7, s1, s2, s3, s4, s5, s6, s7, x=40, y=40, transparency=1, min=4, max=6, con = T)
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


