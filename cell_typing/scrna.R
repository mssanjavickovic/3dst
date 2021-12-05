#### Single cell associations in ST RA ###
#Load libraries
suppressMessages(suppressWarnings(library(plot3D,warn.conflicts = F, quietly = T)))
suppressMessages(suppressWarnings(library(gdata,warn.conflicts = F, quietly = T)))
suppressMessages(suppressWarnings(library(qdapTools,warn.conflicts = F, quietly = T)))
suppressMessages(suppressWarnings(library(data.table,warn.conflicts = F, quietly = T)))
suppressMessages(suppressWarnings(library(akima,warn.conflicts = F, quietly = T)))
suppressMessages(suppressWarnings(library(stringr,warn.conflicts = F, quietly = T)))
suppressMessages(suppressWarnings(library(ggplot2,warn.conflicts = F, quietly = T)))
suppressMessages(suppressWarnings(library(RColorBrewer,warn.conflicts = F, quietly = T)))
suppressMessages(suppressWarnings(library(cowplot,warn.conflicts = F, quietly = T)))
suppressMessages(suppressWarnings(library(reshape2,warn.conflicts = F, quietly = T)))
suppressMessages(suppressWarnings(library(dplyr,warn.conflicts = F, quietly = T)))
suppressMessages(suppressWarnings(library(ggpubr,warn.conflicts = F, quietly = T)))

setwd("/Users/sanjavickovic/Desktop/morphoSPOT/3dst_repo/3dst/cell_typing")

# removes all glob variables from env
rm(list = ls())

# Run sc imputation for ex. RA1 samples
# Which samples do you want to use? 
norm_samples = "RA2"

# Where are you norm expression R files located? 
path_samples = "../data/"

# Set output directory that will contain all cell_typing pdf plots and output gene files
path_output = "../../../plot_outputs/"

# Load data as R objects 
### loads the normalized R obj and xyz spots coordinates R obj
files_norm = list.files(pattern = glob2rx(paste0(norm_samples, "_norm.exp.values.*")), path = path_samples)
norm_mat_loaded_in_env = ""
for (i in 1:length(files_norm)){
  load(paste0(path_samples, files_norm[i]))
  load(paste0(path_samples, norm_samples, "_", i, "_selected_adjusted_spots_3D_manual_app"))
  norm_mat_loaded_in_env = c(norm_mat_loaded_in_env, paste0("m", i) )
}
norm_mat_loaded_in_env = norm_mat_loaded_in_env[-1]

# Source f(x) file
source("../functions/Read_Functions.R") # in ./functions

# Combine all the norm data in one matrix
RA.norm = do.call(mbind, mget(norm_mat_loaded_in_env))
RA.norm[is.na(RA.norm)] <- 0
m = RA.norm # this is a combined matrix for all sections from a single RA patient # precaution

# Read in sc data
degs = read.csv("../data/Stephenson.clusters.csv", header = T) # This is a table with DE genes avaiable from Stephenson et al
dates = read.delim("../data/Genes.Dates.txt", header = F) # This is a table provided by 10X with all "failed gene name conversions from Excel"
inf_all = read.csv(paste0("../data/", norm_samples, "_zstack_Infs.csv"), header = T)

# Change to output directory
setwd(path_output)

# Clean up gene names in the Stephenson table
gn = degs$gene
j = ""
for (i in gn){
  if (i %in% dates$V2) j = c(j, dates[i,]$V1)
  else j = c(j, i)
}
j = j[-1]
degs$gene <- j

# Select top 100 markers from the cleaned table for each cell type
markers = data.table(degs)
m.use = markers[avg_logFC > 1][p_val_adj < .05][order(-avg_logFC)][,gene[1:500], cluster]

# Convert to named list
m.use = tapply(m.use$V1, m.use$cluster, function(a) sort(unique(a)))

# Then score these in the expression matrix you have, e.g.
scores = sapply(m.use, function(genes.use) colSums(m[row.names(m) %in% genes.use,]))

# Name and output if needed
names(m.use) = paste(paste0("scRNAseq_", tolower(norm_samples)), names(m.use), sep="_")
for (i in 1:length(names(m.use))){
  nm = m.use[i]
  write.table(list2df(nm)$X1, file = paste(names(m.use)[i], "_clusters.txt", sep=""), sep="\t", quote =F)
}

# Choose which cluster files you want to use based on names
cluster_files = list.files('./', pattern=glob2rx(paste0(paste0("scRNAseq_", tolower(norm_samples)), "_*_clusters.txt")))

# Perform scoring and centering
genes_selected = ""
rn_ref = "" # make first gene subset
fil_m = matrix(ncol=ncol(m), nrow=1) # make smaller matrix based on gene subsets
all_clusters = "" # collect all cluster names
co = 0.0 # Set cut off for gene-to-gene correlations, recommended 0.2

# Start the loop
for (i in c(1:length(cluster_files))){
  name = cluster_files[i]
  genes = read.delim(name, header=T)
  
  # clean up lists
  gene_list = genes
  gene_list = gene_list[!is.na(gene_list)]
  if (length(gene_list) == 0) next # if no genes found in this dataset, skip
  sample = sapply(strsplit(name,split=".txt"),"[[",1) # get unique sample name
  
  # Subset the large matric based on genes per cluster; tidy up the matrix
  m_sub = m[row.names(m) %in% gene_list,]
  
  # Do scoring and centring
  if (is.null(nrow(m_sub)) == "TRUE") {
    rn = m_sub
    rn_means <- sum(rn)
    rn_cent <- rn_means/max(rn_means)
    rn_cent = as.matrix(t(rn_cent))
    row.names(rn_cent) <- c(paste(sample,gene_list, sep="_"))
    
  } else {
    rn = m_sub[rowSums(m_sub)>0,] # remove possibility of cor between 0s
    
    # Do gene-to-gene correlations and filtering (within each cluster gene list)
    z <- cor(t(rn),  method = c("pearson"))
    if (is.null(nrow(z)) == TRUE) {next}
    z[lower.tri(z,diag=TRUE)]=NA  # Prepare to drop duplicates and meaningless information
    z=as.data.frame(as.table(z))  # Turn into a 3-column table
    z=na.omit(z)  # Get rid of the junk we flagged above
    z=z[order(-abs(z$Freq)),]  # order based on frequencies
    if (nrow(z) < 2) next # if less than 2 genes are coexpressed, skip
    
    # Subest based on  cut off and save genes per cluster
    n <- (subset(z, z [,3] > co)) [,1:2]
    rn_ref <- c(rn_ref, union(n[,1], n[,2]))
    un = union(n[,1], n[,2])
    if (length(un) == 0) next # if no genes passed filter, skip
    save(un, file = paste("cluster_genes_ra1", i, sep=""))
    genes_selected = c(genes_selected, un)
    write.table(un, file = paste0(sample, "_subsetted.txt"), quote = F, row.names = T) # this is the table used as input to GO term analysis
    
    # Subset based on "good" genes per cluster and center the data
    rn <- m_sub[row.names(m_sub) %in% un,]
    rownames(rn) <- c(paste(sample,rownames(rn), sep="_"))
    rn_means <- colSums(rn)
    rn_cent <- rn_means/max(rn_means)
  }
  
  # Add all "good" and "centered" values in a new matrix
  all_clusters = c(all_clusters, paste(sample)) # collect all cluster names
  fil_m = rbind(fil_m, rn_cent)
}

# Clean up
fil_m = fil_m[-1,]
all_clusters = all_clusters[-1]
row.names(fil_m) = all_clusters
colnames(fil_m) = colnames(m)

# Rename matrix
mf1 = fil_m
mf1 = sweep(mf1, 2, colSums(mf1),`/`)
write.table(mf1, file = paste0("scRNAseq_all_sections_", norm_samples, ".csv"), sep =",", quote =F, row.names = T, col.names = T)

# Split mf1 matrix into smaller matrix objects per section for plotting
for (i in 1:length(files_norm)){
  print(i)
  tmp = mf1[,grep(pattern = paste0("X", i), colnames(mf1), value = F)]
  tmp[is.nan(tmp)] <- 0
  print(min(tmp))
  assign(paste0("mn", i),tmp)
  x = sapply(strsplit(colnames(get(paste0("mn", i))), "_"), "[[",2)
  y = sapply(strsplit(colnames(get(paste0("mn", i))), "_"), "[[",3)
  img = paste0(norm_samples, "_", i)
  tmp1 = cbind(img, x, y, t(get(paste0("mn", i))))
  colnames(tmp1) = c("Image_ID", "x", "y", str_replace(str_replace(row.names(get(paste0("mn", i))), "_clusters", ""), paste0("scRNAseq_", tolower(norm_samples), "_"), ""))
  print(sum(colSums(t(get(paste0("mn", i))))))
  write.table(tmp1, file = paste0("scRNAseq",i,".", norm_samples, ".csv"), sep =",", quote =F, row.names = F, col.names = T)
}



# Plot spatial heatmaps per cell type
setwd(path_output)
# for (i in row.names(mf1)){
#   # Plot spatial spots
#   ## (un)comment for specific samples ex. use only plot.gene.2d.4 for RA1
#   if (norm_samples == 'RA1') plot.gene.2d.4("RA1_scRNAseqtest",i, mn1, mn2, mn3, mn4, s1, s2, s3, s4, con = T, x=40, y=20, transparency=1, min=min(mf1), max=max(mf1)) # for RA1
#   if (norm_samples == 'RA2') plot.gene.2d.7("RA2_scRNAseq",i, mn1, mn2, mn3, mn4, mn5, mn6, mn7, s1, s2, s3, s4, s5, s6, s7, con = T, x=40, y=20, transparency=1, min=min(mf1), max=max(mf1)) # for RA2
#   #if (norm_samples == 'RA2') plot.gene.3d.7("RA2_scRNAseq",i, mn1, mn2, mn3, mn4, mn5, mn6, mn7, s1, s2, s3, s4, s5, s6, s7, x=40, y=20, transparency=0.1, min=min(mf1), max=max(mf1)) # for RA2
#   if (norm_samples == 'RA3') plot.gene.2d.4("RA3_scRNAseq",i, mn1, mn2, mn3, mn4, s1, s2, s3, s4, con = T, x=40, y=40, transparency=1, min=min(mf1), max=0.15) # for RA3
#   if (norm_samples == 'RA4') plot.gene.2d.5("RA4_scRNAseq",i, mn1, mn2, mn3, mn4, mn5, s1, s2, s3, s4, s5, con = T, x=40, y=40, transparency=1, min=min(mf1), max=0.25) # for RA4
#   if (norm_samples == 'RA5') plot.gene.2d.3("RA5_scRNAseq",i, mn1, mn2, mn3, s1, s2, s3, con = T, x=80, y=80, transparency=1, min=min(mf1), max=max(mf1)) # for RA5
#   if (norm_samples == 'RA6') plot.gene.2d.4("RA6_scRNAseq",i, mn1, mn2, mn3, mn4, s1, s2, s3, s4, con = T, x=40, y=40, transparency=1, min=min(mf1), max=0.2) # for RA6
# }

# Avg cell signature per infiltrate 
### Take only infiltrates that are present in all sections according to annotations
setwd("/Users/sanjavickovic/Desktop/morphoSPOT/3dst_repo/3dst/cell_typing")
mf1 = fil_m

# Subset per infiltrate as the naming ie. Inf1 presumes we are following the same infiltrate (given 3D alignment) in all present sections
all_inf = ""
for (i in 1:length(unique(inf_all$Inf.))){
  assign(paste0("SC.norm.inf", i), as.matrix(mf1[,colnames(mf1) %in% row.names(subset(inf_all, inf_all$Inf. == paste0("Inf", i)))]))
  all_inf = c(all_inf, colnames(get(paste0("SC.norm.inf", i)))) # Make an object with all annotated infiltrates
}
all_inf = all_inf[-1]

# get all cluster x_y from Cluster*RA*
cl = list.files(pattern = glob2rx(paste0("Cluster*", norm_samples, "*.csv")), path = path_output)
clusters = matrix(nrow = 1, ncol = 2)
colnames(clusters) = c("spcluster", "section")
for (f in cl){
  scabun = read.csv(paste0(path_output, f), sep = ",")
  tmp = cbind(paste(paste0("X", sapply(strsplit(scabun[,1], "_"), "[[", 2)),scabun[,2],scabun[,3], sep = "_"), sapply(strsplit(f, "\\."), "[[", 1), paste0("Section", sapply(strsplit(scabun[,1], "_"), "[[", 2)))
  row.names(tmp) = tmp[,1]
  tmp = as.matrix(tmp[,-1])
  colnames(tmp) = c("spcluster", "section")
  clusters = rbind(clusters, tmp)
}
clusters = as.matrix(clusters[-1,])
colnames(clusters) = c("spcluster", "section")

# get scores per cluster
dfcl = as.data.frame(merge(t(mf1),clusters, by = 0))
rownames(dfcl) = dfcl$Row.names
dfcl = dfcl[,-1]
dfcl = melt(dfcl)
colnames(dfcl) = c("spcluster", "section", "ct", "value")
dfcl$short_ct = paste0(sapply(strsplit(as.character(dfcl$ct), "_"), "[[",3),
                       "_",sapply(strsplit(as.character(dfcl$ct), "_"), "[[",4))

dfclavg = dfcl %>%
  group_by(spcluster, section, short_ct) %>%
  summarise(meanct = mean(value)) %>% 
  mutate(freq = 100*meanct / sum(meanct))%>%
  slice_max(order_by = meanct, n = 3)

set3 <- colorRampPalette(c(brewer.pal('Set3',n=12),c(brewer.pal('Set1',n=9))))
setnm = setNames(set3(16), unique(dfcl$short_ct))
setnm_final = setnm[names(setnm) %in% unique(dfcl$short_ct)]

#plot scores per cluster
myplots <- vector("list", 1)
plotsizer = length(unique(dfclavg$short_ct))*length(unique(dfclavg$section))

myplots[[1]] <-ggplot(data = dfclavg, aes(x = short_ct, y = freq, fill = short_ct, na.rm = TRUE)) + 
  stat_boxplot_custom(qs = c(0.0, 0.25, 0.5, 0.75, 1),size=.3) + 
  facet_grid(~spcluster) + 
  labs(x = "", y = "Cell type [%]") + 
  scale_fill_manual(name = "Cell type", values = setnm_final)+  
  labs(title = paste0(norm_samples, " cell type scores per cluster")) + theme(axis.text.x = element_text(angle = 30, vjust=1, hjust=1), legend.key.size = unit(0.2, "cm"), axis.title.x=element_blank())#,axis.text.x=element_blank(),axis.ticks.x=element_blank())
save_plot(paste0(path_output,"Average_SC_sigs_per_cluster_", norm_samples, ".pdf"), plot_grid(plotlist=myplots), ncol = 1, nrow = 1)
write.table(dfclavg, file = paste0(path_output,"Average_SC_per_spatial_cluster_", norm_samples, ".csv"), sep = ",", quote = F)

 # Calculate avg gene expression per inf and section
inf_mat_loaded_in_env=""
for (i in 1:length(unique(inf_all$Inf.))){
  assign(paste0("inf", i, "_rm"), sum_section(get(paste0("SC.norm.inf", i))))
  assign(paste0("inf", i, "_rm"), 100*sweep(get(paste0("inf", i, "_rm")), 2, colSums(get(paste0("inf", i, "_rm"))), "/"))
  inf_mat_loaded_in_env = c(inf_mat_loaded_in_env, paste0("inf", i, "_rm"))
}
inf_mat_loaded_in_env = inf_mat_loaded_in_env[-1]
Infs = do.call(mbind, mget(inf_mat_loaded_in_env))
Infs[is.na(Infs)] <- 0
colnames(Infs) = str_replace(colnames(Infs), "X", "Section")
colnames(Infs) = paste0("inf", 1:length(inf_mat_loaded_in_env), "_",colnames(Infs))
Infs = t(Infs)

# Calculate average cell type contribution to the rest of the biopsy (non-inf)
rest.norm = mf1[,!colnames(mf1) %in% row.names(inf_all)]

# Calculate avg gene expression for all spots not in inf (so called rest)
assign("avg_rest", sum_section(rest.norm))
avg_rest = 100*sweep(avg_rest,2,colSums(avg_rest), "/")
colnames(avg_rest) = paste0("rest_Section", 1:length(norm_mat_loaded_in_env))

# plot Inf and rest together per section
myplots <- vector("list", 1)
melted <- melt(rbind(Infs, t(avg_rest)))
colnames(melted) = c("area","ct","Avg")
melted$inf <- sapply(strsplit(as.character(melted$area), "_"),"[[",1)
melted$section <- sapply(strsplit(as.character(melted$area), "_"),"[[",2)
set3 <- colorRampPalette(c(brewer.pal('Set3',n=12),c(brewer.pal('Set1',n=9))))
setnm = setNames(set3(16), sapply(strsplit(cluster_files, ".txt",""), "[[",1))
setnm_final = setnm[names(setnm) %in% unique(melted$ct)]

myplots[[1]] <-ggplot(data = melted, aes(x = inf, y = Avg, fill = ct)) + 
  geom_bar(position="fill", stat="identity") + 
  facet_grid(~section) + 
  labs(x = "", y = "Cell type [%]") + 
  scale_fill_manual(name = "Cell types", values = setnm_final) +
  labs(title = paste0("Cell type scores")) + theme(legend.key.size = unit(0.2, "cm"))
print(myplots[[1]])
save_plot(paste0(path_output,"Average_SC_sigs_per_section_", norm_samples, ".pdf"), plot_grid(plotlist=myplots), ncol = 3, nrow = 2)

melted$cluster = substr(melted$inf, 1,3)
melted$ct_short = paste0(sapply(strsplit(as.character(melted$ct), "_"), "[[",3),
                         "_",sapply(strsplit(as.character(melted$ct), "_"), "[[",4))
counter = 1
for (j in unique(melted$ct_short)){
  print(j)
  meltedsub = melted[melted$ct_short == j,]
  myplots[[counter]] <- ggplot(meltedsub, aes(Avg, y = factor(ct, levels = paste0("celltype_", 1:length(unique(melted$cluster)))), fill = cluster)) +
    ggridges::geom_density_ridges(alpha = 0.5) +
    #scale_x_continuous(limits = c(0, 0.2), labels = scales::percent) +
    labs(y = j) +
           #scale_fill_manual(values = c("cluster 1" = "#4477AA", "background" = "#CC6677")) +
           theme_minimal()
         counter = counter + 1
}
#plot_grid(plotlist=myplots)
save_plot(paste0(path_output,"Average_SC_density_", norm_samples, ".pdf"), plot_grid(plotlist=myplots), ncol = 4, nrow = 4)
melted$ra = norm_samples
write.table(melted, file = paste0(path_output,"Average_SC_density_", norm_samples, ".csv"), sep = ",", quote = F)

# plot Inf and rest together per infiltrate
myplots <- vector("list", 1)
myplots[[1]] <- ggplot(data = melted, aes(x = section, y = Avg, fill = ct)) + 
  geom_bar(position="fill", stat="identity") + 
  facet_grid(~inf) + 
  labs(x = "", y = "Cell type [%]") + 
  scale_fill_manual(name = "Cell types", values = setnm_final) +
  labs(title = paste0("Cell type scores")) + theme(legend.key.size = unit(0.2, "cm")) +
  theme(axis.text.x=element_text(angle=30, hjust=1))
print(myplots[[1]])
save_plot(paste0(path_output,"Average_SC_sigs_per_inf_", norm_samples, ".pdf"), plot_grid(plotlist=myplots), ncol = 3, nrow = 2)

# plot correlation plots between cell type abundances per inf
#ct1 = paste0("scRNAseq_", tolower(norm_samples), "_", "macrophage", "_clusters") #echange cell types here
ct1 = paste0("scRNAseq_", tolower(norm_samples), "_", "fibro2b_thy1", "_clusters") #echange cell types here #b_cells
#ct1 = paste0("scRNAseq_", tolower(norm_samples), "_", "dendritic_cells", "_clusters") #echange cell types here
ct2 = paste0("scRNAseq_", tolower(norm_samples), "_", "plasma_cells", "_clusters") #echange cell types here #b_cells

melted_cts_cor1 = melted[(melted$ct == ct1),] # & melted$inf != 'rest'
melted_cts_cor2 = melted[(melted$ct == ct2),]
head(melted_cts_cor1)
head(melted_cts_cor2)

df = data.frame(melted_cts_cor1, melted_cts_cor2)
myplots <- vector("list", 1)
#df = df[df$inf == "rest",]

#write.table(df, file =  "/Users/sanjavickovic/Desktop/morphoSPOT/Manuscript_Nov2021/figshare/SupFig9c_RA3.csv", sep = ",", quote = F) 

myplots[[1]] <- 
ggplot(data = df, aes(x = as.numeric(Avg), y = as.numeric(Avg.1))) + 
  facet_grid(~inf) + 
  geom_point(size = 2) + 
  geom_smooth(method="lm") +
  labs(x = paste0(ct1, " [%]"), y = paste0(ct2, " [%]")) + 
  labs(title = paste0("Cell type correlations")) + theme(legend.key.size = unit(0.2, "cm")) +
  theme(axis.text.x=element_text(angle=0, hjust=1)) +
  geom_smooth(method = "lm", formula = y~x,col="darkgray", se=FALSE) +
  stat_cor(method = "pearson", p.accuracy = 0.05) +
  facet_wrap(~ inf, scales = "free")
print(myplots[[1]])
save_plot(paste0(path_output,"Average_SC_sigs_per_inf_correlations_", ct1, "_", ct2, "_", norm_samples, ".pdf"), plot_grid(plotlist=myplots), ncol = 3, nrow = 2)

# plot correlation plots between cell type abundances in whole biosy
ct1 = paste0("scRNAseq_", tolower(norm_samples), "_", "fibro2b_thy1", "_clusters") #echange cell types here # fibro2b_thy1
#"fibroblasts_HLADRAsublining"
#"fibro1_cd55"
#"fibro2b_thy1"
unique(melted$ct)

ct2 = paste0("scRNAseq_", tolower(norm_samples), "_", "plasma_cells", "_clusters") #echange cell types here #b_cells
melted_cts_cor1 = melted[(melted$ct == ct1),] # & melted$inf != 'rest'
melted_cts_cor2 = melted[(melted$ct == ct2),]
df = data.frame(melted_cts_cor1, melted_cts_cor2)

#write.table(df, file = paste0("/Users/sanjavickovic/Desktop/morphoSPOT/Manuscript_Nov2021/figshare/SupFig9a_", norm_samples, "_plasma_thy1.csv"), sep = ",", quote = F) 

myplots <- vector("list", 1)
myplots[[1]] <- 
  ggplot(data = df, aes(x = as.numeric(Avg), y = as.numeric(Avg.1))) + 
  #facet_grid(~section) + 
  geom_point(size = 2) + 
  geom_smooth(method="lm") +
  labs(x = paste0(ct1, " [%]"), y = paste0(ct2, " [%]")) + 
  labs(title = paste0("Cell type correlations")) + theme(legend.key.size = unit(0.2, "cm")) +
  theme(axis.text.x=element_text(angle=0, hjust=1)) +
  geom_smooth(method = "lm", formula = y~x,col="darkgray", se=FALSE) +
  stat_cor(method = "pearson", p.accuracy = 0.05)
  #facet_wrap(~ inf, scales = "free")
print(myplots[[1]])
save_plot(paste0(path_output,"Average_SC_sigs_correlations_", ct1, "_", ct2, "_", norm_samples, ".pdf"), plot_grid(plotlist=myplots), ncol = 3, nrow = 2)

# save melted file
melted$sample = norm_samples
write.table(melted, paste0(path_output, "Average_sc_sigs_table_", norm_samples, ".tsv"), sep = "\t")

### This part of the code should be run when data on all samples is collected
# check if cell type abundance changes are significant in all samples
setwd("/Users/sanjavickovic/Desktop/morphoSPOT/3dst_repo/3dst/cell_typing")
fl = list.files(pattern = glob2rx("Average_sc_sigs_table_*"), path = path_output)
abun = matrix(nrow = 1, ncol = 9)
colnames(abun) = c("area","ct","Avg","inf","section","cluster" ,"ct_short","ra","sample")
for (f in fl){
  scabun = read.csv(paste0(path_output, f), sep = "\t")
  #print(f)
  #print(colnames(scabun))
  abun = rbind(abun, scabun)
  ncol(scabun)
}
abun = abun[-1,]

# do some renaming
nm = c(abun["ct"])
ctypes = ""
for (i in nm){
  ctypes = c(ctypes, paste(sapply(strsplit(i, "_"), "[[", 3),sapply(strsplit(i, "_"), "[[", 4), sep = "_"))
}
ctypes = ctypes[-1]
abun["ct"] = ctypes

# first checks general changes in abundances in all samples (ex. b cells)
by_vs_am <- abun %>% group_by(sample, ct, section)
by_vs <- by_vs_am %>% summarise(meanct = mean(Avg))
for (j in unique(by_vs$ct)){
  by_vs_ct = by_vs[by_vs$ct == j,]
  print(paste0("checking cell type: ", j))
  for (i in 1:6){
    print(paste0("RA",i))
    if (length(by_vs_ct[by_vs_ct$sample == paste0("RA",i),]$meanct)==0) next
    if (t.test(by_vs_ct[by_vs_ct$sample == paste0("RA",i),]$meanct,by_vs_ct[by_vs_ct$sample != paste0("RA",i),]$meanct ,alternative = "greater", p.adjust.method = "BH")$p.value < 0.05){
      print("significant")
    }
  }
}

# first checks general changes in abundances in infiltrates (ex. b cells)
by_vs_am <- abun %>% group_by(sample, ct, section, inf)
newinf = ""
for (i in by_vs_am$inf){
  newinf = c(newinf,substr(i, start = 1, stop = 3))
}
newinf = newinf[-1]
by_vs_am$inf = newinf
by_vs_am = by_vs_am[by_vs_am$inf == "inf",]
by_vs <- by_vs_am %>% summarise(meanct = mean(Avg))
for (j in unique(by_vs$ct)){
  by_vs_ct = by_vs[by_vs$ct == j,]
  print(paste0("checking cell type: ", j))
  for (i in 1:6){
    print(paste0("RA",i))
    if (length(by_vs_ct[by_vs_ct$sample == paste0("RA",i),]$meanct)==0) {
      print("no cells")
      next}
    if (t.test(by_vs_ct[by_vs_ct$sample == paste0("RA",i) & by_vs_ct$inf == "inf",]$meanct,by_vs_ct[by_vs_ct$sample != paste0("RA",i),]$meanct, alternative = "greater", p.adjust.method = "BH")$p.value < 0.05){
      print("significant")
      }
  }
}


# first checks general changes in abundances in Rest (ex. b cells)
by_vs_am <- abun %>% group_by(sample, ct, section, inf)
newinf = ""
for (i in by_vs_am$inf){
  newinf = c(newinf,substr(i, start = 1, stop = 3))
}
newinf = newinf[-1]
by_vs_am$inf = newinf
by_vs_am = by_vs_am[by_vs_am$inf == "res",]
by_vs <- by_vs_am %>% summarise(meanct = mean(Avg))
for (j in unique(by_vs$ct)){
  by_vs_ct = by_vs[by_vs$ct == j,]
  print(paste0("checking cell type: ", j))
  for (i in 1:5){
    print(paste0("RA",i))
    if (length(by_vs_ct[by_vs_ct$sample == paste0("RA",i),]$meanct)==0) {
      print("no cells")
      next}
    if (t.test(by_vs_ct[by_vs_ct$sample == paste0("RA",i) & by_vs_ct$inf == "res",]$meanct,by_vs_ct[by_vs_ct$sample != paste0("RA",i),]$meanct, alternative = "greater")$p.value < 0.05){
      print("significant")
    }
  }
}

# plot all density sc plots after all samples have been analyzed 
by_vs_am <- abun %>% group_by(sample, ct, section, inf)

by_vs_am = abun%>%
  mutate(serostatus = case_when(
    endsWith(sample, "1") ~ "seropositive",
    endsWith(sample, "2") ~ "seropositive",
    endsWith(sample, "3") ~ "seropositive",
    endsWith(sample, "4") ~ "seronegative",
    endsWith(sample, "5") ~ "seronegative",
    endsWith(sample, "6") ~ "seronegative"
  ))


#write.table(by_vs_am, file =  "/Users/sanjavickovic/Desktop/morphoSPOT/Manuscript_Nov2021/figshare/Fig4c.csv", sep = ",", quote = F)

counter = 1
myplots <- vector("list", length(unique(by_vs_am$ct)))
for (j in unique(by_vs_am$ct_short)){
  print(j)
  meltedsub = by_vs_am[by_vs_am$ct_short == j,]
  meltedsub = meltedsub[meltedsub$cluster == "inf",]
  myplots[[counter]] <- ggplot(meltedsub, aes(Avg, y = factor(ct, levels = paste0("celltype_", 1:length(unique(meltedsub$ct)))), fill = serostatus)) +
    ggridges::geom_density_ridges(alpha = 0.5) +
    scale_x_continuous(limits = c(0, 15)) +
    labs(y = j) +
    scale_fill_manual(values = c("seropositive" = "#CC6677","seronegative" = "#4477AA")) +
    theme_minimal()
  counter = counter + 1
}
save_plot(paste0(path_output,"Average_SC_inf_samples_density.pdf"), plot_grid(plotlist=myplots), ncol = 4, nrow = 4)

counter = 1
for (j in unique(by_vs_am$ct_short)){
  print(j)
  meltedsub = by_vs_am[by_vs_am$ct_short == j,]
  meltedsub = meltedsub[meltedsub$cluster == "res",]
  myplots[[counter]] <- ggplot(meltedsub, aes(Avg, y = factor(ct, levels = paste0("celltype_", 1:length(unique(meltedsub$ct)))), fill = serostatus)) +
    ggridges::geom_density_ridges(alpha = 0.5) +
    scale_x_continuous(limits = c(0, 12.5)) +
    labs(y = j) +
    scale_fill_manual(values = c("seronegative" = "#4477AA", "seropositive" = "#CC6677")) +
    theme_minimal()
  counter = counter + 1
}
save_plot(paste0(path_output,"Average_SC_rest_samples_density.pdf"), plot_grid(plotlist=myplots), ncol = 4, nrow = 4)





