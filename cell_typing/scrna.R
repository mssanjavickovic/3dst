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

# setwd("/Users/svickovi/Library/Mobile Documents/com~apple~CloudDocs/Desktop/morphoSPOT/3dst_repo/3dst/cell_typing")
setwd("/Users/sanjavickovic/Desktop/morphoSPOT/3dst_repo/3dst/cell_typing")

# removes all glob variables from env
rm(list = ls())

# Run sc imputation for ex. RA1 samples
# Which samples do you want to use? 
norm_samples = "RA1"

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
co = 0.2 # Set cut off for gene-to-gene correlations, recommended 0.2

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
  write.table(tmp1, file = paste0("scRNAseq",i,".", norm_samples, ".csv"), sep =",", quote =F, row.names = F, col.names = T)
}

# Plot spatial heatmaps per cell type
# setwd(path_output)
for (i in row.names(mf1)){
  # Plot spatial spots
  ## (un)comment for specific samples ex. use only plot.gene.2d.4 for RA1
  if (norm_samples == 'RA1') plot.gene.2d.4("RA1_scRNAseq",i, mn1, mn2, mn3, mn4, s1, s2, s3, s4, con = T, x=40, y=20, transparency=1, min=min(mf1), max=max(mf1)) # for RA1
  if (norm_samples == 'RA2') plot.gene.2d.7("RA2_scRNAseq",i, mn1, mn2, mn3, mn4, mn5, mn6, mn7, s1, s2, s3, s4, s5, s6, s7, con = T, x=40, y=20, transparency=1, min=min(mf1), max=max(mf1)) # for RA2
  #if (norm_samples == 'RA2') plot.gene.3d.7("RA2_scRNAseq",i, mn1, mn2, mn3, mn4, mn5, mn6, mn7, s1, s2, s3, s4, s5, s6, s7, x=40, y=20, transparency=0.1, min=min(mf1), max=max(mf1)) # for RA2
  if (norm_samples == 'RA3') plot.gene.2d.4("RA3_scRNAseq",i, mn1, mn2, mn3, mn4, s1, s2, s3, s4, con = T, x=40, y=40, transparency=1, min=min(mf1), max=max(mf1)) # for RA3
  if (norm_samples == 'RA4') plot.gene.2d.5("RA4_scRNAseq",i, mn1, mn2, mn3, mn4, mn5, s1, s2, s3, s4, s5, con = T, x=40, y=40, transparency=1, min=min(mf1), max=max(mf1)) # for RA4
  # if (norm_samples == 'RA5') plot.gene.3d.4("RA5_scRNAseq",i, mn1, mn2, mn3, mn4, s1, s2, s3, s4, x=80, y=60, transparency=0.2, min=min(mf1), max=max(mf1)) # for RA5
  if (norm_samples == 'RA5') plot.gene.2d.3("RA5_scRNAseq",i, mn1, mn2, mn3, s1, s2, s3, con = T, x=80, y=80, transparency=1, min=min(mf1), max=max(mf1)) # for RA5
}

# Avg cell signature per infiltrate 
### Take only infiltrates that are present in all sections according to annotations
# setwd("/Users/sanjavickovic/Desktop/morphoSPOT/3dst_repo/3dst/cell_typing")
setwd("/Users/svickovi/Library/Mobile Documents/com~apple~CloudDocs/Desktop/morphoSPOT/3dst_repo/3dst/cell_typing")
inf_all = read.csv(paste0("../data/", norm_samples, "_zstack_Infs.csv"), header = T)
mf1 = fil_m

# Subset per infiltrate as the naming ie. Inf1 presumes we are following the same infiltrate (given 3D alignment) in all present sections
all_inf = ""
for (i in 1:length(unique(inf_all$Inf.))){
  assign(paste0("SC.norm.inf", i), as.matrix(mf1[,colnames(mf1) %in% row.names(subset(inf_all, inf_all$Inf. == paste0("Inf", i)))]))
  all_inf = c(all_inf, colnames(get(paste0("SC.norm.inf", i)))) # Make an object with all annotated infiltrates
}
all_inf = all_inf[-1]

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
library(reshape2) # for melt
melted <- melt(rbind(Infs, t(avg_rest)))
colnames(melted) = c("area","ct","Avg")
melted$inf <- sapply(strsplit(as.character(melted$area), "_"),"[[",1)
melted$section <- sapply(strsplit(as.character(melted$area), "_"),"[[",2)
set3 <- colorRampPalette(c(brewer.pal('Set3',n=12),c(brewer.pal('Set1',n=9))))
setnm = setNames(set3(13), sapply(strsplit(cluster_files, ".txt",""), "[[",1))
setnm_final = setnm[names(setnm) %in% unique(melted$ct)]

myplots[[1]] <-ggplot(data = melted, aes(x = inf, y = Avg, fill = ct)) + 
  geom_bar(position="fill", stat="identity") + 
  facet_grid(~section) + 
  labs(x = "", y = "Cell type [%]") + 
  scale_fill_manual(name = "Cell types", values = setnm_final) +
  labs(title = paste0("Cell type scores")) + theme(legend.key.size = unit(0.2, "cm"))
save_plot(paste0(path_output,"Average_SC_sigs_per_section_", norm_samples, ".pdf"), plot_grid(plotlist=myplots), ncol = 3, nrow = 2)

# plot Inf and rest together per infiltrate
myplots <- vector("list", 1)
myplots[[1]] <- ggplot(data = melted, aes(x = section, y = Avg, fill = ct)) + 
  geom_bar(position="fill", stat="identity") + 
  facet_grid(~inf) + 
  labs(x = "", y = "Cell type [%]") + 
  scale_fill_manual(name = "Cell types", values = setnm_final) +
  labs(title = paste0("Cell type scores")) + theme(legend.key.size = unit(0.2, "cm")) +
  theme(axis.text.x=element_text(angle=30, hjust=1))
save_plot(paste0(path_output,"Average_SC_sigs_per_inf_", norm_samples, ".pdf"), plot_grid(plotlist=myplots), ncol = 3, nrow = 2)

# plot correlation plots between cell type abundances
ct1 = paste0("scRNAseq_", tolower(norm_samples), "_", "dendritic_cells", "_clusters") #echange cell types here
ct2 = paste0("scRNAseq_", tolower(norm_samples), "_", "b_cells", "_clusters") #echange cell types here #fibro2b_thy1
melted$ct
melted_cts_cor1 = melted[(melted$ct == ct1),] # & melted$inf != 'rest'
melted_cts_cor2 = melted[(melted$ct == ct2),]
df = data.frame(melted_cts_cor1, melted_cts_cor2)
myplots <- vector("list", 1)
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
save_plot(paste0(path_output,"Average_SC_sigs_per_inf_correlations_", ct1, "_", ct2, "_", norm_samples, ".pdf"), plot_grid(plotlist=myplots), ncol = 3, nrow = 2)
     

