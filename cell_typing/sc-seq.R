#### Single cell associations in ST RA ###
# Load libraries
library(plot3D)
library(gdata)
library(qdapTools)
library(data.table)
library(akima)

# Load data as R objects for RA1
load('RA1_norm.exp.values.1') # in ./data
load('RA1_norm.exp.values.2') # in ./data
load('RA1_norm.exp.values.3') # in ./data
load('RA1_norm.exp.values.4') # in ./data
load('RA1_1_selected_adjusted_spots_3D_manual_app') # in ./data
load('RA1_2_selected_adjusted_spots_3D_manual_app') # in ./data
load('RA1_3_selected_adjusted_spots_3D_manual_app') # in ./data
load('RA1_4_selected_adjusted_spots_3D_manual_app') # in ./data

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

# Source f(x) file
source("Read_Functions.R") # in ./functions

# Combine in one file
RA.norm = cbind(m1, m2, m3, m4) # RA1
RA.norm = cbind(m1, m2, m3, m4, m5, m6, m7) # RA2
m = RA.norm # this is a combined matrix for all sections from RA1 patient with 4 norm expression values (individual files avaiable in ./data)

# Read in sc data
degs = read.csv("Stephenson.clusters.csv", header = T) # This is a table with DE genes avaiable from Stephenson et al
dates = read.delim("Genes.Dates.txt", header = F) # This is a table provided by 10X with all "fail gene name conversions from Excel"

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
m.use = markers[avg_logFC > 1][p_val_adj < .05][order(-avg_logFC)][,gene[1:200],cluster]

# Convert to named list
m.use = tapply(m.use$V1, m.use$cluster, function(a) sort(unique(a)))

# Then score these in the expression matrix you have, e.g.
scores = sapply(m.use, function(genes.use) colSums(m[row.names(m) %in% genes.use,]))

# Name and output if needed
names(m.use) = paste("scRNAseq_ra1", names(m.use), sep="_")
names(m.use) = paste("scRNAseq_ra2", names(m.use), sep="_")
for (i in 1:length(names(m.use))){
  nm = m.use[i]
  write.table(list2df(nm)$X1, file = paste(names(m.use)[i], "_clusters.txt", sep=""), sep="\t", quote =F)
}

# Choose which cluster files you want to use based on names
cluster_files = list.files('./', pattern=glob2rx("scRNAseq_ra1_*_clusters.txt"))
cluster_files = list.files('./', pattern=glob2rx("scRNAseq_ra2_*_clusters.txt"))

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
    rn_means <- mean(rn)
    rn_cent <- (rn-abs(rn_means))/abs(sd(rn))
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
      if (nrow(z) < 2) next # if less than 3 genes are coexpressed, skip
    
      # Subest based on  cut off and save genes per cluster
      n <- (subset(z, z [,3] > co)) [,1:2]
      rn_ref <- c(rn_ref, union(n[,1], n[,2]))
      un = union(n[,1], n[,2])
      if (length(un) == 0) next # if no genes passed filter, skip
      save(un, file = paste("cluster_genes_ra1", i, sep=""))
      genes_selected = c(genes_selected, un)
      write.table(un, file = paste0(sample, "_subsetted.txt"), quote = F, row.names = T) # this is the table as inpit to GO term analysis
      
      # Subset based on "good" genes per cluster and center the data
      rn <- m_sub[row.names(m_sub) %in% un,]
      rownames(rn) <- c(paste(sample,rownames(rn), sep="_"))
      rn_means <- colMeans(rn)
      rn_cent <- (rn_means-abs(mean(rn_means)))/abs(sd(rn_means))
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
mf1 = fil_m # precaution

# Split mf1 matrix into smaller matrices per section for RA1
mn1 = mf1[,sapply(strsplit(colnames(mf1), split="_"),"[[",1) == "X1"] # use this sections order for RA1
mn2 = mf1[,sapply(strsplit(colnames(mf1), split="_"),"[[",1) == "X2"] # use this sections order for RA1 
mn3 = mf1[,sapply(strsplit(colnames(mf1), split="_"),"[[",1) == "X3"] # use this sections order for RA1
mn4 = mf1[,sapply(strsplit(colnames(mf1), split="_"),"[[",1) == "X4"] # use this sections order for RA1

# Split mf1 matrix into smaller matrices per section for RA2
mn1 = mf1[,sapply(strsplit(colnames(mf1), split="_"),"[[",1) == "X1"] # use this sections order for RA2
mn2 = mf1[,sapply(strsplit(colnames(mf1), split="_"),"[[",1) == "X2"] # use this sections order for RA2 
mn3 = mf1[,sapply(strsplit(colnames(mf1), split="_"),"[[",1) == "X3"] # use this sections order for RA2
mn4 = mf1[,sapply(strsplit(colnames(mf1), split="_"),"[[",1) == "X4"] # use this sections order for RA2
mn5 = mf1[,sapply(strsplit(colnames(mf1), split="_"),"[[",1) == "X5"] # use this sections order for RA2 
mn6 = mf1[,sapply(strsplit(colnames(mf1), split="_"),"[[",1) == "X6"] # use this sections order for RA2
mn7 = mf1[,sapply(strsplit(colnames(mf1), split="_"),"[[",1) == "X7"] # use this sections order for RA2

# Plot spatial heatmaps per cell type
for (i in row.names(mf1)){
  # Plot spatial spots
  plot.gene.2d.4("RA1_scRNAseq_",i, mn1, mn2, mn3, mn4, s1, s2, s3, s4, con = T, x=40, y=20, transparency=1, min=min(mf1), max=max(mf1)) # for RA1
  # plot.gene.2d.9("RA2_scRNAseq_",i, mn1, mn2, mn3, mn4, mn5, mn6, mn7, s1, s2, s3, s4, s5, s6, s7, con = T, x=40, y=20, transparency=1, min=min(mf1), max=max(mf1)) # for RA2
}


