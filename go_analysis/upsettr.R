### Finding common genes between spatial clusters and patient biopsies 
# removes all glob variables from env
rm(list = ls())

# Loading libraries
suppressMessages(suppressWarnings(library(ggplot2,warn.conflicts = F, quietly = T)))
suppressMessages(suppressWarnings(library(stringr,warn.conflicts = F, quietly = T)))
suppressMessages(suppressWarnings(library(UpSetR,warn.conflicts = F, quietly = T)))

# Read the fuctions file 
setwd("/Users/svickovi/Library/Mobile Documents/com~apple~CloudDocs/Desktop/morphoSPOT/3dst_repo/3dst/data")
source('../functions/Read_Functions.R')

# First find common genes for all sections and cluster#1 etc 
# Set output directory that will contain all plots and output gene files
path_output = "/Users/svickovi/Library/Mobile Documents/com~apple~CloudDocs/Desktop/morphoSPOT/plot_outputs/"
setwd(path_output)

# Where are you raw csv expression files located? 
path_samples = list.files("./", pattern = glob2rx("DEGs.Cluster*"))

for (i in path_samples){
  print(i)
  nm = str_replace(str_replace(i, "DEGs.", ""), ".csv", "")
  assign(nm, cbind(as.matrix(row.names(read.csv(i, header = T, row.names = 1)), 1), as.numeric(1)))
  nm_tmp  = get(nm)
  row.names(nm_tmp) = nm_tmp[,1]
  nm_tmp = as.matrix(nm_tmp[,-1])
  assign(nm, nm_tmp)
}
clusters_loaded_in_env <- grep("Cluster", names(.GlobalEnv),value=TRUE)

for (j in clusters_loaded_in_env){
  cluster_counter = str_replace(sapply(str_split(j, ".RA"),"[[",1), "Cluster","")
  print(cluster_counter)
  assign(paste0("Cluster",cluster_counter,".norm"), do.call(mbind, mget(clusters_loaded_in_env[grep(paste0("Cluster", cluster_counter),clusters_loaded_in_env, value = F)])))
  tmp = get(paste0("Cluster",cluster_counter,".norm"))
  tmp[is.na(tmp)] <- 0
  colnames(tmp) = clusters_loaded_in_env[grep(paste0("Cluster", cluster_counter),clusters_loaded_in_env, value = F)]
  assign(paste0("Cluster",cluster_counter,".norm"),tmp)
}
numbered_clusters_loaded_in_env <- grep(glob2rx("Cluster*norm"), names(.GlobalEnv),value=TRUE)

# Plot DE gene names intersects between different biopsy clusters
for (x in numbered_clusters_loaded_in_env){
  gnmat = get(x)
  nms = row.names(gnmat)
  gnmat <- data.frame(apply(gnmat, 2, function(x) as.numeric(as.character(x))))
  row.names(gnmat) = nms
  gnmat = as.data.frame(t(gnmat))
  gnmat$Identifier = row.names(gnmat)
 print(upset(gnmat, sets = c(colnames(gnmat)[1:(length(colnames(gnmat))-1)]), sets.bar.color = "#56B4E9",
        order.by = "freq", nintersects = NA))
}


Cluster4.RA2
