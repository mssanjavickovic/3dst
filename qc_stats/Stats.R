### Script to output seq quality stats for RA samples 
suppressMessages(suppressWarnings(library(stringr,warn.conflicts = F, quietly = T)))
suppressMessages(suppressWarnings(library(ggplot2,warn.conflicts = F, quietly = T)))
suppressMessages(suppressWarnings(library(cowplot,warn.conflicts = F, quietly = T)))
# removes all glob variables from env
rm(list = ls())

# Read the fuctions file 
setwd("/Users/svickovi/Library/Mobile Documents/com~apple~CloudDocs/Desktop/morphoSPOT/3dst_repo/3dst/data")
source('../functions/Read_Functions.R')

# Where are you raw csv expression files located? 
path_samples = list.files("../data/", pattern = glob2rx("RA*Raw.exp*"))

# Set output directory that will contain all plots and output gene files
path_output = "/Users/svickovi/Library/Mobile Documents/com~apple~CloudDocs/Desktop/morphoSPOT/plot_output"


stats = matrix(ncol = 6, nrow = 1)
for (ra in unique(sapply(str_split(path_samples, "_"),"[[",1))){

  # List avaiable raw matrices for normalization 
  files_raw = path_samples[grep(path_samples, pattern = ra)]

  # Load all raw counts matrices in your env
  sample_counter = 1
  trans_means = ""
  genes_means = ""
  ra_ncols = ""
  for (i in files_raw){
    print(i)
    assign(paste0("raw_mat", sample_counter), read.csv(i, header = T, row.names = 1))
    trans_means = c(trans_means, mean(colSums(get(paste0("raw_mat", sample_counter)))))
    genes_means = c(genes_means, mean(colSums(get(paste0("raw_mat", sample_counter))>0)))
    ra_ncols = c(ra_ncols, ncol(get(paste0("raw_mat", sample_counter))))
    sample_counter = sample_counter + 1
  }
  raw_mat_loaded_in_env <- grep("raw_mat",names(.GlobalEnv),value=TRUE)
  trans_means = as.numeric(trans_means[-1])
  genes_means = as.numeric(genes_means[-1])
  ra_ncols = as.numeric(ra_ncols[-1])
  ra_size = as.numeric(length(unique(raw_mat_loaded_in_env)))
  trans_sd = as.numeric(sd(trans_means)/sqrt(ra_size))
  genes_sd = as.numeric(sd(genes_means)/sqrt(ra_size))
  trans_means = mean(trans_means)
  genes_means = mean(genes_means)
  
  if (ra == "RA1"){
    seq = (c(60000000,60000000,60000000,60000000)/ra_ncols)/100
    seq_sd = sd(seq)/sqrt(length(seq))
    seq_mean = mean(seq)
  }
  if (ra == "RA2"){
    seq = (c(69746243,79676431,52784882,52440887,60719029,57316904,71892234)/ra_ncols)/100
    seq_sd = sd(seq)/sqrt(length(seq))
    seq_mean = mean(seq)
  }
  if (ra == "RA3"){
    seq = (c(65765190,51577431,61293801,58890473)/ra_ncols)/100
    seq_sd = sd(seq)/sqrt(length(seq))
    seq_mean = mean(seq)
  }
  if (ra == "RA4"){
    seq = (c(59322039,64288682,47225356,53125282, 63085300)/ra_ncols)/100
    seq_sd = sd(seq)/sqrt(length(seq))
    seq_mean = mean(seq)
  }
  if (ra == "RA5"){
    seq = (c(14999492,45159498,62707152,61266189)/ra_ncols)/100
    seq_sd = sd(seq)/sqrt(length(seq))
    seq_mean = mean(seq)
  }
  
  stats = rbind(stats, cbind(seq_mean, seq_sd, trans_means, trans_sd, genes_means, genes_sd))
  
}
stats = stats[-1,]
row.names(stats) = unique(sapply(str_split(path_samples, "_"),"[[",1))

# Plot avg gene expression for each Inf separately
counter = 1
myplots <- vector("list", 3)
for(i in 1:3){
  print(counter)
  print(counter+1)
  myData = data.frame(stats[,counter:(counter+1)])
  colnames(myData) = c("Avg", "SEM")
  limits <- aes(ymax = Avg + SEM, ymin = Avg - SEM)
  if (i == 1) {y_lab = "Avg seq depth count per ST spot x 10E2"}
  if (i == 2) {y_lab = "Avg UMI count per ST spot"}
  if (i == 3) {y_lab = "Avg gene count per ST spot"}
  myplots[[i]] <- ggplot(data = myData, aes(x = factor(row.names(myData)), y = Avg, fill = row.names(myData))) + 
    geom_bar(stat = "identity", position = position_dodge(0.99)) + 
    geom_errorbar(limits, position = position_dodge(0.99), width = 0.25) + 
    scale_fill_manual(name = "RA patient biopsy", values = c("firebrick4","firebrick","firebrick3","firebrick2","firebrick1")) + 
    labs(x = "", y = "Number of observations") + 
    labs(title = y_lab) + 
    ylab(y_lab)
  counter = counter + 2
}
plot_grid(plotlist=myplots)