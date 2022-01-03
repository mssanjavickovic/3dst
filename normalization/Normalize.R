# Example normalization for all RA samples
#Load libraries
suppressMessages(suppressWarnings(library(scran,warn.conflicts = F, quietly = T)))
suppressMessages(suppressWarnings(library(SingleCellExperiment,warn.conflicts = F, quietly = T)))
suppressMessages(suppressWarnings(library(stringr,warn.conflicts = F, quietly = T)))
suppressMessages(suppressWarnings(library(scater,warn.conflicts = F, quietly = T)))

# removes all glob variables from env
rm(list = ls())

# Read the fuctions file 
setwd("/Users/sanjavickovic/Desktop/morphoSPOT/3dst_repo/3dst/normalization")
source('../functions/Read_Functions.R')

# Which samples do you want to normalize? 
norm_samples = "RA6"

# Where are you raw csv expression files located? 
path_samples = paste0("/Users/sanjavickovic/Desktop/morphoSPOT/Raw counts and images for paper/", norm_samples)

# Set output directory that will contain all cell_typing pdf plots and output gene files
path_output = "/Users/sanjavickovic/Desktop/morphoSPOT/norm_counts_output"

# List avaiable raw matrices for normalization 
files_raw = list.files(pattern =  glob2rx(paste0(norm_samples, "_Raw.exp_*.csv")), path = path_samples)

# Load all raw counts matrices in your env
sample_counter = 1
for (i in files_raw){
  print(i)
  assign(paste0("raw_mat", sample_counter), read.csv(paste0(path_samples, "/", i), header = T, row.names = 1))
  sample_counter = sample_counter + 1
}
raw_mat_loaded_in_env <- grep("raw_mat",names(.GlobalEnv),value=TRUE)

# Combine all the raw data in one matrix
exp.values = do.call(mbind,mget(raw_mat_loaded_in_env))
exp.values[is.na(exp.values)] <- 0

# Label sections as covariates
section.labels = ""
section.data = ""
for (secs in 1:length(raw_mat_loaded_in_env)){
  section.labels = c(section.labels, paste0(paste0("X",secs, "_"), str_replace(colnames(get(paste0("raw_mat", secs))), "X", "")))
  section.data = c(section.data, rep(paste0("X",secs), ncol(get(paste0("raw_mat", secs)))))
}
section.labels = section.labels[-1]
section.data = section.data[-1]

# Rename colnames to reflect sections from a single patient biopsy
colnames(exp.values) = section.labels 

# Make a SCE object required by scran
sce <- SingleCellExperiment(list(counts=as.matrix(exp.values)),
                            colData=DataFrame(section=section.data),
                            rowData=DataFrame(label=row.names(exp.values)))

# Compute some lib metrics
# for more details on scran normalization, please check Aaron's tutorial: https://bioconductor.org/packages/release/bioc/vignettes/scran/inst/doc/scran.html#33_computing_normalized_expression_values
counts <- assay(sce)
libsizes <- colSums(counts)
size.factors <- libsizes/mean(libsizes)
qcstats <- perCellQCMetrics(sce)
qcfilter <- quickPerCellQC(qcstats)
sce <- sce[,!qcfilter$discard]
summary(qcfilter$discard) # true denotes number of discarded ST spots

# Performes normalization
clusters <- quickCluster(sce,  min.size=200) # do not want very small clusters, #20 for most, 200 for RA6
sce <- computeSumFactors(sce, clusters=clusters, positive = TRUE)
summary(sizeFactors(sce))
sce <- logNormCounts(sce) # this represents the log normalized object

# Output normalized counts into separate files and R objects
norm.counts = assay(sce, 'logcounts')

# Change dir to some output dir 
setwd(path_output)

###output norm counts
for (secs in 1:length(raw_mat_loaded_in_env)){
  assign(paste0("m", secs), assay(sce, "logcounts")[,grep(pattern = paste0("X", secs), colnames(norm.counts), value = F)])
  print(paste0("Saving normalized file ... : ", paste0(norm_samples, '_Norm.exp_', secs, '.csv')))
  save(list=paste0("m", secs), file = paste0(norm_samples, '_norm.exp.values.', secs))
  write.table(get(paste0("m", secs)), file = paste0(norm_samples, '_Norm.exp_', secs, '.csv'), sep = ",", quote = F, col.names = T, row.names = T)
}

