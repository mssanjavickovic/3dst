### Go terms analysis
suppressPackageStartupMessages({
  library(openxlsx)
  library(clusterProfiler)
})

# removes all glob variables from env
rm(list = ls())

setwd('/Users/sanjavickovic/Desktop/morphoSPOT/3dst_repo/3dst/de_analysis')
# setwd("/Users/svickovi/Library/Mobile Documents/com~apple~CloudDocs/Desktop/morphoSPOT/3dst_repo/3dst/go_analysis")

# Loads files for ST analysis
path_output = "../../../plot_outputs"

# Which samples do you want to use? 
norm_samples = "RA5"

txt.files <- list.files(pattern = glob2rx(paste0("DEGs.Inf.Cluster*_subsetted*",norm_samples, "*.txt")), path = path_output, full.names = T)
#txt.files <- list.files(pattern = glob2rx(paste0("DEGs.Spatial.Cluster*_subsetted*",norm_samples, "*.txt")), path = path_output, full.names = T)
# txt.files <- list.files(pattern = glob2rx(paste0("*",tolower(norm_samples), "*_clusters_subsetted.txt")), path = path_output, full.names = T)


celltypes <-  sapply(strsplit(sapply(strsplit(txt.files, split = "/"),"[[",5), "_"),"[[",1)
# celltypes <-  str_replace(sapply(strsplit(str_replace(txt.files, "_clusters_subsetted.txt", ""), "/"),"[[",5),"scRNAseq_","")


GeneSets <- setNames(lapply(txt.files, function(cell.type.file) {
  read.table(cell.type.file, header = T, stringsAsFactors = F)$x
}), nm = celltypes)


h <- suppressPackageStartupMessages({read.gmt("../data/c5.bp.v6.2.symbols.gmt")})

bp.eriched.terms <- setNames(lapply(celltypes, function(cl) {
  egmt <- enricher(GeneSets[[cl]], TERM2GENE = h, pvalueCutoff = 0.05)
  as.data.frame(egmt)
}), nm = celltypes)

sh <- createWorkbook(title = paste0(norm_samples, " inf clusters GO enrcihment"))
# sh <- createWorkbook(title = paste0(norm_samples, " celltypes GO enrcihment"))
# sh <- createWorkbook(title = paste0(norm_samples, " spatial clusters GO enrcihment"))

for (cl in celltypes) {
  addWorksheet(wb = sh, sheetName = cl)
  firstrow <- paste0(norm_samples," cluster: ", gsub(pattern = "Cluster", replacement = "", x = cl))
  writeData(wb = sh, x = firstrow, sheet = cl)
  secondrow <- paste0("marker genes: ", paste(GeneSets[[cl]], collapse = ", "))
  writeData(wb = sh, x = secondrow, sheet = cl, startRow = 2)
  writeData(wb = sh, x = "", sheet = cl, startRow = 3)
  writeData(wb = sh, x = bp.eriched.terms[[cl]][, -1], sheet = cl, startRow = 4)
  setColWidths(wb = sh, sheet = cl, cols = 1:ncol(bp.eriched.terms[[cl]]), widths = "auto")
  print(bp.eriched.terms[[cl]][, -1])
}

saveWorkbook(sh, file = paste0(path_output, "/", norm_samples, "_inf_clusters_GO_enrichment.xlsx"), overwrite = TRUE)
# saveWorkbook(sh, file = paste0(path_output, "/", norm_samples, "_cell_types_GO_enrichment.xlsx"), overwrite = TRUE)
#saveWorkbook(sh, file = paste0(path_output, "/", norm_samples, "_spatial_clusters_GO_enrichment.xlsx"), overwrite = TRUE)          




# make one big doc with all top 3 pathways 
# Which samples do you want to use? 
norm_samples = c("RA1", "RA2","RA3","RA4","RA5")
go_infs = matrix(ncol = 3, nrow = 1)
for (samples in norm_samples){
  print(samples)
  txt.files <- list.files(pattern = glob2rx(paste0("DEGs.Inf.Cluster*_subsetted*",samples, "*.txt")), path = path_output, full.names = T)
  # txt.files <- list.files(pattern = glob2rx(paste0("DEGs.Spatial.Cluster*_subsetted*",norm_samples, "*.txt")), path = path_output, full.names = T)
  #txt.files <- list.files(pattern = glob2rx(paste0("*",tolower(norm_samples), "*_clusters_subsetted.txt")), path = path_output, full.names = T)
  
  celltypes <-  sapply(strsplit(sapply(strsplit(txt.files, split = "/"),"[[",5), "_"),"[[",1)
  # celltypes <-  str_replace(sapply(strsplit(str_replace(txt.files, "_clusters_subsetted.txt", ""), "/"),"[[",5),"scRNAseq_","")

  GeneSets <- setNames(lapply(txt.files, function(cell.type.file) {
    read.table(cell.type.file, header = T, stringsAsFactors = F)$x
  }), nm = celltypes)

  h <- suppressPackageStartupMessages({read.gmt("../data/c5.bp.v6.2.symbols.gmt")})

  bp.eriched.terms <- setNames(lapply(celltypes, function(cl) {
    egmt <- enricher(GeneSets[[cl]], TERM2GENE = h, pvalueCutoff = 0.05)
    as.data.frame(egmt)
  }), nm = celltypes)

  for (cl in celltypes) {
    print(cl)
    if (dim(bp.eriched.terms[[cl]][, -1])[1] == 0) next
    tmp =  bp.eriched.terms[[cl]][, -1]
    #tmp = tmp[tmp$GeneRatio %in% sort(tmp$GeneRatio, decreasing = T)[1:20],]
    go_infs = rbind(go_infs, cbind(paste0(samples, "_", cl), row.names(tmp), tmp$geneID))
  }
}
go_infs = go_infs[-1,]
row.names(go_infs) = go_infs[,1]
go_infs = go_infs[,-1]
# go_infs = go_infs[!duplicated(row.names(go_infs)),]
colnames(go_infs) = c("pathway", "genes")
go_infs = data.frame(go_infs)
go_infs$sample = sapply(strsplit(row.names(go_infs), "_"),"[[",1)
go_infs$cluster = sapply(strsplit(row.names(go_infs), "\\."),"[[",3)
#go_infs$cluster = sapply(strsplit(sapply(strsplit(row.names(go_infs), "\\."),"[[",1),"_"),"[[",3)

gn = ""
for (p in names(sort(table(go_infs[,1]), decreasing = T)[1:5])){
  sub = go_infs[go_infs$pathway == p,]
  print(p)
  for (nm in unique(sub$cluster)){
    print(nm)
    gn = c(unlist(c(strsplit(sub$genes, "/"))))
    #gn = c(gn, sub$genes)
    print(sort(table(gn), decreasing = T)[1:5])
  }
}



            



