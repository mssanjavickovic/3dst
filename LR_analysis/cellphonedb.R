library(dplyr)
library(tidyr)
library(magrittr)
library(ggplot2)
library(stringr)

setwd("/Users/sanjavickovic/Desktop/morphoSPOT/3dst_repo/3dst/go_analysis")

merge.means.pvals <- function(f1, f2) {
  means <- read.table(file = f1, sep = "\t", header = TRUE)
  pvals <- read.table(file = f2, sep = "\t", header = TRUE)
  stopifnot("Some columns are missing. Check input paths..." = length(intersect(colnames(means), colnames(pvals))) == 15)
  pvals_significant <- apply(pvals[, 12:15], 1, function(x) {any(x < 0.01)})
  pvals <- pvals[pvals_significant, ]
  print(head(means))
  means_significant <- apply(means[, 12:15], 1, function(x) {any(log10(x) > 0.1)})
  means <- means[means_significant, ]
  print(head(means))
  ggm.means <- gather(data = means, key = "comparison", value = "means", 12:15)
  ggm.pvals <- gather(data = pvals, key = "comparison", value = "pvalue", 12:15)
  ggm <- merge(ggm.means, ggm.pvals)
  ggm$pvalue <- ifelse(ggm$pvalue == 0, 0.0009, ggm$pvalue)
  ggm$comparison <- gsub(pattern = "\\.", replacement = " | ", x = ggm$comparison)
  return(ggm)
}
DotPlot <- function(ggm) {
  jet.colors <-c("#00007F", "blue", "#007FFF", "cyan",
                       "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000")
  print(jet.colors)
  means_spread <- spread(data = ggm[, setdiff(colnames(ggm), c("id_cp_interaction", "partner_a", "partner_b", "gene_a", 
                                                               "gene_b", "secreted", "receptor_a", "receptor_b", 
                                                               "is_integrin", "annotation_strategy", "pvalue"))], 
                         key = "comparison", value = "means") %>%
    data.frame(row.names = 1, check.names = FALSE)
  dMat <- dist(log10(means_spread))
  tree <- hclust(dMat)
  ggm$interacting_pair <- factor(ggm$interacting_pair, levels = rownames(means_spread)[tree$order])
  a = str_replace(sapply(str_split(ggm$comparison, "\\|"), "[[", 1), " ", "")
  b = str_replace(sapply(str_split(ggm$comparison, "\\|"), "[[", 2), " ", "")
  ggm = ggm[a == b,]
  ggplot(ggm, aes(comparison, interacting_pair, size = -log10(pvalue), fill = log10(means))) +
    geom_point(shape = 21) +
    scale_fill_gradientn(colours = jet.colors) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
}


RA1to3_cluster_1_vs_all <- merge.means.pvals(f1 = "../data/RA1to3_cluster_1_vs_all/means.txt",
                                             f2 = "../data/RA1to3_cluster_1_vs_all/pvalues.txt")
p <- DotPlot(RA1to3_cluster_1_vs_all)
p
RA4to6_cluster_1_vs_all <- merge.means.pvals(f1 = "../data/RA4to6_cluster_1_vs_all/means.txt",
                                             f2 = "../data/RA4to6_cluster_1_vs_all/pvalues.txt")
p <- DotPlot(RA4to6_cluster_1_vs_all)
p
seropositive_vs_seronegative <- merge.means.pvals(f1 = "../data/seropositive_vs_seronegative_cluster_1/means.txt",
                                                  f2 = "../data/seropositive_vs_seronegative_cluster_1/pvalues.txt")

write.table(seropositive_vs_seronegative, file = "/Users/sanjavickovic/Desktop/morphoSPOT/Manuscript_Nov2021/figshare/Fig4b.csv", sep = ",", quote = F)
p <- DotPlot(seropositive_vs_seronegative)
p

print(RColorBrewer::brewer.pal(n = 9, name = "Spectral"))
