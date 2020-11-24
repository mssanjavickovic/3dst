### Script for grouping infiltrates into separate spatial groups ###
library(RColorBrewer)

# removes all glob variables from env
rm(list = ls())

# set wd
# setwd('/Users/sanjavickovic/Desktop/morphoSPOT/3dst_repo/3dst/grouping')
setwd('/Users/svickovi/Library/Mobile Documents/com~apple~CloudDocs/Desktop/morphoSPOT/3dst_repo/3dst/grouping')
# Which samples do you want to use? 
norm_samples = "RA5"

# Where are you norm expression R files located? 
path_samples = "../data/"

# Set output directory that will contain all cell_typing pdf plots and output gene files
path_output = "../../../plot_outputs/"

# Load data as R objects for selected RA patient
anns = list.files(pattern = glob2rx(paste0(norm_samples, "*annotations.txt")), path = path_samples)
norm_mat_loaded_in_env = ""
section_labels = ""
all_ann_inf = ""
for (i in 1:length(anns)){
  load(paste0(path_samples, norm_samples, "_", i, "_selected_adjusted_spots_3D_manual_app"))
  s = get(paste0("s",i))
  assign(paste0("mat"), read.delim(paste0(path_samples, norm_samples, "_", i, "_annotations.txt"), header = T)) # files found on SCP

  # makes sure x_y names are the same
  mat$x_y = paste0("X",i,"_", mat$x_y)
  row.names(mat) = mat$x_y
  
  # subset so both have same spots 
  s = s[row.names(s) %in% row.names(mat),]
  mat = mat[row.names(mat) %in% row.names(s),]
  
  # sort by row.names
  s = s[sort(row.names(s)),]
  mat = mat[sort(row.names(mat)),]
  
  # merge into one df
  if (unique(row.names(s) == row.names(mat)) == T){
    mat_xy = cbind(s,mat)
  }
  
  # Create data set with only infiltrates in matrix
  mat_xy = mat_xy[mat_xy$value %in% grep(mat$value, pattern ="Infiltrate", value = T) ,]
  mat_xy$value ="Infiltrate"
  
  # Create example data and transform into projected coordinate system
  x = mat_xy$V1
  y = mat_xy$V2
  chc <- hclust(dist(data.frame(x=x,y=y)), method="ward.D2")
  plot(chc)
  # Distance with a 13 cluster threshold threshold  
  # if (i == 4) {chc.d40 <- cutree(chc, h = 0.5)}
  # else  chc.d40 <- cutree(chc, h = 0.6) 
  chc.d40 <- cutree(chc, h = 0.5) 
  print(paste0("Found ", max(chc.d40), " infiltrate clusters.."))
  
  # Join results to meuse sp points
  xy <- data.frame(data.frame(rownames=mat_xy$x_y, x=x,y=y), Clust=chc.d40)
  
  # Plot resultsxy
  # Add extra space to right of plot area; change clipping to figure
  pdf(paste0(norm_samples, "_",i, "_inf_spatial_groups.pdf"), width = 5, height = 5)
  par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
  set3 <- colorRampPalette(c(brewer.pal('Set3',n=12),c(brewer.pal('Set1',n=9))))
  plot(xy$x, xy$y,col=setNames(set3(max(xy$Clust)), unique(xy$Clust))[xy$Clust], pch=19, xlab = "x coordinates", ylab = "y coordinates")
  text(xy$x, xy$y,labels=xy$rownames, cex=0.3, font=2)
  title(main="Infiltrate Clusters")
  # so turn off clipping:
  legend("topright", inset=c(-0.3,0), legend=paste("Cluster", 1:max(chc.d40),sep=""),
         col=setNames(set3(max(xy$Clust)), unique(xy$Clust)), pch=16, bg="white",cex = 0.5)
  dev.off()
  
  # output *all_infs* files 
  xy$Clust = paste0("Inf",xy$Clust)
  outxy = as.matrix(xy[,4])
  row.names(outxy) = xy$rownames
  colnames(outxy) = "infiltrate"
  write.table(outxy, file = paste0("../data/", norm_samples,"_",i,"_all_inf.csv"), row.names = T, col.names = T, quote = F, sep = ",")
}

# setwd(path_output)
fl = list.files(pattern = glob2rx(paste0(norm_samples, "*_all_inf.csv")), path = path_samples)

infs1 = read.csv(paste0(path_samples,fl[1]), row.names = 1, header = T)
infs2 = read.csv(paste0(path_samples,fl[2]), row.names = 1, header = T)
infs3 = read.csv(paste0(path_samples,fl[3]), row.names = 1, header = T)
infs4 = read.csv(paste0(path_samples,fl[4]), row.names = 1, header = T)

infs_1 = c(row.names(subset(infs1, infs1$infiltrate == "Inf1")),
           row.names(subset(infs2, infs2$infiltrate == "Inf1")),
           row.names(subset(infs3, infs3$infiltrate == "Inf1")),
           row.names(subset(infs4, infs4$infiltrate == "Inf2")))
infs_exp1 = matrix(nrow = length(infs_1), ncol = 1)
infs_exp1[,1] = "Inf1"
row.names(infs_exp1) = infs_1
colnames(infs_exp1) = "Inf"


infs_2 = c(row.names(subset(infs1, infs1$infiltrate == "Inf2")),
           row.names(subset(infs2, infs2$infiltrate == "Inf2")),
           row.names(subset(infs3, infs3$infiltrate == "Inf2")),
           row.names(subset(infs4, infs4$infiltrate == "Inf3")),
           row.names(subset(infs4, infs4$infiltrate == "Inf4")))
infs_exp2 = matrix(nrow = length(infs_2), ncol = 1)
infs_exp2[,1] = "Inf2"
row.names(infs_exp2) = infs_2
colnames(infs_exp2) = "Inf"

infs_3 = c(row.names(subset(infs1, infs1$infiltrate == "Inf3")),
           row.names(subset(infs2, infs2$infiltrate == "Inf4")),
           row.names(subset(infs3, infs3$infiltrate == "Inf4")),
           row.names(subset(infs4, infs4$infiltrate == "Inf5")))
infs_exp3= matrix(nrow = length(infs_3), ncol = 1)
infs_exp3[,1] = "Inf3"
row.names(infs_exp3) = infs_3
colnames(infs_exp3) = "Inf"

infs_4 = c(row.names(subset(infs1, infs1$infiltrate == "Inf4")),
           row.names(subset(infs2, infs2$infiltrate == "Inf3")),
           row.names(subset(infs3, infs3$infiltrate == "Inf6")),
           row.names(subset(infs4, infs4$infiltrate == "Inf7")))
infs_exp4= matrix(nrow = length(infs_4), ncol = 1)
infs_exp4[,1] = "Inf4"
row.names(infs_exp4) = infs_4
colnames(infs_exp4) = "Inf"

infs_5 = c(row.names(subset(infs1, infs1$infiltrate == "Inf5")),
           row.names(subset(infs1, infs1$infiltrate == "Inf6")),
           row.names(subset(infs2, infs2$infiltrate == "Inf6")),
           row.names(subset(infs3, infs3$infiltrate == "Inf3")),
           row.names(subset(infs3, infs3$infiltrate == "Inf5")),
           row.names(subset(infs4, infs4$infiltrate == "Inf1")),
           row.names(subset(infs4, infs4$infiltrate == "Inf6")))
infs_exp5= matrix(nrow = length(infs_5), ncol = 1)
infs_exp5[,1] = "Inf5"
row.names(infs_exp5) = infs_5
colnames(infs_exp5) = "Inf"

infs_6 = c(row.names(subset(infs1, infs1$infiltrate == "Inf7")),
           row.names(subset(infs2, infs2$infiltrate == "Inf5")),
           row.names(subset(infs3, infs3$infiltrate == "Inf8")),
           row.names(subset(infs4, infs4$infiltrate == "Inf8")),
           row.names(subset(infs4, infs4$infiltrate == "Inf10")))
infs_exp6= matrix(nrow = length(infs_6), ncol = 1)
infs_exp6[,1] = "Inf6"
row.names(infs_exp6) = infs_6
colnames(infs_exp6) = "Inf"

infs_7 = c(row.names(subset(infs1, infs1$infiltrate == "Inf8")),
           row.names(subset(infs2, infs2$infiltrate == "Inf8")),
           row.names(subset(infs2, infs2$infiltrate == "Inf9")),
           row.names(subset(infs3, infs3$infiltrate == "Inf7")),
           row.names(subset(infs4, infs4$infiltrate == "Inf9")))
infs_exp7= matrix(nrow = length(infs_7), ncol = 1)
infs_exp7[,1] = "Inf7"
row.names(infs_exp7) = infs_7
colnames(infs_exp7) = "Inf"

infs_8 = c(row.names(subset(infs1, infs1$infiltrate == "Inf9")),
           row.names(subset(infs1, infs1$infiltrate == "Inf10")),
           row.names(subset(infs2, infs2$infiltrate == "Inf7")),
           row.names(subset(infs3, infs3$infiltrate == "Inf9")),
           row.names(subset(infs4, infs4$infiltrate == "Inf11")),
           row.names(subset(infs4, infs4$infiltrate == "Inf12")))
infs_exp8= matrix(nrow = length(infs_8), ncol = 1)
infs_exp8[,1] = "Inf8"
row.names(infs_exp8) = infs_8
colnames(infs_exp8) = "Inf"

infs_9 = c(row.names(subset(infs1, infs1$infiltrate == "Inf11")),
           row.names(subset(infs2, infs2$infiltrate == "Inf10")),
           row.names(subset(infs2, infs2$infiltrate == "Inf11")),
           row.names(subset(infs3, infs3$infiltrate == "Inf10")),
           row.names(subset(infs4, infs4$infiltrate == "Inf13")))
infs_exp9= matrix(nrow = length(infs_9), ncol = 1)
infs_exp9[,1] = "Inf9"
row.names(infs_exp9) = infs_9
colnames(infs_exp9) = "Inf"

inf_all_ra3 = rbind(infs_exp1,infs_exp2,infs_exp3,
                    infs_exp4,infs_exp5,infs_exp6,
                    infs_exp7,infs_exp8,infs_exp9)


length(unique(row.names(inf_all_ra3))) == nrow(inf_all_ra3)
write.table(inf_all_ra3, file = "../data/RA3_zstack_Infs.csv", col.names  = T, sep =",")


# check1= inf_all_ra3[grep(row.names(inf_all_ra3),pattern= "X1", value = T),]
# x=sapply(strsplit(names(check1), "_"),"[[",2)
# y=sapply(strsplit(names(check1), "_"),"[[",3)
# set3 <- colorRampPalette(c(brewer.pal('Set3',n=12),c(brewer.pal('Set1',n=9))))
# plot(x, y,col=setNames(set3(9), unique(check1))[check1], pch=19, xlab = "x coordinates", ylab = "y coordinates")
#
# check1= inf_all_ra3[grep(row.names(inf_all_ra3),pattern= "X4", value = T),]
# x=sapply(strsplit(names(check1), "_"),"[[",2)
# y=sapply(strsplit(names(check1), "_"),"[[",3)
# set3 <- colorRampPalette(c(brewer.pal('Set3',n=12),c(brewer.pal('Set1',n=9))))
# plot(x, y,col=setNames(set3(9), unique(check1))[check1], pch=19, xlab = "x coordinates", ylab = "y coordinates")
# 
# fl = list.files(pattern = glob2rx(paste0(norm_samples, "*_all_inf.csv")), path = path_samples)
# 
# infs1 = read.csv(paste0(path_samples,fl[1]), row.names = 1, header = T)
# infs2 = read.csv(paste0(path_samples,fl[2]), row.names = 1, header = T)
# infs3 = read.csv(paste0(path_samples,fl[3]), row.names = 1, header = T)
# infs4 = read.csv(paste0(path_samples,fl[4]), row.names = 1, header = T)
# infs5 = read.csv(paste0(path_samples,fl[5]), row.names = 1, header = T)
# 
# infs_1 = c(row.names(subset(infs1, infs1$infiltrate == "Inf1")),
#            row.names(subset(infs1, infs1$infiltrate == "Inf2")),
#             row.names(subset(infs2, infs2$infiltrate == "Inf1")),
#             row.names(subset(infs3, infs3$infiltrate == "Inf1")),
#             row.names(subset(infs4, infs4$infiltrate == "Inf4")),
#            row.names(subset(infs5, infs5$infiltrate == "Inf4")),
#            row.names(subset(infs5, infs5$infiltrate == "Inf1")))
# infs_1 = infs_1[infs_1 != "X5_19_25"]
# infs_exp1 = matrix(nrow = length(infs_1), ncol = 1)
# infs_exp1[,1] = "Inf1"
# row.names(infs_exp1) = infs_1
# colnames(infs_exp1) = "Inf"
# 
# infs_2 = c(row.names(subset(infs1, infs1$infiltrate == "Inf5")),
#             row.names(subset(infs2, infs2$infiltrate == "Inf3")),
#             row.names(subset(infs3, infs3$infiltrate == "Inf5")),
#             row.names(subset(infs4, infs4$infiltrate == "Inf3")),
#             row.names(subset(infs5, infs5$infiltrate == "Inf5")))
# infs_exp2 = matrix(nrow = length(infs_2), ncol = 1)
# infs_exp2[,1] = "Inf2"
# row.names(infs_exp2) = infs_2
# colnames(infs_exp2) = "Inf"
# 
# infs_3 = c(row.names(subset(infs1, infs1$infiltrate == "Inf6")),
#             row.names(subset(infs2, infs2$infiltrate == "Inf4")),
#             row.names(subset(infs3, infs3$infiltrate == "Inf6")),
#             row.names(subset(infs4, infs4$infiltrate == "Inf5")),
#            row.names(subset(infs5, infs5$infiltrate == "Inf6")))
# infs_exp3= matrix(nrow = length(infs_3), ncol = 1)
# infs_exp3[,1] = "Inf3"
# row.names(infs_exp3) = infs_3
# 
# inf_all_ra4 = rbind(infs_exp1,infs_exp2,infs_exp3)
# inf_all_ra4 = as.matrix(inf_all_ra4[!row.names(inf_all_ra4) %in% c("X5_45_41", "X4_44_44", "X3_3_19", "X2_45_45", "X1_5_31", "X1_2_42"),])
# colnames(inf_all_ra4) = "Inf"
# write.table(inf_all_ra4, file = "../data/RA4_zstack_Infs.csv", col.names  = T, sep =",")
# #
# check1= inf_all_ra4[grep(row.names(inf_all_ra4),pattern= "X1", value = T),]
# x=sapply(strsplit(names(check1), "_"),"[[",2)
# y=sapply(strsplit(names(check1), "_"),"[[",3)
# set3 <- colorRampPalette(c(brewer.pal('Set3',n=12),c(brewer.pal('Set1',n=9))))
# plot(x, y,col=setNames(set3(9), unique(check1))[check1], pch=19, xlab = "x coordinates", ylab = "y coordinates")

# 
# fl = list.files(pattern = glob2rx(paste0(norm_samples, "*_all_inf.csv")), path = path_samples)
# 
# # infs1 = read.csv(paste0(path_samples,fl[1]), row.names = 1, header = T)
# infs2 = read.csv(paste0(path_samples,fl[1]), row.names = 1, header = T)
# infs3 = read.csv(paste0(path_samples,fl[2]), row.names = 1, header = T)
# infs4 = read.csv(paste0(path_samples,fl[3]), row.names = 1, header = T)
# 
# infs_1 = c(row.names(subset(infs2, infs2$infiltrate == "Inf1")),
#            row.names(subset(infs3, infs3$infiltrate == "Inf1")),
#            row.names(subset(infs4, infs4$infiltrate == "Inf1")))
# infs_exp1 = matrix(nrow = length(infs_1), ncol = 1)
# infs_exp1[,1] = "Inf1"
# row.names(infs_exp1) = infs_1
# colnames(infs_exp1) = "Inf"
# 
# 
# infs_2 = c(row.names(subset(infs2, infs2$infiltrate == "Inf2")),
#            row.names(subset(infs2, infs2$infiltrate == "Inf5")),
#            row.names(subset(infs3, infs3$infiltrate == "Inf4")),
#            row.names(subset(infs4, infs4$infiltrate == "Inf5")),
#            row.names(subset(infs4, infs4$infiltrate == "Inf7")))
# infs_exp2 = matrix(nrow = length(infs_2), ncol = 1)
# infs_exp2[,1] = "Inf2"
# row.names(infs_exp2) = infs_2
# colnames(infs_exp2) = "Inf"
# 
# infs_3 = c(row.names(subset(infs2, infs2$infiltrate == "Inf3")),
#            row.names(subset(infs3, infs3$infiltrate == "Inf2")),
#            row.names(subset(infs4, infs4$infiltrate == "Inf4")))
# infs_exp3= matrix(nrow = length(infs_3), ncol = 1)
# infs_exp3[,1] = "Inf3"
# row.names(infs_exp3) = infs_3
# colnames(infs_exp3) = "Inf"
# 
# infs_4 = c(row.names(subset(infs2, infs2$infiltrate == "Inf4")),
#            row.names(subset(infs3, infs3$infiltrate == "Inf3")),
#            row.names(subset(infs4, infs4$infiltrate == "Inf6")))
# infs_exp4= matrix(nrow = length(infs_4), ncol = 1)
# infs_exp4[,1] = "Inf4"
# row.names(infs_exp4) = infs_4
# colnames(infs_exp4) = "Inf"
# 
# inf_all_ra5 = rbind(infs_exp1,infs_exp2,infs_exp3,infs_exp4)
# inf_all_ra5 = subset(inf_all_ra5,!row.names(inf_all_ra5) %in% c("X3_5_49", "X3_4_48",  "X1_38_52", "X2_36_50","X2_36_51", "X2_35_51", "X3_36_48","X3_37_47"))
# colnames(inf_all_ra5) = "Inf"
# 
# infs_5 = c("X1_38_52", "X2_36_50", "X2_36_51","X2_35_51", "X3_36_48", "X3_37_47")
# infs_exp5= matrix(nrow = length(infs_5), ncol = 1)
# infs_exp5[,1] = "Inf5"
# row.names(infs_exp5) = infs_5
# colnames(infs_exp5) = "Inf"
# inf_all_ra5 = rbind(inf_all_ra5, infs_exp5)
# inf_all_ra5
# 
# length(unique(row.names(inf_all_ra5))) == nrow(inf_all_ra5)
# write.table(inf_all_ra5, file = "../data/RA5_zstack_Infs.csv", col.names  = T, sep =",")
# 
# 
# check1= inf_all_ra5[grep(row.names(inf_all_ra5),pattern= "X3", value = T),]
# x=sapply(strsplit(names(check1), "_"),"[[",2)
# y=sapply(strsplit(names(check1), "_"),"[[",3)
# set3 <- colorRampPalette(c(brewer.pal('Set3',n=12),c(brewer.pal('Set1',n=9))))
# plot(x, y, col=setNames(set3(9), unique(check1))[check1], pch=19, xlab = "x coordinates", ylab = "y coordinates")
# par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
# set3 <- colorRampPalette(c(brewer.pal('Set3',n=12),c(brewer.pal('Set1',n=9))))
# plot(x, y, col=setNames(set3(9), unique(check1))[check1], pch=19, xlab = "x coordinates", ylab = "y coordinates")
# text(x, y,labels=check1, cex=0.3, font=2)
# check1[39]
