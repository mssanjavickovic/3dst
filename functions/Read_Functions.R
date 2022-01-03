### Read in all functions for DEG analysis and plotting in 3D
# DEG test function
bimod.diffExp.test=function(data1,data2,mygenes) {
  p_val=unlist(lapply(mygenes,function(x)diffLRT(as.numeric(data1[x,]),as.numeric(data2[x,]))))
  p_val[is.na(p_val)]=1
  avg_diff=unlist(lapply(mygenes,function(x)(expMean(as.numeric(data1[x,]))-expMean(as.numeric(data2[x,])))))
  toRet=data.frame(cbind(p_val,avg_diff),row.names=mygenes)
  toRet=toRet[order(toRet$p_val),]
  return(toRet)
}

# Ratio test function
diffLRT = function(x,y,xmin=1) {
  lrtX=bimodLikData(x)
  lrtY=bimodLikData(y)
  lrtZ=bimodLikData(c(x,y))
  lrt_diff=2*(lrtX+lrtY-lrtZ)
  return(1-pchisq(lrt_diff,3))
}

# Likelihood function
bimodLikData=function(x,xmin=0) {
  x1=x[x<=xmin]
  x2=x[x>xmin]
  xal=minmax(length(x2)/length(x),min=1e-5,max=(1-1e-5))
  likA=length(x1)*log(1-xal)
  mysd=sd(x2)
  if(length(x2)<2) {
    mysd=1
  }
  likB=length(x2)*log(xal)+sum(dnorm(x2,mean(x2),mysd,log=TRUE))
  return(likA+likB)
}

# Outersect function
outersect <- function(x, y) {
  sort(c(setdiff(x, y),
         setdiff(y, x)))
}

# Minmax function
minmax=function(data,min,max) {
  data2=data
  data2[data2>max]=max
  data2[data2<min]=min
  return(data2)
}

# Exp mean function 
expMean=function(x) {
  return(log(mean(exp(x)-1)+1))
}

# Combine all matrices together by row.names/genes
mbind<-function(...){
  Reduce( function(x,y){cbind(x,y[match(row.names(x),row.names(y)),])}, list(...) )
}

# mergebind<-function(x,y){merge(x,y,by="row.names",all.x=TRUE)}
# Cluster2.RA2
# merge(Cluster2.RA1, Cluster2.RA2)

# Function to calculate rowMeans per section 
rm_section = function (RA.norm.inf) {
  sapply(unique(sapply(strsplit(colnames(RA.norm.inf), split = "_"),"[[",1)), function(x) {
    if (!is.null(ncol(RA.norm.inf[,sapply(strsplit(colnames(RA.norm.inf), split = "_"),"[[",1) == x]))){
      rowMeans(RA.norm.inf[,sapply(strsplit(colnames(RA.norm.inf), split = "_"),"[[",1) == x])
    } else {
      RA.norm.inf[,sapply(strsplit(colnames(RA.norm.inf), split = "_"),"[[",1) == x]
    }
  })
}

rm_section_all = function (RA.norm.inf) {
  sapply(unique(sapply(strsplit(colnames(RA.norm.inf), split = "_"),"[[",1)), function(x) {
    if (!is.null(ncol(RA.norm.inf[,sapply(strsplit(colnames(RA.norm.inf), split = "_"),"[[",1) == x]))){
      RA.norm.inf[,sapply(strsplit(colnames(RA.norm.inf), split = "_"),"[[",1) == x]
    } else {
      RA.norm.inf[,sapply(strsplit(colnames(RA.norm.inf), split = "_"),"[[",1) == x]
    }
  })
}

sum_section = function (RA.norm.inf) {
  sapply(unique(sapply(strsplit(colnames(RA.norm.inf), split = "_"),"[[",1)), function(x) {
    if (!is.null(ncol(RA.norm.inf[,sapply(strsplit(colnames(RA.norm.inf), split = "_"),"[[",1) == x]))){
      rowSums(RA.norm.inf[,sapply(strsplit(colnames(RA.norm.inf), split = "_"),"[[",1) == x])
    } else {
      RA.norm.inf[,sapply(strsplit(colnames(RA.norm.inf), split = "_"),"[[",1) == x]
    }
  })
}

# Function to calculate sd for each row(gene) per section 
sd_section = function (RA.norm.inf) {
  library(matrixStats)
  sapply(unique(sapply(strsplit(colnames(RA.norm.inf), split = "_"),"[[",1)), function(x) {
    if (!is.null(ncol(RA.norm.inf[,sapply(strsplit(colnames(RA.norm.inf), split = "_"),"[[",1) == x]))){
      rowSds(RA.norm.inf[,sapply(strsplit(colnames(RA.norm.inf), split = "_"),"[[",1) == x])
    } else {
      rep(0, nrow(RA.norm.inf))
    }
  })
}

# Function to calculate size for each row(gene) (dataset) per section 
size_section = function (RA.norm.inf) {
  library(matrixStats)
  sapply(unique(sapply(strsplit(colnames(RA.norm.inf), split = "_"),"[[",1)), function(x) {
    if (!is.null(ncol(RA.norm.inf[,sapply(strsplit(colnames(RA.norm.inf), split = "_"),"[[",1) == x]))){
      ncol(RA.norm.inf[,sapply(strsplit(colnames(RA.norm.inf), split = "_"),"[[",1) == x])
    } else {
      1
    }
  })
}


# Run DPT on the whole inf matrix
run_dpt = function(inf_all){
  #inf_all = RA.norm[,colnames(RA.norm) %in% all_inf]
  data = t(as.matrix(inf_all))
  s = find.sigmas(data, step = 0.01, steps = 100, sample_rows = 500)
  ts <- Transitions(data, sigma = s)
  pt <- dpt(ts, branching = TRUE)
  ev <- eigen(as.matrix(ts@transitions), TRUE)$vectors
  dm <- as.data.frame(ev[, -1])
  colnames(dm) <- paste0('RA', seq_len(ncol(dm)))
  # rownames(dm) = rownames(data)
  # rownames(pt) = rownames(data)

  #df = DiffusionMap(data, density.norm = F)

  # ev <-  df@eigenvectors
  # dm <- as.data.frame(ev[, -1])
  # pc = prcomp(inf_all)
  # rot = pc[2]
  # rot = as.data.frame(rot)
  # 
  # dm = rot
  colnames(dm) <- paste0('RA', seq_len(ncol(dm)))
  rownames(dm) = colnames(inf_all)
  #scatter3D(dm$RA1, dm$RA2, dm$RA3)
  # eigs = pc$sdev^2
  
  
  
  # eigs[1]/sum(eigs)
  # eigs[2]/sum(eigs)
  # eigs[3]/sum(eigs)
  
  return(dm)
}

# Clean up some gene names; mat is a matrix with gene names as row names
genes_cleanup = function (mat){
  mat = mat[! row.names(mat) %in% grep(pattern="RPS", x= rownames(mat), value = T),]
  mat = mat[! row.names(mat) %in% grep(pattern="RPL", x= rownames(mat), value = T),]
  mat = mat[! row.names(mat) %in% grep(pattern="ambiguous", x= rownames(mat), value = T),]
  mat = mat[! row.names(mat) %in% grep(pattern="MTRNR", x= rownames(mat), value = T),]
  mat = mat[! row.names(mat) %in% grep(pattern="ENSG", x= rownames(mat), value = T),]
  mat = mat[! row.names(mat) %in% grep(pattern="MALAT1", x= rownames(mat), value = T),]
  mat = mat[! row.names(mat) %in% grep(pattern="B2M", x= rownames(mat), value = T),]
  return(mat)
}

# assign spots with specific cluster colors 
assign_cluster_colors = function(clusters, col.panel, k, m3_ref){
  
  col.panel = col.panel[str_order(col.panel)]
  col.inf.clusters = matrix(ncol = 2, nrow = 1) # before called col.clusters
  # col.panel = c("#f15f48", "#4ebceb", "#8a3795") # add more colors in f>5
  for (i in 1:k){
    # print(i)
    # print(names(clusters[clusters == i,]))
    print(paste0("Cluster", i, ": ", col.panel[i]))
    col.inf.clusters = rbind(col.inf.clusters, cbind(names(clusters[clusters == i,]), col.panel[i]))
  }
  col.inf.clusters = col.inf.clusters[-1,]
  row.names(col.inf.clusters) = col.inf.clusters[,1]
  col.inf.clusters = as.matrix(col.inf.clusters[,-1])
  col.inf.clusters = as.matrix(col.inf.clusters[match(colnames(m3_ref), rownames(col.inf.clusters)),])

  return (col.inf.clusters)
  
}

# Make Barplot for avg gene expression per cluster group 
avg_genes_barplot = function(all_bar, gen_names, cluster_names){
  
  #all_bar = RA.norm
  #### Make barplot of avg expression per interesting marker genes
  col.spatial.clusters = cluster_names
  # gen_names = c("CCL19","CXCL13","LTB","PRG4","MMP3","CD52","MS4A1","FN1","TYROBP") #, "CD79A", "CD79B", "TYROBP"
  col.scale = unique(col.spatial.clusters[,1])
  

  ### arrange matrix per cluster 
  # all_bar = mbind(mg1, mg2, mg3, mg4)
  all_bar = all_bar[row.names(all_bar) %in% gen_names,]
  all_bar = all_bar[,match(colnames(all_bar), rownames(col.spatial.clusters))]
  all_bar = all_bar[match(rownames(all_bar), gen_names),]
  
  # Avg expression and std error per cluster for selected genes
  clusters_rm = ""
  clusters_sd = ""
  cluster_names = ""
  counter = 1
  for (i in col.scale){
      #print(i)
      mat.cluster = all_bar[, colnames(all_bar) %in% names(col.spatial.clusters[col.spatial.clusters == i,])]
      clusters_rm_tmp = as.numeric(rowMeans(mat.cluster))
      clusters_sd_tmp = as.numeric(rowSds(mat.cluster))/sqrt(ncol(mat.cluster))
      clusters_rm = c(clusters_rm, clusters_rm_tmp)
      clusters_sd = c(clusters_sd, clusters_sd_tmp)
      cluster_names = c(cluster_names, rep(paste0("Cluster",counter), length(gen_names)))
      counter=counter+1
  }
  clusters_rm = as.numeric(clusters_rm[-1])
  clusters_sd = as.numeric(clusters_sd[-1])
  cluster_names = cluster_names[-1]
  gene_names = rep(row.names(all_bar), length(col.scale))
  
  #get wilcoxons values cluster1 vs all 
  wilc_p = ""
  for (j in rownames(all_bar)){
    print(j)
    if ((j == "PRG4") | (j == "MMP3")| (j == "FN1")| (j == "TYROBP")) met = "greater"
    else met = "less"

    
    counter = 1
    for (i in c(1:length(col.scale))){
      if (i == 1){
        
        #print(i)
        clusters_rm_tmp1 = ""
        clusters_rm_tmp2 = ""
        
        mat.cluster1 = all_bar[j, colnames(all_bar) %in% names(col.spatial.clusters[col.spatial.clusters == col.scale[1],])]
        tss = unique(sapply(strsplit(names(mat.cluster1),"_"),"[[",1))
        for (nm in tss){
          tmp = mat.cluster1[sapply(strsplit(names(mat.cluster1),"_"),"[[",1) == nm]
          clusters_rm_tmp1 = c(clusters_rm_tmp1,mean(tmp)) 
        }
        mat.cluster2 = all_bar[j, colnames(all_bar) %in% names(col.spatial.clusters[col.spatial.clusters == col.scale[1],])]
        tss = unique(sapply(strsplit(names(mat.cluster2),"_"),"[[",1))
        for (nm in tss){
          tmp = mat.cluster2[sapply(strsplit(names(mat.cluster2),"_"),"[[",1) == nm]
          clusters_rm_tmp2 = c(clusters_rm_tmp2,mean(tmp))
        }
        clusters_rm_tmp1 = as.numeric(clusters_rm_tmp1[-1])
        clusters_rm_tmp2 = as.numeric(clusters_rm_tmp2[-1])
        
        #print(clusters_rm_tmp1)
        #print(clusters_rm_tmp2)
        wilc_p = c(wilc_p, 1)

      } else {

        
        
        #print(i)
        clusters_rm_tmp1 = ""
        clusters_rm_tmp2 = ""
        
        mat.cluster1 = all_bar[j, colnames(all_bar) %in% names(col.spatial.clusters[col.spatial.clusters == col.scale[1],])]
        tss = unique(sapply(strsplit(names(mat.cluster1),"_"),"[[",1))
        for (nm in tss){
          tmp = mat.cluster1[sapply(strsplit(names(mat.cluster1),"_"),"[[",1) == nm]
          clusters_rm_tmp1 = c(clusters_rm_tmp1,mean(tmp)) 
        }
        mat.cluster2 = all_bar[j, colnames(all_bar) %in% names(col.spatial.clusters[col.spatial.clusters == col.scale[i],])]
        tss = unique(sapply(strsplit(names(mat.cluster2),"_"),"[[",1))
        for (nm in tss){
          tmp = mat.cluster2[sapply(strsplit(names(mat.cluster2),"_"),"[[",1) == nm]
          clusters_rm_tmp2 = c(clusters_rm_tmp2,mean(tmp))
        }
        clusters_rm_tmp1 = as.numeric(clusters_rm_tmp1[-1])
        clusters_rm_tmp2 = as.numeric(clusters_rm_tmp2[-1])
        
        print(clusters_rm_tmp1)
        print(clusters_rm_tmp2)
        
        if ((j == "CXCL13") | (j == "TYROBP") | (j == "CCL19") | (j == "LTB")){
          if (j == "CXCL13"){
            sv_cx1 = clusters_rm_tmp1}
          if (j == "TYROBP"){
            sv_ty1 = clusters_rm_tmp1}
          if (j == "CCL19"){
            sv_cc1 = clusters_rm_tmp1}
          if (j == "LTB"){
            sv_ltb1 = clusters_rm_tmp1}}
        #print(met)
        wilc_p = c(wilc_p, wilcox.test(clusters_rm_tmp1, clusters_rm_tmp2, paired = FALSE, alternative = met,exact=FALSE)$p.value)
        #print(wilcox.test(clusters_rm_tmp1, clusters_rm_tmp2, paired = FALSE, exact=FALSE)$p.value)
      }
    }
  }
  


  #sv_ty1 = as.numeric(sv_ty1)
  #sv_cx1 = as.numeric(sv_cx1)
  #sv_cc1 = as.numeric(sv_cc1)
  #sv_ltb1 = as.numeric(sv_ltb1)

  #print("testing cxcl13 vs tyrobp per cluster 1")
  #print(sv_cx1)
  #print(sv_ty1)
  #print(wilcox.test(sv_cx1, sv_ty1, paired = FALSE, alternative = 'less')$p.value)
  
  #print("testing cxcl13 vs ltb per cluster 1")
  #print(sv_cx1)
  #print(sv_ltb1)
  #print(wilcox.test(sv_cx1, sv_ltb1, paired = FALSE, alternative = 'less')$p.value)
  
  #print("testing ccl21 vs tyrobp per cluster 1")
  #print(sv_cc1)
  #print(sv_ltb1)
  #print(wilcox.test(sv_cc1, sv_ltb1, paired = FALSE, alternative = 'less')$p.value)

  
  
  wilc_p = as.numeric(wilc_p[-1])
  reordered_wilc_p = ""
  for (i in 1:length(col.scale)){
    reordered_wilc_p = c(reordered_wilc_p, wilc_p[seq(i, length(wilc_p), length(col.scale))])
  }
  reordered_wilc_p = as.numeric(reordered_wilc_p[-1])
  
  signs_p= ""
  for (i in reordered_wilc_p){
    if (i>0.05) signs_p=c(signs_p, "ns")
    if ((i<=0.05) & (i>0.01)) signs_p=c(signs_p, "*")
    if ((i<=0.01) & (i>0.001)) signs_p=c(signs_p, "**")
    if (i<=0.001) signs_p=c(signs_p, "***")
  }
  signs_p = signs_p[-1]
  group2 = cluster_names
  group1 = rep("Cluster1", length(signs_p))
  # Plot barplot
  
  myData = data.frame(clusters_rm, clusters_sd, cluster_names, gene_names,reordered_wilc_p,signs_p, group1, group2)
  print(myData)
  limits <- aes(ymax = clusters_rm + clusters_sd,
                ymin = clusters_rm - clusters_sd)
  
  p <- ggplot(data = myData, aes(x = factor(as.character(gene_names), levels = gen_names), y = clusters_rm, fill = cluster_names)) + 
    ylab(c("Avg Log(norm counts)")) + ylim(0, max(clusters_rm)+1)+ theme_bw() + 
    geom_bar(stat = "identity", position = position_dodge(width = 0.9)) + 
    scale_fill_manual(values = col.scale, aesthetics = "fill") + 
    geom_errorbar(mapping = limits, stat = "identity", position = position_dodge(width = 0.9), width = 0.3) + 
    #stat_pvalue_manual(reordered_wilc_p, label = signs_p, vjust = -1, bracket.nudge.y = 1) +
    theme(axis.line.y = element_line(colour = "black"), axis.title.x = element_blank(), axis.ticks.x= element_blank(), axis.title.y = element_text(), axis.text.x = element_text(angle = 0),panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())

  return(p)

}

# Make Barplot for avg gene expression per cluster group 
avg_genes_barplot_data = function(all_bar, gen_names, cluster_names){
  
  #all_bar = RA.norm
  #### Make barplot of avg expression per interesting marker genes
  col.spatial.clusters = cluster_names
  # gen_names = c("CCL19","CXCL13","LTB","PRG4","MMP3","CD52","MS4A1","FN1","TYROBP") #, "CD79A", "CD79B", "TYROBP"
  col.scale = unique(cluster_names[,1])
  
  ### arrange matrix per cluster 
  # all_bar = mbind(mg1, mg2, mg3, mg4)
  all_bar = all_bar[row.names(all_bar) %in% gen_names,]
  all_bar = all_bar[,match(colnames(all_bar), rownames(cluster_names))]
  all_bar = all_bar[match(rownames(all_bar), gen_names),]
  
  # Avg expression and std error per cluster for selected genes
  clusters_rm = ""
  clusters_sd = ""
  cluster_names = ""
  counter = 1
  for (i in col.scale){
    #print(i)
    mat.cluster = all_bar[, colnames(all_bar) %in% names(col.spatial.clusters[col.spatial.clusters == i,])]
    clusters_rm_tmp = as.numeric(rowMeans(mat.cluster))
    clusters_sd_tmp = as.numeric(rowSds(mat.cluster))/sqrt(ncol(mat.cluster))
    clusters_rm = c(clusters_rm, clusters_rm_tmp)
    clusters_sd = c(clusters_sd, clusters_sd_tmp)
    cluster_names = c(cluster_names, rep(paste0("Cluster",counter), length(gen_names)))
    counter=counter+1
  }
  clusters_rm = as.numeric(clusters_rm[-1])
  clusters_sd = as.numeric(clusters_sd[-1])
  cluster_names = cluster_names[-1]
  gene_names = rep(row.names(all_bar), length(col.scale))
  
  #get wilcoxons values cluster1 vs all 
  wilc_p = ""
  for (j in rownames(all_bar)){
    print(j)
    if ((j == "PRG4") | (j == "MMP3")| (j == "FN1")| (j == "TYROBP")) met = "less"
    else met = "greater"
    
    
    counter = 1
    for (i in c(1:length(col.scale))){
      if (i == 1){
        
        #print(i)
        clusters_rm_tmp1 = ""
        clusters_rm_tmp2 = ""
        
        mat.cluster1 = all_bar[j, colnames(all_bar) %in% names(col.spatial.clusters[col.spatial.clusters == col.scale[1],])]
        tss = unique(sapply(strsplit(names(mat.cluster1),"_"),"[[",1))
        for (nm in tss){
          tmp = mat.cluster1[sapply(strsplit(names(mat.cluster1),"_"),"[[",1) == nm]
          clusters_rm_tmp1 = c(clusters_rm_tmp1,mean(tmp)) 
        }
        mat.cluster2 = all_bar[j, colnames(all_bar) %in% names(col.spatial.clusters[col.spatial.clusters == col.scale[1],])]
        tss = unique(sapply(strsplit(names(mat.cluster2),"_"),"[[",1))
        for (nm in tss){
          tmp = mat.cluster2[sapply(strsplit(names(mat.cluster2),"_"),"[[",1) == nm]
          clusters_rm_tmp2 = c(clusters_rm_tmp2,mean(tmp))
        }
        clusters_rm_tmp1 = as.numeric(clusters_rm_tmp1[-1])
        clusters_rm_tmp2 = as.numeric(clusters_rm_tmp2[-1])
        
        #print(clusters_rm_tmp1)
        #print(clusters_rm_tmp2)
        wilc_p = c(wilc_p, 1)
        
      } else {
        
        
        
        #print(i)
        clusters_rm_tmp1 = ""
        clusters_rm_tmp2 = ""
        
        mat.cluster1 = all_bar[j, colnames(all_bar) %in% names(col.spatial.clusters[col.spatial.clusters == col.scale[1],])]
        tss = unique(sapply(strsplit(names(mat.cluster1),"_"),"[[",1))
        for (nm in tss){
          tmp = mat.cluster1[sapply(strsplit(names(mat.cluster1),"_"),"[[",1) == nm]
          clusters_rm_tmp1 = c(clusters_rm_tmp1,mean(tmp)) 
        }
        mat.cluster2 = all_bar[j, colnames(all_bar) %in% names(col.spatial.clusters[col.spatial.clusters == col.scale[i],])]
        tss = unique(sapply(strsplit(names(mat.cluster2),"_"),"[[",1))
        for (nm in tss){
          tmp = mat.cluster2[sapply(strsplit(names(mat.cluster2),"_"),"[[",1) == nm]
          clusters_rm_tmp2 = c(clusters_rm_tmp2,mean(tmp))
        }
        clusters_rm_tmp1 = as.numeric(clusters_rm_tmp1[-1])
        clusters_rm_tmp2 = as.numeric(clusters_rm_tmp2[-1])
        
        #print(clusters_rm_tmp1)
        #print(clusters_rm_tmp2)
        
        if ((j == "CXCL13") | (j == "TYROBP") | (j == "CCL19") | (j == "LTB")){
          if (j == "CXCL13"){
            sv_cx1 = clusters_rm_tmp1}
          if (j == "TYROBP"){
            sv_ty1 = clusters_rm_tmp1}
          if (j == "CCL19"){
            sv_cc1 = clusters_rm_tmp1}
          if (j == "LTB"){
            sv_ltb1 = clusters_rm_tmp1}}
        #print(met)
        wilc_p = c(wilc_p, wilcox.test(clusters_rm_tmp1, clusters_rm_tmp2, paired = FALSE, alternative = met,exact=FALSE)$p.value)
        #print(wilcox.test(clusters_rm_tmp1, clusters_rm_tmp2, paired = FALSE, exact=FALSE)$p.value)
      }
    }
  }
  
  
  
  #sv_ty1 = as.numeric(sv_ty1)
  #sv_cx1 = as.numeric(sv_cx1)
  #sv_cc1 = as.numeric(sv_cc1)
  #sv_ltb1 = as.numeric(sv_ltb1)
  
  #print("testing cxcl13 vs tyrobp per cluster 1")
  #print(sv_cx1)
  #print(sv_ty1)
  #print(wilcox.test(sv_cx1, sv_ty1, paired = FALSE, alternative = 'less')$p.value)
  
  #print("testing cxcl13 vs ltb per cluster 1")
  #print(sv_cx1)
  #print(sv_ltb1)
  #print(wilcox.test(sv_cx1, sv_ltb1, paired = FALSE, alternative = 'less')$p.value)
  
  #print("testing ccl21 vs tyrobp per cluster 1")
  #print(sv_cc1)
  #print(sv_ltb1)
  #print(wilcox.test(sv_cc1, sv_ltb1, paired = FALSE, alternative = 'less')$p.value)
  
  
  
  wilc_p = as.numeric(wilc_p[-1])
  reordered_wilc_p = ""
  for (i in 1:length(col.scale)){
    reordered_wilc_p = c(reordered_wilc_p, wilc_p[seq(i, length(wilc_p), length(col.scale))])
  }
  reordered_wilc_p = as.numeric(reordered_wilc_p[-1])
  
  signs_p= ""
  for (i in reordered_wilc_p){
    if (i>0.05) signs_p=c(signs_p, "ns")
    if ((i<=0.05) & (i>0.01)) signs_p=c(signs_p, "*")
    if ((i<=0.01) & (i>0.001)) signs_p=c(signs_p, "**")
    if (i<=0.001) signs_p=c(signs_p, "***")
  }
  signs_p = signs_p[-1]
  group2 = cluster_names
  group1 = rep("Cluster1", length(signs_p))
  # Plot barplot
  
  myData = data.frame(clusters_rm, clusters_sd, cluster_names, gene_names,reordered_wilc_p,signs_p, group1, group2)

  return(myData)
  
}

# Check how many inf and clusters correspond
inf_per_cluster = function(all_bar, all_inf, col.spatial.clusters){
  
  col.scale = unique(col.spatial.clusters[,1])

  counter = 1
   for (i in col.scale){
     
      mat.cluster = all_bar[, colnames(all_bar) %in% names(col.spatial.clusters[col.spatial.clusters == i,])]
      mat.cluster = mat.cluster[, colnames(mat.cluster) %in% all_inf]
    
      if (is.null(ncol(mat.cluster)) == T){
        numbers_cl = 1} 
      else {numbers_cl=ncol(mat.cluster)}
     
      print(paste0("Cluster: ",counter," contains ", round(100*numbers_cl/length(all_inf),2) ,"% of annotated infiltrates"))
      counter = counter +1 
  }
}

# Make regular heatmap with single or double annotations
plot_heatmap = function(m3_ref, col.inf.clusters){
 #  source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")
  
  # m3_ref = m3_ref[genes,]
  quantile.range <- quantile(m3_ref, probs = seq(0, 1, 0.05))
  palette.breaks <- seq(quantile.range["0%"], quantile.range["95%"], 0.1)
  
  # use http://colorbrewer2.org/ to find optimal divergent color palette (or set own)
  col.pal = colorRampPalette(c("black","purple","darkorchid1","gold1"))
  color.palette  <- col.pal(length(palette.breaks) - 1)
  
  if (ncol(col.inf.clusters) == 1){
  
  col.inf.clusters2 = as.matrix(col.inf.clusters[str_order(col.inf.clusters[,1]),])
  
  m3_ref = m3_ref[,match(row.names(col.inf.clusters2), colnames(m3_ref))]
  df_col = data.frame(col.inf.clusters2[match(row.names(col.inf.clusters2), colnames(m3_ref)),])
  colnames(df_col) = c("Clusters")
    
  reps=""
  for (i in 1:length(unique(names(table(df_col))))){
    reps = c(reps, rep(paste0("Cluster",i), table(df_col)[i]))
  }
  reps = reps[-1]
  
  annotdf <- data.frame(row.names = rownames(df_col), Clusters = reps)
  newCols <- colorRampPalette(unique(df_col$Clusters))
  mycolors <- newCols(length(unique(annotdf$Clusters)))
  names(mycolors) <- unique(annotdf$Clusters)
  mycolors <- list(Clusters = mycolors)
  
  heat = pheatmap(m3_ref, show_colnames = F , cluster_rows = T, cluster_cols = F, 
           col=color.palette, breaks = palette.breaks,clustering_distance_cols = "euclidean",
            scale="none", trace = "none",density.info = "none", cexRow = 0.05,
           annotation_col = annotdf, annotation_colors = mycolors)
  } else {
    
    # col.inf.clusters = double_ann
    col.inf.clusters2 = as.matrix(col.inf.clusters[str_order(col.inf.clusters[,2]),])
    
    m3_ref = m3_ref[,match(row.names(col.inf.clusters2), colnames(m3_ref))]
    df_col = data.frame(col.inf.clusters2[match(row.names(col.inf.clusters2), colnames(m3_ref)),])
    colnames(df_col) = c("Infiltrates", "Clusters")
    
    df_col_cope = df_col
    counter = 1
    for (i in names(table(df_col[,2]))){
      df_col_cope[,2] = str_replace(df_col_cope[,2], i, paste0("Cluster", counter))
      counter = counter + 1
    }

    for (i in names(table(df_col[,1]))){
      if (i == "red"){
        df_col_cope[,1] = str_replace(df_col_cope[,1], i, paste0("Infiltrate"))}
      else df_col_cope[,1] = str_replace(df_col_cope[,1], i, paste0("Other"))
    } 

    annotdf <- data.frame(row.names = rownames(df_col_cope), Clusters = df_col_cope$Clusters, Infiltrates = df_col_cope$Infiltrates)
    newCols_clus <- colorRampPalette(unique(df_col$Clusters))
    mycolors_clus <- newCols_clus(length(unique(annotdf$Clusters)))
    names(mycolors_clus) <- unique(annotdf$Clusters)
        
    newCols_inf <- colorRampPalette(unique(df_col$Infiltrates))
    mycolors_inf <- newCols_inf(length(unique(annotdf$Infiltrates)))
    names(mycolors_inf) <- unique(annotdf$Infiltrates)
    
    mycolors <- list(Clusters = mycolors_clus, Infiltrates = mycolors_inf)
    
    heat <- pheatmap(m3_ref, show_colnames = F , cluster_rows = T, cluster_cols = F, 
             col=color.palette, breaks = palette.breaks,clustering_distance_cols = "euclidean",
             scale="none", trace = "none",density.info = "none", cexRow = 0.1,
             annotation_col = annotdf, annotation_colors = mycolors, fontsize_row = 7)
    
    return(heat)
    
  }
}

# Run analysis of annotated infiltrates 
run_inf_analysis = function(inf_all, k, norm_samples, path_output){
  
  # Run DPT on the whole inf matrix
  print("Running PC on infiltrates ...")
  #inf_all = RA.norm[,colnames(RA.norm) %in% all_inf]
  #k= 3
  dm = run_dpt(inf_all)
  
  # Do hierarchical clustering
  print("Running clustering on infiltrates ...")
  hc = hclust(dist(as.matrix(dm$RA1, dm$RA2, dm$RA3)),method = "ward.D2")
 
  k = k
  clusters = cutree(hc, k = k) # 2 clusters
  names(clusters) = rownames(dm) 
  clusters = as.matrix(clusters)
  
  # Assign colors to clusters
  col.panel = c("#4ebceb","#8a3795","#f15f48")
  col.inf.clusters = assign_cluster_colors(clusters, col.panel = col.panel, k = k, inf_all)
  col.inf.clusters = as.matrix(col.inf.clusters[match(row.names(col.inf.clusters), row.names(dm)),])
  
  ### Do DEG analysis based on Nature Biotechnology 33, 495–502 (2015) doi:10.1038/nbt.3192
  inf_all = genes_cleanup(inf_all)
  m3_ref = inf_all # Precaution
  
  # Split into Inf cluster groups
  for (i in 1:k){
    assign(paste0("mg", i), m3_ref[,clusters == i])
  }
  
  # Do bimodial test where degg[,1] is the p-value and degg[,2] is the relative difference between the groups
  genes = ""
  for (i in 1:k){
    print(paste0("Peforming DE analysis for inf cluster: ", i))
    degg <- bimod.diffExp.test(get(paste0("mg", i)), m3_ref[,!colnames(m3_ref) %in% colnames(get(paste0("mg", i)))], row.names(get(paste0("mg", i))))
    assign(paste0("clus", i, ".genes"), rownames(subset(degg, degg[,1] < 0.001 & degg[,2] > 0.5)))
    assign(paste0("clus", i, ".genes.values"), degg[get(paste0("clus", i, ".genes")),])
    genes = c(genes, get(paste0("clus", i, ".genes")))
    write.table(get(paste0("clus", i, ".genes")), file = paste0(path_output, "DEGs.", "Inf.Cluster",i,"_subsetted.", norm_samples, ".txt"), quote = F, row.names = T)
    write.table(get(paste0("clus", i, ".genes.values")), file = paste0(path_output, "DEGs.", "Inf.Cluster",i,".", norm_samples, ".csv"), quote = F, row.names = T, col.names = T, sep =",")
  }
  
  # Clean up gene names (remove house keepers)
  genes = genes[-1]
  house.keeping.genes = row.names(read.delim("../data/house_keeping_final.txt", row.names = 1)) # found in ./data
  genes <- intersect(genes, outersect(genes,row.names(house.keeping.genes)))
  
  # Plot Heatmap
  # pdf("PCAs_Inf.pdf")
  print("Generating plots on infiltrates ...")
  pdf(paste0(path_output,"PCA_Infiltrates_", norm_samples, ".pdf"))
  scatter3D(dm$RA1, dm$RA2, dm$RA3, main = "Infiltrate PCA", pch=16, col = col.inf.clusters, colkey = list(plot = FALSE))
  write.table(dm[,1:3], file = paste0("/Users/sanjavickovic/Desktop/morphoSPOT/Manuscript_Nov2021/figshare/run_inf_analysis_pca_", norm_samples, ".csv"), sep = ",", quote = F)
  dev.off()
  hc_plot = plot(hc, labels = F, sub="", xlab = "", main = "Infiltrate Dendogram")
  rect.hclust(hc , k = k, border = col.panel[1:k])
  # dev.off()
  
  # pdf("Heatmap_Inf.pdf")
  col.inf.clusters = as.matrix(col.inf.clusters[match(colnames(m3_ref),row.names(col.inf.clusters)),])
  
  heat = plot_heatmap(m3_ref[genes,], col.inf.clusters)
  mmm = rbind(t(as.matrix(col.inf.clusters)), m3_ref[genes,])
  write.table(mmm, file = paste0("/Users/sanjavickovic/Desktop/morphoSPOT/Manuscript_Nov2021/figshare/run_inf_analysis_heatmap_", norm_samples, ".csv"), sep = ",", quote = F)
  
  pdf(paste0(path_output, "Heatmap_Infiltrates_", norm_samples, ".pdf"))
  genes1 = c("CXCL13","RAC2","LCP1", "MS4A1", "CD52", "FTL", "FTH1", "MZB1", "PIM2", "IGLL1", "CCL21", "CCL19", "XBP1", "FN1", "PCOLCE", "PRG4", "IL7R", "CRTAC1", "FOSB", "PRDX4", "IGLL5", "MMP3", "SSR4", "TYROBP", "CD79A", "CD79B", "LTB")
  add.flag(heat, kept.labels = genes1, repel.degree = 0)
  dev.off()
  
  return(col.inf.clusters)
}

# Run analysis of annotated infiltrates 
run_spatial_cluster_analysis = function(RA.norm, k, read_tsne_from_memory, norm_samples, path_output){
  
  f = k

  # Make tSNE clustering on all spots from the RA1 biopsy 
  all.exp.values.norm = RA.norm
  
  # this takes out per section variance
  # make a sce object
  sce <- SingleCellExperiment(list(logcounts=as.matrix(all.exp.values.norm)), # this is based on norm counts so object names is logcounts
                              colData=DataFrame(section=section_labels),
                              rowData=DataFrame(label=row.names(all.exp.values.norm)))

  # get variable genes
  print("Getting variable genes ...")
  dec3 <- modelGeneVar(sce, block=sce$section)
  
  # some QC plots to evaluate variance per section
  # per.block <- dec3$per.block
  # par(mfrow=c(3, 2))
  # for (i in seq_along(per.block)) { # check that variacle is correctly computed when taking into account sections/experiments
  #   decX <- per.block[[i]]
  #   plot(decX$mean, decX$total, xlab="Mean log-expression", 
  #        ylab="Variance", main=names(per.block)[i])
  #   curve(metadata(decX)$trend(x), col="blue", add=TRUE)
  # }
  # dev.off()

  # Get the top genes
  if (norm_samples == "RA1") co = 2000
  if (norm_samples == "RA2") co = 2000
  if (norm_samples == "RA3") co = 500
  if (norm_samples == "RA4") co = 100
  if (norm_samples == "RA5") co = 500
  if (norm_samples == "RA6") co = 200
  top.hvgs <- getTopHVGs(dec3, n=co)
  
  # Subset on only variable genes
  all.exp.values.norm = all.exp.values.norm[rownames(all.exp.values.norm) %in% top.hvgs,]
  all.exp.values.norm = genes_cleanup(all.exp.values.norm)
  
  # Color points based on only subsetting inflitrates
  tsne_inf = all.exp.values.norm[,colnames(all.exp.values.norm) %in% all_inf]
  tsne_rest = all.exp.values.norm[,!colnames(all.exp.values.norm) %in% all_inf]
  all_matrix = cbind(tsne_inf, tsne_rest)

  if (read_tsne_from_memory == "no"){
  
  # How many components are important? 
  print("How many PCs are important ...")
  per_out = permutationPA(t(all_matrix), B = 10, threshold = 0.05, verbose = TRUE, seed = NULL)
  pesr = per_out$r
  
  # library(Seurat)
  # all_seurat_obj = CreateSeuratObject(counts=all_matrix)
  # all_seurat_obj = ScaleData(all_seurat_obj)
  # all_seurat_obj <- RunPCA(object = all_seurat_obj,  npcs = 50, verbose = FALSE, features=row.names(all_matrix), assay= 'RNA')
  # all_seurat_obj <- JackStraw(object = all_seurat_obj, reduction = "pca", dims = 50, num.replicate = 10,  prop.freq = 0.1, verbose = TRUE)
  # all_seurat_obj <- ScoreJackStraw(all_seurat_obj, reduction = "pca", dims = 1:20, score.thresh = 0.01, do.plot = FALSE)
  # ElbowPlot(all_seurat_obj)

  print("Running tSNE ...")
  # Run tSNE again
  if (is.na(pesr)) {
    pesr = 20}
  tsne_out <- Rtsne(t(all_matrix), dims = 3, initial_dims = pesr, theta = 0.5, check_duplicates = FALSE, pca = TRUE, perplexity = 20, max_iter = 1000, verbose = TRUE) # Run TSNE; initial_dims = 10 in the manuscript
  
  # Make matrix of tsne data output
  tsne <- as.data.frame(tsne_out$Y) # saved matrix availalbe in ./data
  print("Saving tSNE matrix to data ...")
  rownames(tsne) = colnames(all_matrix)
  write.table(tsne, file = paste0("../data/tsne_", norm_samples, ".txt"), quote = F, sep ="\t", col.names = T, row.names=T)
  
  } else {tsne = read.delim(paste0("../data/tsne_", norm_samples, ".txt"))}
  
  t1 <- as.matrix(tsne[1])
  t2 <- as.matrix(tsne[2])
  t3 <- as.matrix(tsne[3])
  tsne_m <- cbind (t1,t2,t3)

  # Hierachical clustering on tsne matrix
  hc = hclust(dist(as.matrix(tsne_m[,1], tsne_m[,2])),method = "ward.D2")
  #f=4
  clusters = as.matrix(cutree(hc, k = f))
  col.panel = c("#FED8B1", "gray85", "gray76",  "gray61", "gray90", "gray99")
  row.names(clusters) = rownames(tsne_m)
  
  print("Plotting dendogram ... Please adjust k and rerun if needed")
  hc_plot = plot(hc, labels = F, sub="", xlab = "", main = "Spatial Cluster Dendogram")
  rect.hclust(hc, k = f, border = col.panel[1:f])

  # Assign colors to clusters
  print("Plotting tSNEs ...")
  
  col.spatial.clusters = assign_cluster_colors(clusters, col.panel, k = f, all_matrix)
  scatter3D(tsne_m[,1],tsne_m[,2],tsne_m[,3], phi = 0,theta = -180, axes=FALSE, ann=FALSE,pch=16, col=col.spatial.clusters)
  
  ann.inf.clusters = col.spatial.clusters
  ann.inf.clusters[row.names(ann.inf.clusters) %in% colnames(tsne_inf),] <- "red"
  ann.inf.clusters[!row.names(ann.inf.clusters) %in% colnames(tsne_inf),] <- "grey51"
  
  scatter3D(tsne_m[,1], tsne_m[,2],tsne_m[,3], phi = 0,theta = -180, axes=FALSE, ann=FALSE,pch=16, col=ann.inf.clusters)

  #Do DEG analysis
  m3_ref = all_matrix
  #print(head(all_matrix))
  
  # Split into spatial cluster groups
  for (i in 1:f){
    assign(paste0("mg", i), m3_ref[,clusters == i])
    secs = sapply(strsplit(colnames(get(paste0("mg", i))), "_"), "[[",1)
    
    for (j in unique(secs)){
      mtmp =  get(paste0("mg", i))[,grep(j, colnames(get(paste0("mg", i))), value = T)]
      x = sapply(strsplit(colnames(mtmp), "_"), "[[",2)
      y = sapply(strsplit(colnames(mtmp), "_"), "[[",3)
      img = paste0(norm_samples, "_", str_replace(j, "X",""))
      tmp1 = cbind(img, x, y, t(mtmp))
      colnames(tmp1) = c("Image_ID", "x", "y", row.names(mtmp))
      write.table(tmp1, file = paste0(path_output, "Cluster",i,".",str_replace(j, "X",""),".", norm_samples, ".csv"), sep =",", quote =F, row.names = F, col.names = T)
    }

  }

  # Do bimodial test where degg[,1] is the p-value and degg[,2] is the relative difference between the groups
  genes = ""
  for (i in 1:f){
    print(paste0("Peforming DE analysis for cluster: ", i))
    degg <- bimod.diffExp.test(get(paste0("mg", i)), m3_ref[,!colnames(m3_ref) %in% colnames(get(paste0("mg", i)))], row.names(get(paste0("mg", i))))
    assign(paste0("clus", i, ".genes"), rownames(subset(degg, degg[,1] < 0.001 & degg[,2] > 0.5)))
    assign(paste0("clus", i, ".genes.values"), degg[get(paste0("clus", i, ".genes")),])
    genes = c(genes, get(paste0("clus", i, ".genes")))
    write.table(get(paste0("clus", i, ".genes")), file = paste0(path_output, "DEGs.", "Spatial.Cluster",i,"_subsetted.", norm_samples, ".txt"), quote = F, row.names = T)
    write.table(get(paste0("clus", i, ".genes.values")), file = paste0(path_output, "DEGs.", "Cluster",i,".", norm_samples, ".csv"), quote = F, row.names = T, col.names = T, sep =",")
  }

  # Clean up gene names (remove house keepers)
  genes = genes[-1]
  house.keeping.genes = row.names(read.delim("../data/house_keeping_final.txt", row.names = 1)) # found in ./data
  genes <- intersect(genes, outersect(genes, row.names(house.keeping.genes)))

  # Plot Heatmap
  #pdf("Heatmap_tSNE.pdf")
  # make double annotation colors
  double_ann = mbind(ann.inf.clusters, col.spatial.clusters)
  heat = plot_heatmap(m3_ref[genes,], double_ann)
  pdf(paste0(path_output, "Heatmap_Cluster_", norm_samples, ".pdf"))
  add.flag(heat, kept.labels = c("TIMP1","CLU", "PIM2", "DERL3","HTRA1", "IGLL5", "CCL19","CXCL13","LTB","PRG4","MMP3","CD52","MS4A1","FN1", "TYROBP", "SSR4", "VCAM1", "MARCO", "SEPP1", "COL1A2", "CXCL12", "CXCR4", "CCL21","IL32", "RAC2", "PIM2"), repel.degree = 0)
  dev.off()
  #dev.off()
  
  # save files for figshare 
  m_output = rbind(t(double_ann), m3_ref[genes,])
  write.table(m_output, file = paste0("/Users/sanjavickovic/Desktop/morphoSPOT/Manuscript_Nov2021/figshare/Heatmap_Cluster_Matrix_", norm_samples, ".csv"), sep = ",", quote = F)
  
  return(col.spatial.clusters)
}

assign_col_inf_numbers = function(RA.norm,col.inf.clusters){
  
  ann.col.inf = matrix(nrow=ncol(RA.norm), ncol = 1)
  rownames(ann.col.inf) = colnames(RA.norm)
  for (sam in colnames(RA.norm)){
    if ((!sam %in% row.names(col.inf.clusters)) == TRUE) ann.col.inf[sam,1] = 4
    else {
    if (col.inf.clusters[sam,1] == "#f15f48") { ann.col.inf[sam,1] = 3}
    if (col.inf.clusters[sam,1] == "#4ebceb") { ann.col.inf[sam,1] = 1}
    if (col.inf.clusters[sam,1] == "#8a3795") { ann.col.inf[sam,1] = 2}
    }
  }
  return (ann.col.inf)
}

assign_col_cluster_numbers = function(col.spatial.clusters,all_inf){

  ann.all_inf = col.spatial.clusters
  ann.all_inf[row.names(ann.all_inf) %in% all_inf,] <- "#FF7F00"
  ann.all_inf[!row.names(ann.all_inf) %in% all_inf,] <- "grey51"
  ann.cluster = matrix(nrow=nrow(col.spatial.clusters), ncol = 1)
  rownames(ann.cluster) = rownames(col.spatial.clusters)
  double_ann = mbind(col.spatial.clusters, ann.all_inf)
  
  for (sam in rownames(double_ann)){
    if ((double_ann[sam,1] == "#FED8B1") & (double_ann[sam,2] == "#FF7F00") | (double_ann[sam,1] == "gray61") & (double_ann[sam,2] == "#FF7F00") | (double_ann[sam,1] == "gray76") & (double_ann[sam,2] == "#FF7F00") | (double_ann[sam,1] == "gray85") & (double_ann[sam,2] == "#FF7F00")){ ann.cluster[sam,1] = 5}
    else {
      if (double_ann[sam,1] == "#FED8B1") { ann.cluster[sam,1] = 4}
      if (double_ann[sam,1] == "gray61") { ann.cluster[sam,1] = 3}
      if (double_ann[sam,1] == "gray76") { ann.cluster[sam,1] = 2}
      if (double_ann[sam,1] == "gray85") { ann.cluster[sam,1] = 1}}
  }
  
  return(ann.cluster)
}

### Plot Heatmaps
### Plot gene specific stuff
plot.gene.3d.4 = function(sample, gene, mn1, mn2, mn3, mn4, s1, s2, s3, s4, x, y, transparency, min, max){
  
  for (i in c(1:4)) {
    if (i == "1") {
      spots = s1
      exp.values = mn1}
    if (i == "2") {
      spots = s2
      exp.values = mn2}
    if (i == "3") {
      spots = s3
      exp.values = mn3}
    if (i == "4") {
      spots = s4
      exp.values = mn4}
    
    genes.barcodes = exp.values[rownames(exp.values) == gene,]
    exp.values = exp.values[,colnames(exp.values) %in% rownames(spots)]
    exp.values= exp.values[rowSums(exp.values) != 0,]
    spots = spots[colnames(exp.values),]
    genes.barcodes.1 = cbind(spots, genes.barcodes)
    
    x1 = as.numeric(genes.barcodes.1[,1])
    x1 = x1[!is.na(x1)]
    y1 = as.numeric(genes.barcodes.1[,2])
    y1 = y1[!is.na(y1)]
    z = as.numeric(genes.barcodes.1[,3])
    dis = rep(0.25 * as.numeric(i), nrow(genes.barcodes.1))
    w = as.numeric(dis)

    s =  interp(x1,y1, z, nx = x, ny = y)
    if (i == "1") {mat.1 = s$z}
    if (i == "2") {mat.2 = s$z}
    if (i == "3") {mat.3 = s$z}
    if (i == "4") {mat.4 = s$z}}
  
  png(paste(sample,"_",gene,"_3d_clusters.png", sep=""))
  par(mar=c(2.1, 6.1, 15.1, 2.1), xpd=TRUE)

  image3D(z = 0, phi=25, theta = 10, colvar = mat.1, smooth = TRUE, alpha = transparency*1, box = FALSE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F)
  image3D(z = 4, phi=25,  theta = 10, colvar = mat.2, smooth = TRUE, alpha = transparency*1, add = TRUE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F)
  image3D(z = 8, phi=25,  theta = 10, colvar = mat.3, smooth = TRUE, alpha = transparency*1, add = TRUE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F)
  image3D(z = 11, phi=25,  theta =10, colvar = mat.4, smooth = TRUE, alpha = transparency*2, add = TRUE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F)
  dev.off()
}

plot.gene.3d.7 = function(sample, gene, mn1, mn2, mn3, mn4, mn5, mn6, mn7, s1, s2, s3, s4, s5, s6, s7,x, y, transparency, min, max){
  
  suppressWarnings({ for (i in c(1:7)) {
    if (i == "1") {
      spots = s1
      exp.values = mn1}
    if (i == "2") {
      spots = s2
      exp.values = mn2}
    if (i == "3") {
      spots = s3
      exp.values = mn3}
    if (i == "4") {
      spots = s4
      exp.values = mn4}
    if (i == "5") {
      spots = s5
      exp.values = mn5}
    if (i == "6") {
      spots = s6
      exp.values = mn6}
    if (i == "7") {
      spots = s7
      exp.values = mn7}
    
    genes.barcodes = exp.values[rownames(exp.values) == gene,]
    exp.values = exp.values[,colnames(exp.values) %in% rownames(spots)]
    exp.values= exp.values[rowSums(exp.values) != 0,]
    spots = spots[colnames(exp.values),]
    genes.barcodes.1 = cbind(spots, genes.barcodes)
    
    x1 = as.numeric(genes.barcodes.1[,1])
    x1 = x1[!is.na(x1)]
    y1 = as.numeric(genes.barcodes.1[,2])
    y1 = y1[!is.na(y1)]
    z = as.numeric(genes.barcodes.1[,3])
    dis = rep(0.25 * as.numeric(i), nrow(genes.barcodes.1))
    w = as.numeric(dis)
    
    s =  interp(x1,y1, z, nx = x, ny = y)
    if (i == "1") {mat.1 = s$z}
    if (i == "2") {mat.2 = s$z}
    if (i == "3") {mat.3 = s$z}
    if (i == "4") {mat.4 = s$z}
    if (i == "5") {mat.5 = s$z}
    if (i == "6") {mat.6 = s$z}
    if (i == "7") {mat.7 = s$z}
    }})
  
  png(paste(sample,"_",gene,"_3d_clusters.png", sep=""))
  par(mar=c(2.1, 6.1, 15.1, 2.1), xpd=TRUE)
  
  image3D(z = 0, phi=25, theta = 10, colvar = mat.1, smooth = TRUE, alpha = transparency*1, box = FALSE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F)
  image3D(z = 4, phi=25,  theta = 10, colvar = mat.2, smooth = TRUE, alpha = transparency*1, add = TRUE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F)
  image3D(z = 8, phi=25,  theta = 10, colvar = mat.3, smooth = TRUE, alpha = transparency*1, add = TRUE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F)
  image3D(z = 11, phi=25,  theta =10, colvar = mat.4, smooth = TRUE, alpha = transparency*1, add = TRUE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F)
  image3D(z = 14, phi=25,  theta = 10, colvar = mat.5, smooth = TRUE, alpha = transparency*1, add = TRUE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F)
  image3D(z = 17, phi=25,  theta = 10, colvar = mat.6, smooth = TRUE, alpha = transparency*2, add = TRUE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F)
  image3D(z = 21, phi=25,  theta =10, colvar = mat.7, smooth = TRUE, alpha = transparency*2, add = TRUE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F)
  
  dev.off()
}


plot.gene.3d.5 = function(sample, gene, m1, m2, m3, m4, m5, s1, s2, s3, s4, s5, x, y, transparency, min, max){
  
  for (i in c(1:5)) {
    if (i == "1") {
      spots = s1
      exp.values = m1
      }
    if (i == "2") {
      spots = s2
      exp.values = m2}
    if (i == "3") {
      spots = s3
      exp.values = m3}
    if (i == "4") {
      spots = s4
      exp.values = m4}
    if (i == "5") {
      spots = s5
      exp.values = m5}
    
    genes.barcodes = exp.values[rownames(exp.values) == gene,]
    exp.values = exp.values[,colnames(exp.values) %in% rownames(spots)]
    exp.values= exp.values[rowSums(exp.values) != 0,]
    spots = spots[colnames(exp.values),]
    genes.barcodes.1 = cbind(spots, genes.barcodes)
    
    x1 = as.numeric(genes.barcodes.1[,1])
    x1 = x1[!is.na(x1)]
    y1 = as.numeric(genes.barcodes.1[,2])
    y1 = y1[!is.na(y1)]
    z = as.numeric(genes.barcodes.1[,3])
    dis = rep(0.25 * as.numeric(i), nrow(genes.barcodes.1))
    w = as.numeric(dis)

        s =  interp(x1,y1, z, nx = x, ny = y)
    if (i == "1") {mat.1 = s$z}
    if (i == "2") {mat.2 = s$z}
    if (i == "3") {mat.3 = s$z}
    if (i == "4") {mat.4 = s$z}
    if (i == "5") {mat.5 = s$z}}
  
  png(paste(sample,"_",gene,"_3d_clusters.png", sep=""))
  par(mar=c(2.1, 6.1, 15.1, 2.1), xpd=TRUE)
  
  image3D(z = 0, phi=25, theta = 10, colvar = mat.1, smooth = TRUE, alpha = transparency*1, box = FALSE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F)
  image3D(z = 4, phi=25,  theta = 10, colvar = mat.2, smooth = TRUE, alpha = transparency*1, add = TRUE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F)
  image3D(z = 8, phi=25,  theta = 10, colvar = mat.3, smooth = TRUE, alpha = transparency*1, add = TRUE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F)
  image3D(z = 11, phi=25,  theta = 10, colvar = mat.4, smooth = TRUE, alpha = transparency*2, add = TRUE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F)
  image3D(z = 14, phi=25,  theta = 10, colvar = mat.5, smooth = TRUE, alpha = transparency*2.1, add = TRUE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F)
  dev.off()
}


### Plot cluster specific stuff
plot.gene.3d.cluster.4 = function(sample, cluster, m1, m2, m3, m4, s1, s2, s3, s4, x, y, transparency, min, max){
  
    for (i in c(1:4)) {
    if (i == "1") {
      spots = s1
      exp.values = m1}
    if (i == "2") {
      spots = s2
      exp.values = m2}
    if (i == "3") {
      spots = s3
      exp.values = m3}
    if (i == "4") {
      spots = s4
      exp.values = m4}

    clust.spots = rbind(s1,s2,s3,s4)
    rownames(clust.spots) = paste("X",rownames(clust.spots),sep="")
    genes.barcodes = clust.spots[rownames(clust.spots) %in% colnames(exp.values),]
    exp.values = exp.values[,colnames(exp.values) %in% rownames(clust.spots)]
    exp.values= exp.values[rowSums(exp.values) != 0,]
    cluster = ann.cluster
    cluster = cluster[colnames(exp.values),]
    genes.barcodes = cbind(genes.barcodes, cluster)
    
    x1 = as.numeric(genes.barcodes[,1])
    x1 = x1[!is.na(x1)]
    y1 = as.numeric(genes.barcodes[,2])
    y1 = y1[!is.na(y1)]
    z = as.numeric(genes.barcodes[,3])
    dis = rep(0.25 * as.numeric(i), nrow(genes.barcodes))
    w = as.numeric(dis)

    s =  interp(x1,y1, z, nx = x, ny = y)
    if (i == "1") {mat.1 = s$z}
    if (i == "2") {mat.2 = s$z}
    if (i == "3") {mat.3 = s$z}
    if (i == "4") {mat.4 = s$z}}

  
  png(paste(sample,"_3d_clusters.png", sep=""))
  par(mar=c(2.1, 6.1, 15.1, 2.1), xpd=TRUE)
  colorvar = colorRampPalette(c("#FF7F00","#FED8B1","gray61", "gray76", "gray85"))(5)
  image3D(z = 0, phi=25, theta = 10, colvar = mat.1, col=colorvar, smooth = TRUE, alpha = transparency*0.25, box = FALSE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F)
  image3D(z = 4, phi=25,  theta = 10, colvar = mat.2,  col=colorvar, smooth = TRUE, alpha = transparency*0.5, add = TRUE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F)
  image3D(z = 8, phi=25,  theta = 10, colvar = mat.3,  col=colorvar, smooth = TRUE, alpha = transparency*0.75, add = TRUE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F)
  image3D(z = 11, phi=25,  theta = 10, colvar = mat.4,  col=colorvar, smooth = TRUE, alpha = transparency*2, add = TRUE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F)
  dev.off()
  
}

plot.gene.3d.cluster.5 = function(sample, cluster, m1, m2, m3, m4, m5, s1, s2, s3, s4, s5, x, y, transparency, min, max){
  
  for (i in c(1:5)) {
    if (i == "1") {
      spots = s1
      exp.values = m1}
    if (i == "2") {
      spots = s2
      exp.values = m2}
    if (i == "3") {
      spots = s3
      exp.values = m3}
    if (i == "4") {
      spots = s4
      exp.values = m4}
    if (i == "5") {
      spots = s5
      exp.values = m5}
    
    clust.spots = rbind(s1,s2,s3,s4,s5)
    rownames(clust.spots) = paste("X",rownames(clust.spots),sep="")
    genes.barcodes = clust.spots[rownames(clust.spots) %in% colnames(exp.values),]
    exp.values = exp.values[,colnames(exp.values) %in% rownames(clust.spots)]
    exp.values= exp.values[rowSums(exp.values) != 0,]
    cluster = ann.cluster
    cluster = cluster[colnames(exp.values),]
    genes.barcodes = cbind(genes.barcodes, cluster)
    
    x1 = as.numeric(genes.barcodes[,1])
    x1 = x1[!is.na(x1)]
    y1 = as.numeric(genes.barcodes[,2])
    y1 = y1[!is.na(y1)]
    z = as.numeric(genes.barcodes[,3])
    dis = rep(0.25 * as.numeric(i), nrow(genes.barcodes))
    w = as.numeric(dis)
    
    s =  interp(x1,y1, z, nx = x, ny = y)
    if (i == "1") {mat.1 = s$z}
    if (i == "2") {mat.2 = s$z}
    if (i == "3") {mat.3 = s$z}
    if (i == "4") {mat.4 = s$z}
    if (i == "5") {mat.5 = s$z}}
  
  
  png(paste(sample,"_3d_clusters.png", sep=""))
  par(mar=c(2.1, 6.1, 15.1, 2.1), xpd=TRUE)
  colorvar = colorRampPalette(c("#FF7F00","#FED8B1","gray61", "gray76", "gray85"))(5)
  image3D(z = 0, phi=25, theta = 10, colvar = mat.1, col=colorvar, smooth = TRUE, alpha = transparency*0.25, box = FALSE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F)
  image3D(z = 4, phi=25,  theta = 10, colvar = mat.2,  col=colorvar, smooth = TRUE, alpha = transparency*0.5, add = TRUE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F)
  image3D(z = 8, phi=25,  theta = 10, colvar = mat.3,  col=colorvar, smooth = TRUE, alpha = transparency*0.75, add = TRUE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F)
  image3D(z = 11, phi=25,  theta = 10, colvar = mat.4,  col=colorvar, smooth = TRUE, alpha = transparency*2, add = TRUE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F)
  image3D(z = 14, phi=25,  theta = 10, colvar = mat.5,  col=colorvar, smooth = TRUE, alpha = transparency*2.1, add = TRUE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F)
  dev.off()
  
}

### Plot 2D clusters
#Funciton
plot.gene.2d.cluster.4 = function(sample, ann.cluster, m1, m2, m3, m4, s1, s2, s3, s4, x, y, transparency, min, max){
  mats_collection = matrix(ncol = 5, nrow = 1)
  colnames(mats_collection) = c("x", "y", "score", "sample", "section")
  
  for (i in c(1:4)) {
    print(i)
    if (i == "1") {
      spots = s1
      exp.values = m1}
    if (i == "2") {
      spots = s2
      exp.values = m2}
    if (i == "3") {
      spots = s3
      exp.values = m3}
    if (i == "4") {
      spots = s4
      exp.values = m4}
    
    clust.spots = rbind(s1,s2,s3,s4)
    # rownames(clust.spots) = paste("X",rownames(clust.spots),sep="")
    genes.barcodes = clust.spots[rownames(clust.spots) %in% colnames(exp.values),]
    exp.values = exp.values[,colnames(exp.values) %in% rownames(clust.spots)]
    exp.values = exp.values[rowSums(exp.values) != 0,]
    cluster = ann.cluster
    cluster = cluster[colnames(exp.values),]
    genes.barcodes = cbind(genes.barcodes, cluster)
    
    x1 = as.numeric(genes.barcodes[,1])
    x1 = x1[!is.na(x1)]
    y1 = as.numeric(genes.barcodes[,2])
    y1 = y1[!is.na(y1)]
    z = as.numeric(genes.barcodes[,3])
    dis = rep(0.25 * as.numeric(i), nrow(genes.barcodes))
    w = as.numeric(dis)
    
    #x= 40
    #y= 40
    #transparency=1
    # min=1
    # max=5
    s =  interp(x1, y1, z, nx = x, ny = y)
    mats_collection = rbind(mats_collection, make_output_tissue_matrix(s$z, norm_samples, i))
    if (i == "1") {mat.1 = s$z}
    if (i == "2") {mat.2 = s$z}
    if (i == "3") {mat.3 = s$z}
    if (i == "4") {mat.4 = s$z}}

  colorvar = colorRampPalette(c("gray85", "gray76","gray61","#FED8B1","#F49ACC"))(5)
  pdf(paste(sample,"_1_Contour_Cluster",".pdf", sep=""))
  image2D(z = mat.1, contour = T, smooth = TRUE, col=colorvar, alpha = transparency, box = FALSE, inttype = 1, NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F)
  dev.off()
  pdf(paste(sample,"_2_Contour_Cluster",".pdf", sep=""))
  image2D(z = mat.2, contour = T, smooth = TRUE, col=colorvar, alpha = transparency, box = FALSE, inttype = 1, NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F)
  dev.off()
  pdf(paste(sample,"_3_Contour_Cluster",".pdf", sep=""))
  image2D(z = mat.3, contour = T, smooth = TRUE, col=colorvar, alpha = transparency, box = FALSE, inttype = 1, NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F)
  dev.off()
  pdf(paste(sample,"_4_Contour_Cluster",".pdf", sep=""))
  image2D(z = mat.4, contour = T, smooth = TRUE, col=colorvar, alpha = transparency, box = FALSE, inttype = 1, NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F)
  dev.off()
  
  mats_collection = mats_collection[-1,]
  write.table(mats_collection, file = paste0("Spatial_clustering_table_", norm_samples, ".csv"), quote = F,sep = ",")
  
}
  
plot.gene.2d.4 = function(sample, gene,  mn1, mn2, mn3, mn4, s1, s2, s3, s4, x, y, transparency, min, max, con){
  
  mats_collection = matrix(ncol = 5, nrow = 1)
  colnames(mats_collection) = c("x", "y", "score", "sample", "section")
  
  for (i in c(1:4)) {
    if (i == "1") {
      spots = s1
      exp.values = mn1}
    if (i == "2") {
      spots = s2
      exp.values = mn2}
    if (i == "3") {
      spots = s3
      exp.values = mn3}
    if (i == "4") {
      spots = s4
      exp.values = mn4}

    genes.barcodes = exp.values[rownames(exp.values) == gene,]
    exp.values = exp.values[,colnames(exp.values) %in% rownames(spots)]
    exp.values= exp.values[rowSums(exp.values) != 0,]
    spots = spots[colnames(exp.values),]
    genes.barcodes.1 = cbind(spots, genes.barcodes)
    
    print(i)
    
    x1 = as.numeric(genes.barcodes.1[,1])
    x1 = x1[!is.na(x1)]
    y1 = as.numeric(genes.barcodes.1[,2])
    y1 = y1[!is.na(y1)]
    z = as.numeric(genes.barcodes.1[,3])
    dis = rep(0.25 * as.numeric(i), nrow(genes.barcodes.1))
    w = as.numeric(dis)
    
    s =  interp(x1,y1, z, nx = x, ny = y)
    mats_collection = rbind(mats_collection, make_output_tissue_matrix(s$z, norm_samples, i))
    if (i == "1") {mat.1 = s$z}
    if (i == "2") {mat.2 = s$z}
    if (i == "3") {mat.3 = s$z}
    if (i == "4") {mat.4 = s$z}}
  
  pdf(paste(sample,"_1_Contour_",gene,".pdf", sep=""))
  image2D(z = mat.1, contour = con, smooth = TRUE, alpha = transparency, box = FALSE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F)
  dev.off()
  pdf(paste(sample,"_2_Contour_",gene,".pdf", sep=""))
  image2D(z = mat.2, contour = con, smooth = TRUE, alpha = transparency, box = FALSE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F)
  dev.off()
  pdf(paste(sample,"_3_Contour_",gene,".pdf", sep=""))
  image2D(z = mat.3, contour = con, smooth = TRUE, alpha = transparency, box = FALSE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F)
  dev.off()
  pdf(paste(sample,"_4_Contour_",gene,".pdf", sep=""))
  image2D(z = mat.4, contour = con, smooth = TRUE, alpha = transparency, box = FALSE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F)
  dev.off()
  mats_collection = mats_collection[-1,]
  write.table(mats_collection, file = paste0("Gene_expression_tissue_table_", gene, "_", norm_samples, ".csv"), quote = F,sep = ",")
  
}

plot.gene.2d.5 = function(sample, gene,  m1, m2, m3, m4, m5, s1, s2, s3, s4, s5, x, y, transparency, min, max, con){
  mats_collection = matrix(ncol = 5, nrow = 1)
  colnames(mats_collection) = c("x", "y", "score", "sample", "section")
  
  for (i in c(1:5)) {
    if (i == "1") {
      spots = s1
      exp.values = m1}
    if (i == "2") {
      spots = s2
      exp.values = m2}
    if (i == "3") {
      spots = s3
      exp.values = m3}
    if (i == "4") {
      spots = s4
      exp.values = m4}
    if (i == "5") {
      spots = s5
      exp.values = m5}
    
    genes.barcodes = exp.values[rownames(exp.values) == gene,]
    exp.values = exp.values[,colnames(exp.values) %in% rownames(spots)]
    exp.values= exp.values[rowSums(exp.values) != 0,]
    spots = spots[colnames(exp.values),]
    genes.barcodes.1 = cbind(spots, genes.barcodes)
    
    print(i)
    
    x1 = as.numeric(genes.barcodes.1[,1])
    x1 = x1[!is.na(x1)]
    y1 = as.numeric(genes.barcodes.1[,2])
    y1 = y1[!is.na(y1)]
    z = as.numeric(genes.barcodes.1[,3])
    dis = rep(0.25 * as.numeric(i), nrow(genes.barcodes.1))
    w = as.numeric(dis)
    
    s =  interp(x1,y1, z, nx = x, ny = y)
    mats_collection = rbind(mats_collection, make_output_tissue_matrix(s$z, norm_samples, i))
    
    if (i == "1") {mat.1 = s$z}
    if (i == "2") {mat.2 = s$z}
    if (i == "3") {mat.3 = s$z}
    if (i == "4") {mat.4 = s$z}
    if (i == "5") {mat.5 = s$z}}
  
  pdf(paste(sample,"_1_Contour_",gene,".pdf", sep=""))
  image2D(z = mat.1, contour = con, smooth = TRUE, alpha = transparency, box = FALSE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F)
  dev.off()
  pdf(paste(sample,"_2_Contour_",gene,".pdf", sep=""))
  image2D(z = mat.2, contour = con, smooth = TRUE, alpha = transparency, box = FALSE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F)
  dev.off()
  pdf(paste(sample,"_3_Contour_",gene,".pdf", sep=""))
  image2D(z = mat.3, contour = con, smooth = TRUE, alpha = transparency, box = FALSE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F)
  dev.off()
  pdf(paste(sample,"_4_Contour_",gene,".pdf", sep=""))
  image2D(z = mat.4, contour = con, smooth = TRUE, alpha = transparency, box = FALSE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F)
  dev.off()
  pdf(paste(sample,"_5_Contour_",gene,".pdf", sep=""))
  image2D(z = mat.5, contour = con, smooth = TRUE, alpha = transparency, box = FALSE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F)
  dev.off()
  mats_collection = mats_collection[-1,]
  write.table(mats_collection, file = paste0("Gene_expression_tissue_table_", gene, "_", norm_samples, ".csv"), quote = F,sep = ",")
  
}

plot.gene.2d.cluster.7 = function(sample, cluster, m1, m2, m3, m4, m5, m6, m7, s1, s2, s3, s4, s5, s6, s7, x, y, transparency, min, max){
  
  mats_collection = matrix(ncol = 5, nrow = 1)
  colnames(mats_collection) = c("x", "y", "score", "sample", "section")
  
  for (i in c(1:7)) {
    if (i == "1") {
      spots = s1
      exp.values = m1}
    if (i == "2") {
      spots = s2
      exp.values = m2}
    if (i == "3") {
      spots = s3
      exp.values = m3}
    if (i == "4") {
      spots = s4
      exp.values = m4}
    if (i == "5") {
      spots = s5
      exp.values = m5}
    if (i == "6") {
      spots = s6
      exp.values = m6}
    if (i == "7") {
      spots = s7
      exp.values = m7}
    
    clust.spots = rbind(s1,s2,s3,s4,s5,s6,s7)
    # rownames(clust.spots) = paste("X",rownames(clust.spots),sep="")
    genes.barcodes = clust.spots[rownames(clust.spots) %in% colnames(exp.values),]
    exp.values = exp.values[,colnames(exp.values) %in% rownames(clust.spots)]
    exp.values = exp.values[rowSums(exp.values) != 0,]
    cluster = ann.cluster[colnames(exp.values),]
    genes.barcodes = cbind(genes.barcodes, cluster)
    
    x1 = as.numeric(genes.barcodes[,1])
    x1 = x1[!is.na(x1)]
    y1 = as.numeric(genes.barcodes[,2])
    y1 = y1[!is.na(y1)]
    z = as.numeric(genes.barcodes[,3])
    dis = rep(0.25 * as.numeric(i), nrow(genes.barcodes))
    w = as.numeric(dis)
    
    # x= 40
    # y= 40
    # transparency=1
    # min=1
    # max=5
    s =  interp(x1, y1, z, nx = x, ny = y)
    mats_collection = rbind(mats_collection, make_output_tissue_matrix(s$z, norm_samples, i))
    
    if (i == "1") {mat.1 = s$z}
    if (i == "2") {mat.2 = s$z}
    if (i == "3") {mat.3 = s$z}
    if (i == "4") {mat.4 = s$z}
    if (i == "5") {mat.5 = s$z}
    if (i == "6") {mat.6 = s$z}
    if (i == "7") {mat.7 = s$z}
    }
  
  colorvar = colorRampPalette(c("gray85", "gray76","gray61","#FED8B1","#F49ACC"))(5)
  pdf(paste(sample,"_1_Contour_Cluster",".pdf", sep=""))
  image2D(z = mat.1, contour = T, smooth = TRUE, col=colorvar, alpha = transparency, box = FALSE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F, colkey = list(plot = FALSE))
  dev.off()
  pdf(paste(sample,"_2_Contour_Cluster",".pdf", sep=""))
  image2D(z = mat.2, contour = T, smooth = TRUE, col=colorvar, alpha = transparency, box = FALSE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F, colkey = list(plot = FALSE))
  dev.off()
  pdf(paste(sample,"_3_Contour_Cluster",".pdf", sep=""))
  image2D(z = mat.3, contour = T, smooth = TRUE, col=colorvar, alpha = transparency, box = FALSE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F, colkey = list(plot = FALSE))
  dev.off()
  pdf(paste(sample,"_4_Contour_Cluster",".pdf", sep=""))
  image2D(z = mat.4, contour = T, smooth = TRUE, col=colorvar, alpha = transparency, box = FALSE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F, colkey = list(plot = FALSE))
  dev.off()
  pdf(paste(sample,"_5_Contour_Cluster",".pdf", sep=""))
  image2D(z = mat.5, contour = T, smooth = TRUE, col=colorvar, alpha = transparency, box = FALSE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F, colkey = list(plot = FALSE))
  dev.off()
  pdf(paste(sample,"_6_Contour_Cluster",".pdf", sep=""))
  image2D(z = mat.6, contour = T, smooth = TRUE, col=colorvar, alpha = transparency, box = FALSE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F, colkey = list(plot = FALSE))
  dev.off()
  pdf(paste(sample,"_7_Contour_Cluster",".pdf", sep=""))
  image2D(z = mat.7, contour = T, smooth = TRUE, col=colorvar, alpha = transparency, box = FALSE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F, colkey = list(plot = FALSE))
  dev.off()
  
  mats_collection = mats_collection[-1,]
  write.table(mats_collection, file = paste0("Spatial_clustering_table_", norm_samples, ".csv"), quote = F,sep = ",")
  
}

plot.gene.2d.cluster.5 = function(sample, ann.cluster, m1, m2, m3, m4, m5, s1, s2, s3, s4, s5, x, y, transparency, min, max){
  mats_collection = matrix(ncol = 5, nrow = 1)
  colnames(mats_collection) = c("x", "y", "score", "sample", "section")
  
  for (i in c(1:5)) {
    if (i == "1") {
      spots = s1
      exp.values = m1}
    if (i == "2") {
      spots = s2
      exp.values = m2}
    if (i == "3") {
      spots = s3
      exp.values = m3}
    if (i == "4") {
      spots = s4
      exp.values = m4}
    if (i == "5") {
      spots = s5
      exp.values = m5}

    cluster = ann.cluster
    clust.spots = rbind(s1,s2,s3,s4,s5)
    genes.barcodes = clust.spots[rownames(clust.spots) %in% colnames(exp.values),]
    exp.values = exp.values[,colnames(exp.values) %in% rownames(clust.spots)]
    exp.values = exp.values[rowSums(exp.values) != 0,]
    cluster = cluster[colnames(exp.values),]
    genes.barcodes = cbind(genes.barcodes, cluster)
    
    x1 = as.numeric(genes.barcodes[,1])
    x1 = x1[!is.na(x1)]
    y1 = as.numeric(genes.barcodes[,2])
    y1 = y1[!is.na(y1)]
    z = as.numeric(genes.barcodes[,3])
    dis = rep(0.25 * as.numeric(i), nrow(genes.barcodes))
    w = as.numeric(dis)
    
    # x=40
    # y=40
    s =  interp(x1, y1, z, nx = x, ny = y)
    mats_collection = rbind(mats_collection, make_output_tissue_matrix(s$z, norm_samples, i))
    
    if (i == "1") {mat.1 = s$z}
    if (i == "2") {mat.2 = s$z}
    if (i == "3") {mat.3 = s$z}
    if (i == "4") {mat.4 = s$z}
    if (i == "5") {mat.5 = s$z}
  }
  
  library(rgl)
  colorvar = colorRampPalette(c("gray85", "gray76","gray61","#FED8B1","#F49ACC"))(5)
  pdf(paste(sample,"_1_Contour_Cluster",".pdf", sep=""))
  image2D(z = mat.1, contour = T, smooth = TRUE, col=colorvar, alpha = transparency, box = FALSE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F)
  dev.off()
  pdf(paste(sample,"_2_Contour_Cluster",".pdf", sep=""))
  image2D(z = mat.2, contour = T, smooth = TRUE, col=colorvar, alpha = transparency, box = FALSE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F)
  dev.off()
  pdf(paste(sample,"_3_Contour_Cluster",".pdf", sep=""))
  image2D(z = mat.3, contour = T, smooth = TRUE, col=colorvar, alpha = transparency, box = FALSE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F)
  dev.off()
  pdf(paste(sample,"_4_Contour_Cluster",".pdf", sep=""))
  image2D(z = mat.4, contour = T, smooth = TRUE, col=colorvar, alpha = transparency, box = FALSE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F)
  dev.off()
  pdf(paste(sample,"_5_Contour_Cluster",".pdf", sep=""))
  image2D(z = mat.5, contour = T, smooth = TRUE, col=colorvar, alpha = transparency, box = FALSE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F)
  dev.off()
  
  mats_collection = mats_collection[-1,]
  write.table(mats_collection, file = paste0("Spatial_clustering_table_", norm_samples, ".csv"), quote = F,sep = ",")
  
}

plot.gene.2d.7 = function(sample, gene, m1, m2, m3, m4, m5, m6, m7, s1, s2, s3, s4, s5, s6, s7, x, y, transparency, min, max, con){ #, inf_all){
  mats_collection = matrix(ncol = 5, nrow = 1)
  colnames(mats_collection) = c("x", "y", "score", "sample", "section" ) #, "Inf.")
  #inf_all$names = row.names(inf_all)
  #print(inf_all)
  suppressWarnings({for (i in c(1:7)) {
    if (i == "1") {
      spots = s1
      exp.values = m1}
    if (i == "2") {
      spots = s2
      exp.values = m2}
    if (i == "3") {
      spots = s3
      exp.values = m3}
    if (i == "4") {
      spots = s4
      exp.values = m4}
    if (i == "5") {
      spots = s5
      exp.values = m5}
    if (i == "6") {
      spots = s6
      exp.values = m6}
    if (i == "7") {
      spots = s7
      exp.values = m7}
    
    # spots = s1
    # exp.values = mn1
    # gene = row.names(mf1)[1]
    # x = 40
    # y = 40
    
    genes.barcodes = exp.values[rownames(exp.values) == gene,]
    exp.values = exp.values[,colnames(exp.values) %in% rownames(spots)]
    exp.values= exp.values[rowSums(exp.values) != 0,]
    spots = spots[colnames(exp.values),]
    genes.barcodes.1 = cbind(spots, genes.barcodes)
    
    x1 = as.numeric(genes.barcodes.1[,1])
    x1 = x1[!is.na(x1)]
    y1 = as.numeric(genes.barcodes.1[,2])
    y1 = y1[!is.na(y1)]
    z = as.numeric(genes.barcodes.1[,3])

    s =  interp(x1, y1, z, nx = x, ny = y)
    # mat_tmp = make_output_tissue_matrix(s$z, norm_samples, i)
    # row.names(mat_tmp) = paste(str_replace_all(mat_tmp[,"section"], "Section", "X"), mat_tmp[,"x"], mat_tmp[,"y"], sep = "_")
    # mat_tmp.1 = merge(mat_tmp, inf_all, by = 0)
    # mat_tmp.1 = mat_tmp.1[,-1]
    # mat_tmp.1 = mat_tmp.1[,-7]
    # print(mat_tmp.1)
    # mats_collection = rbind(mats_collection, mat_tmp.1)
    mats_collection = rbind(mats_collection, make_output_tissue_matrix(s$z, norm_samples, i))
    
    # genes.barcodes.small = genes.barcodes.1[row.names(genes.barcodes.1) %in% row.names(inf_all[inf_all$Inf. == "Inf6",]),]
    # x1.small = as.numeric(genes.barcodes.small[,1])
    # x1.small = x1.small[!is.na(x1.small)]
    # y1.small = as.numeric(genes.barcodes.small[,2])
    # y1.small = y1.small[!is.na(y1.small)]
    # z.small = as.numeric(genes.barcodes.small[,3])
    # dis.small = rep(0.25 * as.numeric(i), nrow(genes.barcodes.small))
    # w = as.numeric(dis.small)
    # 
    # x.small = 100
    # y.small = 100
    # # transparency = 1
    # s.small =  interp(x1.small, y1.small, z.small, nx = x.small, ny = y.small)
    
    if (i == "1") {
      mat.1 = s$z}
      #mat.1.small = s.small$z}
    if (i == "2") {
      mat.2 = s$z}
      #mat.2.small = s.small$z}
    if (i == "3") {
      mat.3 = s$z}
      #mat.3.small = s.small$z}
    if (i == "4") {
      mat.4 = s$z}
      #mat.4.small = s.small$z}
    if (i == "5") {
      mat.5 = s$z}
      #mat.5.small = s.small$z}
    if (i == "6") {
      mat.6 = s$z}
      #mat.6.small = s.small$z}
    if (i == "7") {
      mat.7 = s$z}
      #mat.7.small = s.small$z}
  }
  })
  pdf(paste(sample,"_1_Contour_",gene,".pdf", sep=""))
  image2D(z = mat.1, contour = con, smooth = TRUE, alpha = transparency, box = FALSE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F)
  dev.off()
  pdf(paste(sample,"_2_Contour_",gene,".pdf", sep=""))
  image2D(z = mat.2, contour = con, smooth = TRUE, alpha = transparency, box = FALSE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F)
  dev.off()
  pdf(paste(sample,"_3_Contour_",gene,".pdf", sep=""))
  image2D(z = mat.3, contour = con, smooth = TRUE, alpha = transparency, box = FALSE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F)
  dev.off()
  pdf(paste(sample,"_4_Contour_",gene,".pdf", sep=""))
  image2D(z = mat.4, contour = con, smooth = TRUE, alpha = transparency, box = FALSE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F)
  dev.off()
  pdf(paste(sample,"_5_Contour_",gene,".pdf", sep=""))
  image2D(z = mat.5, contour = con, smooth = TRUE, alpha = transparency, box = FALSE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F)
  dev.off()
  pdf(paste(sample,"_6_Contour_",gene,".pdf", sep=""))
  image2D(z = mat.6, contour = con, smooth = TRUE, alpha = transparency, box = FALSE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F)
  dev.off()
  pdf(paste(sample,"_7_Contour_",gene,".pdf", sep=""))
  image2D(z = mat.7, contour = con, smooth = TRUE, alpha = transparency, box = FALSE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F)
  dev.off()
  mats_collection = mats_collection[-1,]
  write.table(mats_collection, file = paste0("Gene_expression_tissue_table_", gene, "_", norm_samples, ".csv"), quote = F,sep = ",")
  
  
  # pdf(paste(sample,"_1_Contour_Inf6_",gene,".pdf", sep=""))
  # image2D(z = mat.1.small, contour = T, smooth = TRUE, alpha = transparency, box = FALSE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x.small), y = c(1:y.small), scale = F, colkey = list(plot = FALSE))
  # dev.off()
  # pdf(paste(sample,"_2_Contour_Inf6_",gene,".pdf", sep=""))
  # image2D(z = mat.2.small, contour = T, smooth = TRUE, alpha = transparency, box = FALSE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x.small), y = c(1:y.small), scale = F, colkey = list(plot = FALSE))
  # dev.off()
  # pdf(paste(sample,"_3_Contour_Inf6_",gene,".pdf", sep=""))
  # image2D(z = mat.3.small, contour = T, smooth = TRUE, alpha = transparency, box = FALSE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x.small), y = c(1:y.small), scale = F, colkey = list(plot = FALSE))
  # dev.off()
  # pdf(paste(sample,"_4_Contour_Inf6_",gene,".pdf", sep=""))
  # image2D(z = mat.4.small, contour = T, smooth = TRUE, alpha = transparency, box = FALSE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x.small), y = c(1:y.small), scale = F, colkey = list(plot = FALSE))
  # dev.off()
  # pdf(paste(sample,"_5_Contour_Inf6_",gene,".pdf", sep=""))
  # image2D(z = mat.5.small, contour = T, smooth = TRUE, alpha = transparency, box = FALSE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x.small), y = c(1:y.small), scale = F, colkey = list(plot = FALSE))
  # dev.off()
  # pdf(paste(sample,"_6_Contour_Inf6_",gene,".pdf", sep=""))
  # image2D(z = mat.6.small, contour = T, smooth = TRUE, alpha = transparency, box = FALSE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x.small), y = c(1:y.small), scale = F, colkey = list(plot = FALSE))
  # dev.off()
  # pdf(paste(sample,"_7_Contour_Inf6_",gene,".pdf", sep=""))
  # image2D(z = mat.7.small, contour = T, smooth = TRUE,  alpha = transparency, box = FALSE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x.small), y = c(1:y.small), scale = F, colkey = list(plot = FALSE))
  # dev.off()
  
  # plot inf6 for figshare
  #print(head(mats_collection))
  # row.names(mats_collection) = paste(str_replace_all(mats_collection[,"section"], "Section", "X"), mats_collection[,"x"], mats_collection[,"y"], sep = "_")
  # mats_collection = mats_collection[row.names(mats_collection) %in% row.names(inf_all[inf_all$Inf. == "Inf6",]),]
  # write.table(mats_collection, file = paste0("Gene_expression_tissue_table_Inf6_", gene, "_", norm_samples, ".csv"), quote = F,sep = ",")
  
  
}

plot.gene.2d.inf.7 = function(sample, cluster, m1, m2, m3, m4, m5, m6, m7, s1, s2, s3, s4, s5, s6, s7, x, y, transparency, min, max){ #, inf_all)
  
  #mats_collection = matrix(ncol = 5, nrow = 1)
  #colnames(mats_collection) = c("x", "y", "score", "sample", "section")
  #inf_all$names = row.names(inf_all)
  
  for (i in c(1:7)) {
    if (i == "1") {
      spots = s1
      exp.values = m1}
    if (i == "2") {
      spots = s2
      exp.values = m2}
    if (i == "3") {
      spots = s3
      exp.values = m3}
    if (i == "4") {
      spots = s4
      exp.values = m4}
    if (i == "5") {
      spots = s5
      exp.values = m5}
    if (i == "6") {
      spots = s6
      exp.values = m6}
    if (i == "7") {
      spots = s7
      exp.values = m7}

    
    clust.spots = rbind(s1,s2,s3,s4,s5,s6,s7)
    genes.barcodes = clust.spots[rownames(clust.spots) %in% colnames(exp.values),]
    exp.values = exp.values[,colnames(exp.values) %in% rownames(clust.spots)]
    exp.values = exp.values[rowSums(exp.values) != 0,]
    cluster = ann.col.inf
    cluster = cluster[colnames(exp.values),]
    genes.barcodes = cbind(genes.barcodes, cluster)
    
    x1 = as.numeric(genes.barcodes[,1])
    x1 = x1[!is.na(x1)]
    y1 = as.numeric(genes.barcodes[,2])
    y1 = y1[!is.na(y1)]
    z = as.numeric(genes.barcodes[,3])
    dis = rep(0.25 * as.numeric(i), nrow(genes.barcodes))
    w = as.numeric(dis)
    
    # x = 40
    # y = 40
    # transparency = 1
    s =  interp(x1, y1, z, nx = x, ny = y)
    #mats_collection = rbind(mats_collection, make_output_tissue_matrix(s$z, norm_samples, i))
    
    
    #genes.barcodes.small = genes.barcodes[row.names(genes.barcodes) %in% row.names(inf_all[inf_all$Inf. == "Inf6",]),]
    #x1.small = as.numeric(genes.barcodes.small[,1])
    #x1.small = x1.small[!is.na(x1.small)]
    #y1.small = as.numeric(genes.barcodes.small[,2])
    #y1.small = y1.small[!is.na(y1.small)]
    #z.small = as.numeric(genes.barcodes.small[,3])
    #dis.small = rep(0.25 * as.numeric(i), nrow(genes.barcodes.small))
    #w = as.numeric(dis.small)

    #x.small = 100
    #y.small = 100
    # transparency = 1
    #s.small =  interp(x1.small, y1.small, z.small, nx = x.small, ny = y.small)

    if (i == "1") {
      mat.1 = s$z}
      #mat.1.small = s.small$z}
    if (i == "2") {
      mat.2 = s$z}
      #mat.2.small = s.small$z}
    if (i == "3") {
      mat.3 = s$z}
      #mat.3.small = s.small$z}
    if (i == "4") {
      mat.4 = s$z}
      #mat.4.small = s.small$z}
    if (i == "5") {
      mat.5 = s$z}
      #mat.5.small = s.small$z}
    if (i == "6") {
      mat.6 = s$z}
      #mat.6.small = s.small$z}
    if (i == "7") {
      mat.7 = s$z}
      #mat.7.small = s.small$z}
  }

  colorvar = colorRampPalette(c("#8a3795", "#4ebceb","#f15f48","gray80"))(4)
  pdf(paste(sample,"_1_Contour_Inf",".pdf", sep=""))
  image2D(z = mat.1, contour = T, smooth = TRUE, col=colorvar, alpha = transparency, box = FALSE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F, colkey = list(plot = FALSE))
  dev.off()
  pdf(paste(sample,"_2_Contour_Inf",".pdf", sep=""))
  image2D(z = mat.2, contour = T, smooth = TRUE, col=colorvar, alpha = transparency, box = FALSE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F, colkey = list(plot = FALSE))
  dev.off()
  pdf(paste(sample,"_3_Contour_Inf",".pdf", sep=""))
  image2D(z = mat.3, contour = T, smooth = TRUE, col=colorvar, alpha = transparency, box = FALSE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F, colkey = list(plot = FALSE))
  dev.off()
  pdf(paste(sample,"_4_Contour_Inf",".pdf", sep=""))
  image2D(z = mat.4, contour = T, smooth = TRUE, col=colorvar, alpha = transparency, box = FALSE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F, colkey = list(plot = FALSE))
  dev.off()
  pdf(paste(sample,"_5_Contour_Inf",".pdf", sep=""))
  image2D(z = mat.5, contour = T, smooth = TRUE, col=colorvar, alpha = transparency, box = FALSE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F, colkey = list(plot = FALSE))
  dev.off()
  pdf(paste(sample,"_6_Contour_Inf",".pdf", sep=""))
  image2D(z = mat.6, contour = T, smooth = TRUE, col=colorvar, alpha = transparency, box = FALSE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F, colkey = list(plot = FALSE))
  dev.off()
  pdf(paste(sample,"_7_Contour_Inf",".pdf", sep=""))
  image2D(z = mat.7, contour = T, smooth = TRUE, col=colorvar, alpha = transparency, box = FALSE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F, colkey = list(plot = FALSE))
  dev.off()
  # pdf(paste(sample,"_1_Contour_Inf6",".pdf", sep=""))
  # image2D(z = mat.1.small, contour = T, smooth = TRUE, col=colorvar, alpha = transparency, box = FALSE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x.small), y = c(1:y.small), scale = F, colkey = list(plot = FALSE))
  # dev.off()
  # pdf(paste(sample,"_2_Contour_Inf6",".pdf", sep=""))
  # image2D(z = mat.2.small, contour = T, smooth = TRUE, col=colorvar, alpha = transparency, box = FALSE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x.small), y = c(1:y.small), scale = F, colkey = list(plot = FALSE))
  # dev.off()
  # pdf(paste(sample,"_3_Contour_Inf6",".pdf", sep=""))
  # image2D(z = mat.3.small, contour = T, smooth = TRUE, col=colorvar, alpha = transparency, box = FALSE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x.small), y = c(1:y.small), scale = F, colkey = list(plot = FALSE))
  # dev.off()
  # pdf(paste(sample,"_4_Contour_Inf6",".pdf", sep=""))
  # image2D(z = mat.4.small, contour = T, smooth = TRUE, col=colorvar, alpha = transparency, box = FALSE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x.small), y = c(1:y.small), scale = F, colkey = list(plot = FALSE))
  # dev.off()
  # pdf(paste(sample,"_5_Contour_Inf6",".pdf", sep=""))
  # image2D(z = mat.5.small, contour = T, smooth = TRUE, col=colorvar, alpha = transparency, box = FALSE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x.small), y = c(1:y.small), scale = F, colkey = list(plot = FALSE))
  # dev.off()
  # pdf(paste(sample,"_6_Contour_Inf6",".pdf", sep=""))
  # image2D(z = mat.6.small, contour = T, smooth = TRUE, col=colorvar, alpha = transparency, box = FALSE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x.small), y = c(1:y.small), scale = F, colkey = list(plot = FALSE))
  # dev.off()
  # pdf(paste(sample,"_7_Contour_Inf6",".pdf", sep=""))
  # image2D(z = mat.7.small, contour = T, smooth = TRUE, col=colorvar, alpha = transparency, box = FALSE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x.small), y = c(1:y.small), scale = F, colkey = list(plot = FALSE))
  # dev.off()
  
  #mats_collection = mats_collection[-1,]
  #write.table(mats_collection, file = paste0("Inf_clustering_table_", norm_samples, ".csv"), quote = F,sep = ",")
  
  # plot inf6 for figshare
  #print(head(mats_collection))
  #row.names(mats_collection) = paste(str_replace_all(mats_collection[,"section"], "Section", "X"), mats_collection[,"x"], mats_collection[,"y"], sep = "_")
  #mats_collection = mats_collection[row.names(mats_collection) %in% row.names(inf_all[inf_all$Inf. == "Inf6",]),]
  #write.table(mats_collection, file = paste0("Inf_clustering_table_Inf6_", norm_samples, ".csv"), quote = F,sep = ",")
  
}

make_output_tissue_matrix = function(mat.1, norm_samples, i) {

  xs = ""
  ys = ""
  scores = ""
  for (xxx in 1:nrow(mat.1)){
    for (yyy in 1:ncol(mat.1)){
      xs = c(xs,xxx)
      ys = c(ys,yyy)
      scores = c(scores, mat.1[xxx,yyy])
    }
  }

  mat_output = cbind(as.numeric(xs), as.numeric(ys),as.numeric(scores),rep(norm_samples, length(xs)),rep(paste0("Section",i), length(xs)))
  mat_output = mat_output[-1,]
  colnames(mat_output) = c("x", "y", "score", "sample", "section")

  return(mat_output)
}


plot.gene.2d.inf.4 = function(sample, ann.col.inf, m1, m2, m3, m4, s1, s2, s3, s4, x, y, transparency, min, max){
    
    mats_collection = matrix(ncol = 5, nrow = 1)
    colnames(mats_collection) = c("x", "y", "score", "sample", "section")
    for (i in c(1:4)) {
      print(i)
        if (i == "1") {
            spots = s1
            exp.values = m1}
        if (i == "2") {
            spots = s2
            exp.values = m2}
        if (i == "3") {
            spots = s3
            exp.values = m3}
        if (i == "4") {
            spots = s4
            exp.values = m4}
        
        clust.spots = rbind(s1,s2,s3,s4)
        genes.barcodes = clust.spots[rownames(clust.spots) %in% colnames(exp.values),]
        exp.values = exp.values[,colnames(exp.values) %in% rownames(clust.spots)]
        exp.values = exp.values[rowSums(exp.values) != 0,]
        cluster = ann.col.inf
        cluster = cluster[colnames(exp.values),]
        genes.barcodes = cbind(genes.barcodes, cluster)
        
        x1 = as.numeric(genes.barcodes[,1])
        x1 = x1[!is.na(x1)]
        y1 = as.numeric(genes.barcodes[,2])
        y1 = y1[!is.na(y1)]
        z = as.numeric(genes.barcodes[,3])
        dis = rep(0.25 * as.numeric(i), nrow(genes.barcodes))
        w = as.numeric(dis)
        
        #x = 40
        #y = 40
        #transparency = 1
        s =  interp(x1, y1, z, nx = x, ny = y)
        mats_collection = rbind(mats_collection, make_output_tissue_matrix(s$z, norm_samples, i))

        if (i == "1") {mat.1 = s$z}
        if (i == "2") {mat.2 = s$z}
        if (i == "3") {mat.3 = s$z}
        if (i == "4") {mat.4 = s$z}
    }
  

    if (max(cluster) == 3) {colorvar = colorRampPalette(c("#4ebceb", "#8a3795","gray80"))(3)}
    if (max(cluster) == 4) {colorvar = colorRampPalette(c("#4ebceb", "#8a3795", "#f15f48","gray80"))(4)}

    pdf(paste(sample,"_1_Contour_Inf",".pdf", sep=""))
    image2D(z = mat.1, contour = T, smooth = TRUE, col=colorvar, alpha = transparency, box = FALSE, inttype = 1, NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F)
    dev.off()
    pdf(paste(sample,"_2_Contour_Inf",".pdf", sep=""))
    image2D(z = mat.2, contour = T, smooth = TRUE, col=colorvar, alpha = transparency, box = FALSE, inttype = 1, NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F)
    dev.off()
    pdf(paste(sample,"_3_Contour_Inf",".pdf", sep=""))
    image2D(z = mat.3, contour = T, smooth = TRUE, col=colorvar, alpha = transparency, box = FALSE, inttype = 1, NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F)
    dev.off()
    pdf(paste(sample,"_4_Contour_Inf",".pdf", sep=""))
    image2D(z = mat.4, contour = T, smooth = TRUE, col=colorvar, alpha = transparency, box = FALSE, inttype = 1, NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F)
    dev.off()
    
    mats_collection = mats_collection[-1,]
    write.table(mats_collection, file = paste0("Inf_clustering_table_", norm_samples, ".csv"), quote = F,sep = ",")
    
}

plot.gene.2d.inf.5 = function(sample, cluster, m1, m2, m3, m4, m5, s1, s2, s3, s4, s5, x, y, transparency, min, max){
  mats_collection = matrix(ncol = 5, nrow = 1)
  colnames(mats_collection) = c("x", "y", "score", "sample", "section")
  
  for (i in c(1:5)) {
    if (i == "1") {
      spots = s1
      exp.values = m1}
    if (i == "2") {
      spots = s2
      exp.values = m2}
    if (i == "3") {
      spots = s3
      exp.values = m3}
    if (i == "4") {
      spots = s4
      exp.values = m4}
    if (i == "5") {
      spots = s5
      exp.values = m5}

    clust.spots = rbind(s1,s2,s3,s4,s5)
    genes.barcodes = clust.spots[rownames(clust.spots) %in% colnames(exp.values),]
    exp.values = exp.values[,colnames(exp.values) %in% rownames(clust.spots)]
    exp.values = exp.values[rowSums(exp.values) != 0,]
    cluster = ann.col.inf
    cluster = cluster[colnames(exp.values),]
    genes.barcodes = cbind(genes.barcodes, cluster)
    
    x1 = as.numeric(genes.barcodes[,1])
    x1 = x1[!is.na(x1)]
    y1 = as.numeric(genes.barcodes[,2])
    y1 = y1[!is.na(y1)]
    z = as.numeric(genes.barcodes[,3])
    dis = rep(0.25 * as.numeric(i), nrow(genes.barcodes))
    w = as.numeric(dis)
    
    # x = 40
    # y = 40
    # transparency = 1

    s =  interp(x1, y1, z, nx = x, ny = y)
    mats_collection = rbind(mats_collection, make_output_tissue_matrix(s$z, norm_samples, i))
    
    if (i == "1") {mat.1 = s$z}
    if (i == "2") {mat.2 = s$z}
    if (i == "3") {mat.3 = s$z}
    if (i == "4") {mat.4 = s$z}
    if (i == "5") {mat.5 = s$z}
  }
  
  library(rgl)
  
  colorvar = colorRampPalette(c("#8a3795", "#4ebceb","#f15f48","gray80"))(4)
  pdf(paste(sample,"_1_Contour_Inf",".pdf", sep=""))
  image2D(z = mat.1, contour = T, smooth = TRUE, col=colorvar, alpha = transparency, box = FALSE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F)
  dev.off()
  pdf(paste(sample,"_2_Contour_Inf",".pdf", sep=""))
  image2D(z = mat.2, contour = T, smooth = TRUE, col=colorvar, alpha = transparency, box = FALSE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F)
  dev.off()
  pdf(paste(sample,"_3_Contour_Inf",".pdf", sep=""))
  image2D(z = mat.3, contour = T, smooth = TRUE, col=colorvar, alpha = transparency, box = FALSE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F)
  dev.off()
  pdf(paste(sample,"_4_Contour_Inf",".pdf", sep=""))
  image2D(z = mat.4, contour = T, smooth = TRUE, col=colorvar, alpha = transparency, box = FALSE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F)
  dev.off()
  pdf(paste(sample,"_5_Contour_Inf",".pdf", sep=""))
  image2D(z = mat.5, contour = T, smooth = TRUE, col=colorvar, alpha = transparency, box = FALSE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F)
  dev.off()
  mats_collection = mats_collection[-1,]
  write.table(mats_collection, file = paste0("Inf_clustering_table_", norm_samples, ".csv"), quote = F,sep = ",")
  

}

plot.gene.2d.inf.3 = function(sample, cluster, m1, m2, m3, s1, s2, s3, x, y, transparency, min, max){
  mats_collection = matrix(ncol = 5, nrow = 1)
  colnames(mats_collection) = c("x", "y", "score", "sample", "section")
  
  for (i in c(1:3)) {
    if (i == "1") {
      spots = s1
      exp.values = m1}
    if (i == "2") {
      spots = s2
      exp.values = m2}
    if (i == "3") {
      spots = s3
      exp.values = m3}
    
    clust.spots = rbind(s1,s2,s3)
    genes.barcodes = clust.spots[rownames(clust.spots) %in% colnames(exp.values),]
    exp.values = exp.values[,colnames(exp.values) %in% rownames(clust.spots)]
    exp.values = exp.values[rowSums(exp.values) != 0,]
    cluster = ann.col.inf
    cluster = cluster[colnames(exp.values),]
    genes.barcodes = cbind(genes.barcodes, cluster)
    
    x1 = as.numeric(genes.barcodes[,1])
    x1 = x1[!is.na(x1)]
    y1 = as.numeric(genes.barcodes[,2])
    y1 = y1[!is.na(y1)]
    z = as.numeric(genes.barcodes[,3])
    dis = rep(0.25 * as.numeric(i), nrow(genes.barcodes))
    w = as.numeric(dis)
    
    # x = 40
    # y = 40
    # transparency = 1
    s =  interp(x1, y1, z, nx = x, ny = y)
    mats_collection = rbind(mats_collection, make_output_tissue_matrix(s$z, norm_samples, i))
    
    if (i == "1") {mat.1 = s$z}
    if (i == "2") {mat.2 = s$z}
    if (i == "3") {mat.3 = s$z}

  }
  
  colorvar = colorRampPalette(c("#8a3795", "#4ebceb","#f15f48","gray80"))(4)
  pdf(paste(sample,"_1_Contour_Inf",".pdf", sep=""))
  image2D(z = mat.1, contour = T, smooth = TRUE, col=colorvar, alpha = transparency, box = FALSE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F)
  dev.off()
  pdf(paste(sample,"_2_Contour_Inf",".pdf", sep=""))
  image2D(z = mat.2, contour = T, smooth = TRUE, col=colorvar, alpha = transparency, box = FALSE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F)
  dev.off()
  pdf(paste(sample,"_3_Contour_Inf",".pdf", sep=""))
  image2D(z = mat.3, contour = T, smooth = TRUE, col=colorvar, alpha = transparency, box = FALSE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F)
  dev.off()
  mats_collection = mats_collection[-1,]
  write.table(mats_collection, file = paste0("Inf_clustering_table_", norm_samples, ".csv"), quote = F,sep = ",")
  

  
}


plot.gene.2d.cluster.3 = function(sample, ann.cluster, m1, m2, m3, s1, s2, s3, x, y, transparency, min, max){
  mats_collection = matrix(ncol = 5, nrow = 1)
  colnames(mats_collection) = c("x", "y", "score", "sample", "section")
  
  for (i in c(1:3)) {
    if (i == "1") {
      spots = s1
      exp.values = m1}
    if (i == "2") {
      spots = s2
      exp.values = m2}
    if (i == "3") {
      spots = s3
      exp.values = m3}
    
    cluster = ann.cluster
    clust.spots = rbind(s1,s2,s3)
    genes.barcodes = clust.spots[rownames(clust.spots) %in% colnames(exp.values),]
    exp.values = exp.values[,colnames(exp.values) %in% rownames(clust.spots)]
    exp.values = exp.values[rowSums(exp.values) != 0,]
    cluster = cluster[colnames(exp.values),]
    genes.barcodes = cbind(genes.barcodes, cluster)
    
    x1 = as.numeric(genes.barcodes[,1])
    x1 = x1[!is.na(x1)]
    y1 = as.numeric(genes.barcodes[,2])
    y1 = y1[!is.na(y1)]
    z = as.numeric(genes.barcodes[,3])
    dis = rep(0.25 * as.numeric(i), nrow(genes.barcodes))
    w = as.numeric(dis)
    

    # x=40
    # y=40
    s =  interp(x1, y1, z, nx = x, ny = y)
    mats_collection = rbind(mats_collection, make_output_tissue_matrix(s$z, norm_samples, i))
    
    if (i == "1") {mat.1 = s$z}
    if (i == "2") {mat.2 = s$z}
    if (i == "3") {mat.3 = s$z}

  }
  
  library(rgl)
  colorvar = colorRampPalette(c("gray85", "gray76","gray61","#FED8B1","#F49ACC"))(5)
  pdf(paste(sample,"_1_Contour_Cluster",".pdf", sep=""))
  image2D(z = mat.1, contour = T, smooth = TRUE, col=colorvar, alpha = transparency, box = FALSE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F)
  dev.off()
  pdf(paste(sample,"_2_Contour_Cluster",".pdf", sep=""))
  image2D(z = mat.2, contour = T, smooth = TRUE, col=colorvar, alpha = transparency, box = FALSE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F)
  dev.off()
  pdf(paste(sample,"_3_Contour_Cluster",".pdf", sep=""))
  image2D(z = mat.3, contour = T, smooth = TRUE, col=colorvar, alpha = transparency, box = FALSE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F)
  dev.off()
  
  mats_collection = mats_collection[-1,]
  write.table(mats_collection, file = paste0("Spatial_clustering_table_", norm_samples, ".csv"), quote = F,sep = ",")
  

}

plot.gene.2d.3 = function(sample, gene,  m1, m2, m3,  s1, s2, s3,  x, y, transparency, min, max, con){
  mats_collection = matrix(ncol = 5, nrow = 1)
  colnames(mats_collection) = c("x", "y", "score", "sample", "section")
  
  for (i in c(1:3)) {
    if (i == "1") {
      spots = s1
      exp.values = m1}
    if (i == "2") {
      spots = s2
      exp.values = m2}
    if (i == "3") {
      spots = s3
      exp.values = m3}
    
    genes.barcodes = exp.values[rownames(exp.values) == gene,]
    exp.values = exp.values[,colnames(exp.values) %in% rownames(spots)]
    exp.values= exp.values[rowSums(exp.values) != 0,]
    spots = spots[colnames(exp.values),]
    genes.barcodes.1 = cbind(spots, genes.barcodes)
    
    print(i)
    
    x1 = as.numeric(genes.barcodes.1[,1])
    x1 = x1[!is.na(x1)]
    y1 = as.numeric(genes.barcodes.1[,2])
    y1 = y1[!is.na(y1)]
    z = as.numeric(genes.barcodes.1[,3])
    dis = rep(0.25 * as.numeric(i), nrow(genes.barcodes.1))
    w = as.numeric(dis)
    
    s =  interp(x1,y1, z, nx = x, ny = y)
    mats_collection = rbind(mats_collection, make_output_tissue_matrix(s$z, norm_samples, i))
    if (i == "1") {mat.1 = s$z}
    if (i == "2") {mat.2 = s$z}
    if (i == "3") {mat.3 = s$z}
  }
  
  pdf(paste(sample,"_1_Contour_",gene,".pdf", sep=""))
  image2D(z = mat.1, contour = con, smooth = TRUE, alpha = transparency, box = FALSE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F)
  dev.off()
  pdf(paste(sample,"_2_Contour_",gene,".pdf", sep=""))
  image2D(z = mat.2, contour = con, smooth = TRUE, alpha = transparency, box = FALSE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F)
  dev.off()
  pdf(paste(sample,"_3_Contour_",gene,".pdf", sep=""))
  image2D(z = mat.3, contour = con, smooth = TRUE, alpha = transparency, box = FALSE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F)
  dev.off()
  
  mats_collection = mats_collection[-1,]
  write.table(mats_collection, file = paste0("Gene_expression_tissue_table_", gene, "_", norm_samples, ".csv"), quote = F,sep = ",")
  

}

add.flag <- function(pheatmap,
                     kept.labels,
                     repel.degree) {
  
  # repel.degree = number within [0, 1], which controls how much 
  #                space to allocate for repelling labels.
  ## repel.degree = 0: spread out labels over existing range of kept labels
  ## repel.degree = 1: spread out labels over the full y-axis
  
  heatmap <- pheatmap$gtable
  
  new.label <- heatmap$grobs[[which(heatmap$layout$name == "row_names")]] 
  
  # keep only labels in kept.labels, replace the rest with ""
  new.label$label <- ifelse(new.label$label %in% kept.labels, 
                            new.label$label, "")
  
  # calculate evenly spaced out y-axis positions
  repelled.y <- function(d, d.select, k = repel.degree){
    # d = vector of distances for labels
    # d.select = vector of T/F for which labels are significant
    
    # recursive function to get current label positions
    # (note the unit is "npc" for all components of each distance)
    strip.npc <- function(dd){
      if(!"unit.arithmetic" %in% class(dd)) {
        return(as.numeric(dd))
      }
      
      d1 <- strip.npc(dd$arg1)
      d2 <- strip.npc(dd$arg2)
      fn <- dd$fname
      return(lazyeval::lazy_eval(paste(d1, fn, d2)))
    }
    
    full.range <- sapply(seq_along(d), function(i) strip.npc(d[i]))
    selected.range <- sapply(seq_along(d[d.select]), function(i) strip.npc(d[d.select][i]))
    
    return(unit(seq(from = max(selected.range) + k*(max(full.range) - max(selected.range)),
                    to = min(selected.range) - k*(min(selected.range) - min(full.range)), 
                    length.out = sum(d.select)), 
                "npc"))
  }
  new.y.positions <- repelled.y(new.label$y,
                                d.select = new.label$label != "")
  new.flag <- segmentsGrob(x0 = new.label$x,
                           x1 = new.label$x + unit(0.15, "npc"),
                           y0 = new.label$y[new.label$label != ""],
                           y1 = new.y.positions)
  
  # shift position for selected labels
  new.label$x <- new.label$x + unit(0.2, "npc")
  new.label$y[new.label$label != ""] <- new.y.positions
  
  # add flag to heatmap
  heatmap <- gtable::gtable_add_grob(x = heatmap,
                                     grobs = new.flag,
                                     t = 4, 
                                     l = 4
  )
  
  # replace label positions in heatmap
  heatmap$grobs[[which(heatmap$layout$name == "row_names")]] <- new.label
  
  # plot result
  grid.newpage()
  grid.draw(heatmap)
  
  # return a copy of the heatmap invisibly
  invisible(heatmap)
}

stat_boxplot_custom <- function(mapping = NULL, data = NULL,
                                geom = "boxplot", position = "dodge",
                                ...,
                                qs = c(.05, .25, 0.5, 0.75, 0.95),
                                na.rm = FALSE,
                                show.legend = NA,
                                inherit.aes = TRUE) {
  layer(
    data = data,
    mapping = mapping,
    stat = StatBoxplotCustom,
    geom = geom,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      na.rm = na.rm,
      qs = qs,
      ...
    )
  )
}
StatBoxplotCustom <- ggproto("StatBoxplotCustom", Stat,
                             required_aes = c("x", "y"),
                             non_missing_aes = "weight",
                             
                             setup_params = function(data, params) {
                               params$width <- ggplot2:::"%||%"(
                                 params$width, (resolution(data$x) * 0.75)
                               )
                               
                               if (is.double(data$x) && !ggplot2:::has_groups(data) && any(data$x != data$x[1L])) {
                                 warning(
                                   "Continuous x aesthetic -- did you forget aes(group=...)?",
                                   call. = FALSE
                                 )
                               }
                               
                               params
                             },
                             
                             compute_group = function(data, scales, width = NULL, na.rm = FALSE, qs = c(.05, .25, 0.5, 0.75, 0.95)) {
                               
                               if (!is.null(data$weight)) {
                                 mod <- quantreg::rq(y ~ 1, weights = weight, data = data, tau = qs)
                                 stats <- as.numeric(stats::coef(mod))
                               } else {
                                 stats <- as.numeric(stats::quantile(data$y, qs))
                               }
                               names(stats) <- c("ymin", "lower", "middle", "upper", "ymax")
                               iqr <- diff(stats[c(2, 4)])
                               
                               outliers <- (data$y < stats[1]) | (data$y > stats[5])
                               
                               if (length(unique(data$x)) > 1)
                                 width <- diff(range(data$x)) * 0.9
                               
                               df <- as.data.frame(as.list(stats))
                               df$outliers <- list(data$y[outliers])
                               
                               if (is.null(data$weight)) {
                                 n <- sum(!is.na(data$y))
                               } else {
                                 # Sum up weights for non-NA positions of y and weight
                                 n <- sum(data$weight[!is.na(data$y) & !is.na(data$weight)])
                               }
                               
                               df$notchupper <- df$middle + 1.58 * iqr / sqrt(n)
                               df$notchlower <- df$middle - 1.58 * iqr / sqrt(n)
                               
                               df$x <- if (is.factor(data$x)) data$x[1] else mean(range(data$x))
                               df$width <- width
                               df$relvarwidth <- sqrt(n)
                               df
                             }
)


avg_genes_box_old = function(all_bar, gen_names, cluster_names){
  
  #### Make boxplot of avg expression per interesting marker genes
  col.spatial.clusters = cluster_names
  col.scale = unique(col.spatial.clusters[,1])
  
  ### arrange matrix per cluster 
  #all_bar = RA.norm
  all_bar = all_bar[row.names(all_bar) %in% gen_names,]
  all_bar = all_bar[,match(colnames(all_bar), rownames(col.spatial.clusters))]
  all_bar = all_bar[match(rownames(all_bar), gen_names),]
  
  # Avg expression and std error per cluster for selected genes
  clusters_rm = ""
  cluster_names = ""
  gene_matrix_names = ""
  counter = 1
  for (i in col.scale){
    
    mat.cluster = all_bar[, colnames(all_bar) %in% names(col.spatial.clusters[col.spatial.clusters == i,])]
    mat.clusters = matrix(ncol = 1, nrow = length(gen_names))
    for (x in unique(sapply(strsplit(colnames(mat.cluster), "_"), "[[",1))){
      mat.clusters = cbind(mat.clusters, rowMeans(mat.cluster[,sapply(strsplit(colnames(mat.cluster), "_"), "[[",1) == x]))
    }
    mat.clusters = mat.clusters[,-1]
    colnames(mat.clusters) = unique(sapply(strsplit(colnames(mat.cluster), "_"), "[[",1))
    # mat.clusters = as.data.frame(rowMeans(mat.clusters))
    clusters_rm_tmp = ""
    gene_matrix = ""
    for (j in colnames(mat.clusters)){
      clusters_rm_tmp = c(clusters_rm_tmp, mat.clusters[,j])
      gene_matrix = c(gene_matrix, row.names(mat.clusters))
    }
    clusters_rm_tmp = clusters_rm_tmp[-1]
    gene_matrix = gene_matrix[-1]
    clusters_rm = c(clusters_rm, clusters_rm_tmp)
    gene_matrix_names = c(gene_matrix_names, gene_matrix)
    cluster_names = c(cluster_names, rep(paste0("Cluster",counter), length(row.names(mat.clusters))*ncol(mat.clusters)))
    
    counter=counter+1
  }
  
  clusters_rm = as.numeric(clusters_rm[-1])
  cluster_names = cluster_names[-1]
  gene_matrix_names = gene_matrix_names[-1]
  myData = data.frame(clusters_rm, cluster_names, gene_matrix_names)

  anno_df = compare_means(clusters_rm ~ cluster_names, data = myData,
                          group.by = "gene_matrix_names", ref.group = "Cluster1",
                          method = "t.test", p.adjust.method = "BH",
                          paired = F, var.equal = TRUE, exact=FALSE) %>% mutate(y_pos = 8)
  anno_df = anno_df[anno_df$group1 == "Cluster1",]
  anno_df$y_position = as.numeric(rep(seq(round(max(myData$clusters_rm)), round(max(myData$clusters_rm))+0.25*(length(unique(cluster_names))-2), by = 0.25), length(unique(gene_matrix_names))))
  anno_df$group = rep(1:length(unique(gene_matrix_names)), each = length(unique(cluster_names))-1)

  signs_p= ""
  for (i in anno_df$p.adj){
    if (i>0.05) signs_p=c(signs_p, "ns")
    if ((i<=0.05) & (i>0.01)) signs_p=c(signs_p, "*")
    if ((i<=0.01) & (i>0.001)) signs_p=c(signs_p, "**")
    if (i<=0.001) signs_p=c(signs_p, "***")
  }
  signs_p = signs_p[-1]
  anno_df$p.signif.adj = signs_p

  myData$gene_matrix_names <- factor(myData$gene_matrix_names, levels = unique(gen_names)) 
  
  anno_df_f = matrix(nrow = 1, ncol = ncol(anno_df))
  colnames(anno_df_f) = colnames(anno_df)
  for (gene in gen_names){
    anno_df_f = rbind(anno_df_f, anno_df[anno_df$gene_matrix_names == gene,])
  }
  anno_df_f = anno_df_f[-1,]

  g <- ggplot(myData, aes(x = gene_matrix_names, y = clusters_rm, fill = cluster_names), group = dummy) + 
    #geom_boxplot() + 
    scale_fill_manual(values = col.scale, aesthetics = "fill")  + 
    stat_boxplot_custom(qs = c(0.0, 0.25, 0.5, 0.75, 1))
  
  xmins = as.numeric(ggplot_build(g)$data[[1]]$xmin) + (as.numeric(ggplot_build(g)$data[[1]]$xmax[1])-as.numeric(ggplot_build(g)$data[[1]]$xmin[1]))/2
  xmaxs = as.numeric(ggplot_build(g)$data[[1]]$xmax) + (as.numeric(ggplot_build(g)$data[[1]]$xmax[1])-as.numeric(ggplot_build(g)$data[[1]]$xmin[1]))/2
  xmins = xmins[!xmins %in%  xmins[seq(1, length(xmins), length(unique(cluster_names)))]]-(as.numeric(ggplot_build(g)$data[[1]]$xmax[1])-as.numeric(ggplot_build(g)$data[[1]]$xmin[1]))
  xmaxs = xmaxs[!xmaxs %in%  xmaxs[seq(1, length(xmaxs), length(unique(cluster_names)))]]-(as.numeric(ggplot_build(g)$data[[1]]$xmax[1])-as.numeric(ggplot_build(g)$data[[1]]$xmin[1]))
  
  p <- g +
    geom_point(shape=16, position = position_jitterdodge(0.01)) + 
    geom_signif(
    xmin = rep(xmins[seq(1, length(xmins), (length(unique(cluster_names))-1))], each = (length(unique(cluster_names))-1)), xmax = xmaxs, y_position = anno_df_f$y_position, annotations = anno_df_f$p.signif.adj,
    tip_length = 0.01) +
    theme(axis.line.y = element_line(colour = "black"), axis.title.x = element_blank(), axis.ticks.x= element_blank(), axis.title.y = element_text(), axis.text.x = element_text(angle = 0),panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())

  return(p)
}

avg_genes_box = function(all_bar, gen_names, cluster_names){
  
  
  #### Make barplot of avg expression per interesting marker genes
  col.scale = unique(cluster_names[,1])
  col.spatial.clusters = cluster_names
  
  ### arrange matrix per cluster 
  
  all_bar = all_bar[row.names(all_bar) %in% gen_names,]
  all_bar = all_bar[,match(colnames(all_bar), rownames(col.spatial.clusters))]
  all_bar = all_bar[match(rownames(all_bar), gen_names),]
  
  # Avg expression and std error per cluster for selected genes
  clusters_rm = ""
  cluster_names = ""
  gene_matrix_names = ""
  counter = 1
  for (i in col.scale){
    
    mat.cluster = all_bar[, colnames(all_bar) %in% names(col.spatial.clusters[col.spatial.clusters == i,])]
    mat.clusters = matrix(ncol = 1, nrow = length(gen_names))
    for (x in unique(sapply(strsplit(colnames(mat.cluster), "_"), "[[",1))){
      mat.clusters = cbind(mat.clusters, rowMeans(mat.cluster[,sapply(strsplit(colnames(mat.cluster), "_"), "[[",1) == x]))
    }
    mat.clusters = mat.clusters[,-1]
    colnames(mat.clusters) = unique(sapply(strsplit(colnames(mat.cluster), "_"), "[[",1))
    # mat.clusters = as.data.frame(rowMeans(mat.clusters))
    clusters_rm_tmp = ""
    gene_matrix = ""
    for (j in colnames(mat.clusters)){
      clusters_rm_tmp = c(clusters_rm_tmp, mat.clusters[,j])
      gene_matrix = c(gene_matrix, row.names(mat.clusters))
    }
    clusters_rm_tmp = clusters_rm_tmp[-1]
    gene_matrix = gene_matrix[-1]
    clusters_rm = c(clusters_rm, clusters_rm_tmp)
    gene_matrix_names = c(gene_matrix_names, gene_matrix)
    cluster_names = c(cluster_names, rep(paste0("Cluster",counter), length(row.names(mat.clusters))*ncol(mat.clusters)))
    
    counter=counter+1
  }
  
  clusters_rm = as.numeric(clusters_rm[-1])
  cluster_names = cluster_names[-1]
  gene_matrix_names = gene_matrix_names[-1]
  myData = data.frame(clusters_rm, cluster_names, gene_matrix_names)

  anno_df = avg_genes_barplot_data(RA.norm, gen_names, col.spatial.clusters)
  anno_df = anno_df[anno_df$group1 != anno_df$group2,] 
  anno_df = anno_df[order(anno_df$gene_names),]
  anno_df$y_position = as.numeric(rep(seq(round(max(myData$clusters_rm)), round(max(myData$clusters_rm))+0.25*(length(unique(cluster_names))-2), by = 0.25), length(unique(gene_matrix_names))))
  anno_df$group = rep(1:length(unique(gene_matrix_names)), each = length(unique(cluster_names))-1)
  
  myData$gene_matrix_names <- factor(myData$gene_matrix_names, levels = unique(anno_df$gene_names)) 
  
  g <- ggplot(myData, aes(x = gene_matrix_names, y = clusters_rm, fill = cluster_names)) + 
    #geom_boxplot() + 
    scale_fill_manual(values = col.scale, aesthetics = "fill")  + 
    stat_boxplot_custom(qs = c(0.0, 0.25, 0.5, 0.75, 1))
  
  xmins = as.numeric(ggplot_build(g)$data[[1]]$xmin) + (as.numeric(ggplot_build(g)$data[[1]]$xmax[1])-as.numeric(ggplot_build(g)$data[[1]]$xmin[1]))/2
  xmaxs = as.numeric(ggplot_build(g)$data[[1]]$xmax) + (as.numeric(ggplot_build(g)$data[[1]]$xmax[1])-as.numeric(ggplot_build(g)$data[[1]]$xmin[1]))/2
  xmins = xmins[!xmins %in%  xmins[seq(1, length(xmins), 4)]]-(as.numeric(ggplot_build(g)$data[[1]]$xmax[1])-as.numeric(ggplot_build(g)$data[[1]]$xmin[1]))
  xmaxs = xmaxs[!xmaxs %in%  xmaxs[seq(1, length(xmaxs), 4)]]-(as.numeric(ggplot_build(g)$data[[1]]$xmax[1])-as.numeric(ggplot_build(g)$data[[1]]$xmin[1]))
  
  p <- g +
    geom_jitter(shape=16, position=position_jitter(0.01)) + 
    geom_signif(
    xmin = rep(xmins[seq(1, length(xmins), 3)], each = 3), xmax = xmaxs, y_position = anno_df$y_position, annotations = anno_df$signs_p,
    tip_length = 0.01) +
    theme(axis.line.y = element_line(colour = "black"), axis.title.x = element_blank(), axis.ticks.x= element_blank(), axis.title.y = element_text(), axis.text.x = element_text(angle = 0),panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  return(p)
  print(mydata)
}

avg_genes_box_data = function(all_bar, gen_names, cluster_names){
  
  
  #### Make barplot of avg expression per interesting marker genes
  col.scale = unique(cluster_names[,1])
  col.spatial.clusters = cluster_names
  
  ### arrange matrix per cluster 
  #all_bar = RA.norm
  all_bar = all_bar[row.names(all_bar) %in% gen_names,]
  all_bar = all_bar[,match(colnames(all_bar), rownames(col.spatial.clusters))]
  all_bar = all_bar[match(rownames(all_bar), gen_names),]
  
  # Avg expression and std error per cluster for selected genes
  clusters_rm = ""
  cluster_names = ""
  gene_matrix_names = ""
  counter = 1
  for (i in col.scale){
    
    mat.cluster = all_bar[, colnames(all_bar) %in% names(col.spatial.clusters[col.spatial.clusters == i,])]
    mat.clusters = matrix(ncol = 1, nrow = length(gen_names))
    for (x in unique(sapply(strsplit(colnames(mat.cluster), "_"), "[[",1))){
      mat.clusters = cbind(mat.clusters, rowMeans(mat.cluster[,sapply(strsplit(colnames(mat.cluster), "_"), "[[",1) == x]))
    }
    mat.clusters = mat.clusters[,-1]
    colnames(mat.clusters) = unique(sapply(strsplit(colnames(mat.cluster), "_"), "[[",1))
    # mat.clusters = as.data.frame(rowMeans(mat.clusters))
    clusters_rm_tmp = ""
    gene_matrix = ""
    for (j in colnames(mat.clusters)){
      clusters_rm_tmp = c(clusters_rm_tmp, mat.clusters[,j])
      gene_matrix = c(gene_matrix, row.names(mat.clusters))
    }
    clusters_rm_tmp = clusters_rm_tmp[-1]
    gene_matrix = gene_matrix[-1]
    clusters_rm = c(clusters_rm, clusters_rm_tmp)
    gene_matrix_names = c(gene_matrix_names, gene_matrix)
    cluster_names = c(cluster_names, rep(paste0("Cluster",counter), length(row.names(mat.clusters))*ncol(mat.clusters)))
    
    counter=counter+1
  }
  
  clusters_rm = as.numeric(clusters_rm[-1])
  cluster_names = cluster_names[-1]
  gene_matrix_names = gene_matrix_names[-1]
  myData = data.frame(clusters_rm, cluster_names, gene_matrix_names)
  
  anno_df = avg_genes_barplot_data(RA.norm, gen_names, col.spatial.clusters)
  anno_df = anno_df[anno_df$group1 != anno_df$group2,] 
  anno_df = anno_df[order(anno_df$gene_names),]
  anno_df$y_position = as.numeric(rep(seq(round(max(myData$clusters_rm)), round(max(myData$clusters_rm))+0.25*(length(unique(cluster_names))-2), by = 0.25), length(unique(gene_matrix_names))))
  anno_df$group = rep(1:length(unique(gene_matrix_names)), each = length(unique(cluster_names))-1)
  
  myData$gene_matrix_names <- factor(myData$gene_matrix_names, levels = unique(anno_df$gene_names)) 
  
  return(mydata)

}

avg_genes_box_old_data = function(all_bar, gen_names, cluster_names){
  
  #### Make barplot of avg expression per interesting marker genes
  col.scale = unique(cluster_names[,1])
  col.spatial.clusters = cluster_names
  
  ### arrange matrix per cluster 
  all_bar = all_bar[row.names(all_bar) %in% gen_names,]
  all_bar = all_bar[,match(colnames(all_bar), rownames(cluster_names))]
  all_bar = all_bar[match(rownames(all_bar), gen_names),]
  
  # Avg expression and std error per cluster for selected genes
  clusters_rm = ""
  cluster_names = ""
  gene_matrix_names = ""
  counter = 1
  for (i in col.scale){
    
    mat.cluster = all_bar[, colnames(all_bar) %in% names(col.spatial.clusters[col.spatial.clusters == i,])]
    mat.clusters = matrix(ncol = 1, nrow = length(gen_names))
    for (x in unique(sapply(strsplit(colnames(mat.cluster), "_"), "[[",1))){
      mat.clusters = cbind(mat.clusters, rowMeans(mat.cluster[,sapply(strsplit(colnames(mat.cluster), "_"), "[[",1) == x]))
    }
    mat.clusters = mat.clusters[,-1]
    colnames(mat.clusters) = unique(sapply(strsplit(colnames(mat.cluster), "_"), "[[",1))
    # mat.clusters = as.data.frame(rowMeans(mat.clusters))
    clusters_rm_tmp = ""
    gene_matrix = ""
    for (j in colnames(mat.clusters)){
      clusters_rm_tmp = c(clusters_rm_tmp, mat.clusters[,j])
      gene_matrix = c(gene_matrix, row.names(mat.clusters))
    }
    clusters_rm_tmp = clusters_rm_tmp[-1]
    gene_matrix = gene_matrix[-1]
    clusters_rm = c(clusters_rm, clusters_rm_tmp)
    gene_matrix_names = c(gene_matrix_names, gene_matrix)
    cluster_names = c(cluster_names, rep(paste0("Cluster",counter), length(row.names(mat.clusters))*ncol(mat.clusters)))
    
    counter=counter+1
  }
  
  clusters_rm = as.numeric(clusters_rm[-1])
  cluster_names = cluster_names[-1]
  gene_matrix_names = gene_matrix_names[-1]
  myData = data.frame(clusters_rm, cluster_names, gene_matrix_names)
  # for (clus in unique(cluster_names)){
  #   for (gen in unique(gene_matrix_names)){
  #     print(paste0("checking normal distr for: ", clus, " and ", gen))
  #     print(with(myData, shapiro.test(clusters_rm[(cluster_names == clus) & (gene_matrix_names == gen)])))
  #   }
  # }
  # m = myData[(gene_matrix_names == "CD52") & ((cluster_names == "Cluster1") | (cluster_names == "Cluster2")), ]
  # print(m)
  # r = var.test(clusters_rm ~ cluster_names, data = m)
  # print(r)
  
  anno_df = compare_means(clusters_rm ~ cluster_names, data = myData,
                          group.by = "gene_matrix_names", ref.group = "Cluster1",
                          method = "t.test", p.adjust.method = "BH",
                          paired = F, var.equal = TRUE, exact=FALSE) %>% mutate(y_pos = 8) #, alternative = c("two-sided")
  anno_df = anno_df[anno_df$group1 == "Cluster1",]
  anno_df$y_position = as.numeric(rep(seq(round(max(myData$clusters_rm)), round(max(myData$clusters_rm))+0.25*(length(unique(cluster_names))-2), by = 0.25), length(unique(gene_matrix_names))))
  anno_df$group = rep(1:length(unique(gene_matrix_names)), each = length(unique(cluster_names))-1)

  signs_p= ""
  for (i in anno_df$p.adj){
    if (i>0.05) signs_p=c(signs_p, "ns")
    if ((i<=0.05) & (i>0.01)) signs_p=c(signs_p, "*")
    if ((i<=0.01) & (i>0.001)) signs_p=c(signs_p, "**")
    if (i<=0.001) signs_p=c(signs_p, "***")
  }
  signs_p = signs_p[-1]
  anno_df$p.signif.adj = signs_p
  return(anno_df)

}

scale_this <-  function(x, na.rm = TRUE) {
  return((x- min(x)) /(max(x)-min(x)))
}
