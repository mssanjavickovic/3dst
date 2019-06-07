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
  data = t(as.matrix(inf_all))
  s = find_sigmas(data, step = 0.01, steps = 100, sample_rows = 500)
  opt.s = optimal_sigma(s)
  ts <- Transitions(data, sigma = opt.s)
  ev <- base::eigen(as.matrix(ts@transitions), TRUE)$vectors
  dm <- as.data.frame(ev[, -1])
  colnames(dm) <- paste0('RA', seq_len(ncol(dm)))
  rownames(dm) = rownames(data)
  return(dm)
}

# Clean up some gene names; mat is a matrix with gene names as row names
genes_cleanup = function (mat){
  mat = mat[-grep(pattern="RPS", x= rownames(mat), value = F),]
  mat = mat[-grep(pattern="RPL", x= rownames(mat), value = F),]
  mat = mat[-grep(pattern="ambiguous", x= rownames(mat), value = F),]
  mat = mat[-grep(pattern="MTRNR", x= rownames(mat), value = F),]
  mat = mat[-grep(pattern="ENSG", x= rownames(mat), value = F),]
  return(mat)
}

### Plot 3D Heatmaps
# Plot gene specific stuff based on number of sections in each sample
plot.gene.3d.9 = function(sample, gene, m1, m2, m3, m4, m5, m6, m7, m8, m9, s1, s2, s3, s4, s5, s6, s7, s8, s9, x, y, transparency, min, max){
  
  for (i in c(1:9)) {
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
    if (i == "8") {
      spots = s8
      exp.values = m8}
    if (i == "9") {
      spots = s9
      exp.values = m9}
    
    # spots = s10
    # exp.values = m6
    # gene = "CD4"
    print(i)
    
    genes.barcodes = exp.values[rownames(exp.values) == gene,]
    #rownames(spots) = paste("X",rownames(spots),sep="")
    exp.values = exp.values[,colnames(exp.values) %in% rownames(spots)]
    exp.values = exp.values[rowSums(exp.values) != 0,]
    spots = spots[colnames(exp.values),]
    genes.barcodes = cbind(spots, genes.barcodes)
    
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
    if (i == "5") {mat.5 = s$z}
    if (i == "6") {mat.6 = s$z}
    if (i == "7") {mat.7 = s$z}
    if (i == "8") {mat.8 = s$z}
    if (i == "9") {mat.9 = s$z}}

  png(paste(sample,"_",gene,"_clusters.png", sep=""))
  par(mar=c(2.1, 6.1, 15.1, 2.1), xpd=TRUE)
  # use http://colorbrewer2.org/ to find optimal divergent color palette (or set own)

  image3D(z = 0, phi=37, theta = -185, colvar = mat.1, smooth = TRUE, alpha = transparency*1, box = FALSE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F)
  image3D(z = 5, phi=37,  theta = -185, colvar = mat.2,   smooth = TRUE, alpha = transparency*1, add = TRUE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F)
  image3D(z = 10, phi=37,  theta = -185, colvar = mat.3,  smooth = TRUE, alpha = transparency*1, add = TRUE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F)
  image3D(z = 14, phi=37,  theta = -185, colvar = mat.4,  smooth = TRUE, alpha = transparency*1, add = TRUE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F)
  image3D(z = 18, phi=37,  theta = -185, colvar = mat.5,  smooth = TRUE, alpha = transparency*1, add = TRUE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F)
  image3D(z = 22, phi=37,  theta = -185, colvar = mat.6,  smooth = TRUE, alpha = transparency*1, add = TRUE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F)
  image3D(z = 25, phi=37,  theta = -185, colvar = mat.7,  smooth = TRUE, alpha = transparency*1, add = TRUE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F)
  image3D(z = 28, phi=37,  theta = -185, colvar = mat.8,  smooth = TRUE, alpha = transparency*1, add = TRUE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F)
  image3D(z = 30.5, phi=37,  theta = -185, colvar = mat.9,  smooth = TRUE, alpha = transparency*1, add = TRUE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F)
  dev.off()
}

### Plot cluster specific stuff
plot.gene.3d.cluster.9 = function(sample, cluster, m1, m2, m3, m4, m5, m6, m7, m8, m9, s1, s2, s3, s4, s5, s6, s7, s8, s9, x, y, transparency, min, max){
  
  for (i in c(1:9)) {
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
    if (i == "8") {
      spots = s8
      exp.values = m8}
    if (i == "9") {
      spots = s9
      exp.values = m9}
    
    clust.spots = rbind(s1,s2,s3,s4,s5,s6,s7,s8,s9)
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
    if (i == "5") {mat.5 = s$z}
    if (i == "6") {mat.6 = s$z}
    if (i == "7") {mat.7 = s$z}
    if (i == "8") {mat.8 = s$z}
    if (i == "9") {mat.9 = s$z}}
  
  
  png(paste(sample,"_clusters.png", sep=""))
  par(mar=c(2.1, 6.1, 15.1, 2.1), xpd=TRUE)
  colorvar = colorRampPalette(c("#FF7F00","#FED8B1","gray61", "gray76", "gray85"))(5)
  image3D(z = 0, phi=37, theta = -185, colvar = mat.1, col=colorvar, smooth = TRUE, alpha = transparency*1, box = FALSE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F)
  image3D(z = 5, phi=37,  theta = -185, colvar = mat.2,  col=colorvar, smooth = TRUE, alpha = transparency*1, add = TRUE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F)
  image3D(z = 10, phi=37,  theta = -185, colvar = mat.3,  col=colorvar, smooth = TRUE, alpha = transparency*1, add = TRUE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F)
  image3D(z = 14, phi=37,  theta = -185, colvar = mat.4,  col=colorvar, smooth = TRUE, alpha = transparency*1, add = TRUE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F)
  image3D(z = 18, phi=37,  theta = -185, colvar = mat.5,  col=colorvar, smooth = TRUE, alpha = transparency*1, add = TRUE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F)
  image3D(z = 22, phi=37,  theta = -185, colvar = mat.6,  col=colorvar, smooth = TRUE, alpha = transparency*1, add = TRUE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F)
  image3D(z = 25, phi=37,  theta = -185, colvar = mat.7,  col=colorvar, smooth = TRUE, alpha = transparency*1, add = TRUE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F)
  image3D(z = 28, phi=37,  theta = -185, colvar = mat.8,  col=colorvar, smooth = TRUE, alpha = transparency*1, add = TRUE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F)
  image3D(z = 30.5, phi=37,  theta = -185, colvar = mat.9,  col=colorvar, smooth = TRUE, alpha = transparency*1, add = TRUE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F)
  dev.off()

}

### Plot 3D Heatmaps
### Plot gene specific stuff
plot.gene.3d.4 = function(sample, gene, m1, m2, m3, m4, s1, s2, s3, s4,x, y, transparency, min, max){
  
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
    
    rownames(spots) = paste("X",rownames(spots),sep="")
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
    i=4
    s =  interp(x1,y1, z, nx = x, ny = y)
    if (i == "1") {mat.1 = s$z}
    if (i == "2") {mat.2 = s$z}
    if (i == "3") {mat.3 = s$z}
    if (i == "4") {mat.4 = s$z}}
  
  png(paste(sample,"_",gene,"_clusters.png", sep=""))
  par(mar=c(2.1, 6.1, 15.1, 2.1), xpd=TRUE)

  image3D(z = 0, phi=25, theta = 10, colvar = mat.1, smooth = TRUE, alpha = transparency*1, box = FALSE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F)
  image3D(z = 4, phi=25,  theta = 10, colvar = mat.2, smooth = TRUE, alpha = transparency*1, add = TRUE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F)
  image3D(z = 8, phi=25,  theta = 10, colvar = mat.3, smooth = TRUE, alpha = transparency*1, add = TRUE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F)
  image3D(z = 11, phi=25,  theta =10, colvar = mat.4, smooth = TRUE, alpha = transparency*2, add = TRUE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F)
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

  
  png(paste(sample,"_clusters.png", sep=""))
  par(mar=c(2.1, 6.1, 15.1, 2.1), xpd=TRUE)
  colorvar = colorRampPalette(c("#FF7F00","#FED8B1","gray61", "gray76", "gray85"))(5)
  image3D(z = 0, phi=25, theta = 10, colvar = mat.1, col=colorvar, smooth = TRUE, alpha = transparency*0.25, box = FALSE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F)
  image3D(z = 4, phi=25,  theta = 10, colvar = mat.2,  col=colorvar, smooth = TRUE, alpha = transparency*0.5, add = TRUE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F)
  image3D(z = 8, phi=25,  theta = 10, colvar = mat.3,  col=colorvar, smooth = TRUE, alpha = transparency*0.75, add = TRUE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F)
  image3D(z = 11, phi=25,  theta = 10, colvar = mat.4,  col=colorvar, smooth = TRUE, alpha = transparency*2, add = TRUE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F)
  dev.off()
  
}


### Plot 2D clusters
#Funciton
plot.gene.2d.cluster.4 = function(sample, cluster, m1, m2, m3, m4, s1, s2, s3, s4, x, y, transparency, min, max){
  
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
    
    library(akima)
    # x= 40
    # y= 40
    # transparency=1
    # min=1
    # max=5
    s =  interp(x1, y1, z, nx = x, ny = y)
    if (i == "1") {mat.1 = s$z}
    if (i == "2") {mat.2 = s$z}
    if (i == "3") {mat.3 = s$z}
    if (i == "4") {mat.4 = s$z}}

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
}
  
plot.gene.2d.4 = function(sample, gene,  m1, m2, m3, m4, s1, s2, s3, s4, x, y, transparency, min, max, con){
  
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
}


plot.gene.2d.cluster.7 = function(sample, cluster, m1, m2, m3, m4, m5, m6, m7, s1, s2, s3, s4, s5, s6, s7, x, y, transparency, min, max){

  
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
    rownames(clust.spots) = paste("X",rownames(clust.spots),sep="")
    genes.barcodes = clust.spots[rownames(clust.spots) %in% colnames(exp.values),]
    exp.values = exp.values[,colnames(exp.values) %in% rownames(clust.spots)]
    exp.values = exp.values[rowSums(exp.values) != 0,]
    cluster = ann_to_plot
    cluster = cluster[colnames(exp.values),]
    genes.barcodes = cbind(genes.barcodes, cluster)
    
    x1 = as.numeric(genes.barcodes[,1])
    x1 = x1[!is.na(x1)]
    y1 = as.numeric(genes.barcodes[,2])
    y1 = y1[!is.na(y1)]
    z = as.numeric(genes.barcodes[,3])
    dis = rep(0.25 * as.numeric(i), nrow(genes.barcodes))
    w = as.numeric(dis)
    
    library(akima)
    s =  interp(x1, y1, z, nx = x, ny = y)
    if (i == "1") {mat.1 = s$z}
    if (i == "2") {mat.2 = s$z}
    if (i == "3") {mat.3 = s$z}
    if (i == "4") {mat.4 = s$z}
    if (i == "5") {mat.5 = s$z}
    if (i == "6") {mat.6 = s$z}
    if (i == "7") {mat.7 = s$z}
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
  pdf(paste(sample,"_6_Contour_Cluster",".pdf", sep=""))
  image2D(z = mat.6, contour = T, smooth = TRUE, col=colorvar, alpha = transparency, box = FALSE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F)
  dev.off()
  pdf(paste(sample,"_7_Contour_Cluster",".pdf", sep=""))
  image2D(z = mat.7, contour = T, smooth = TRUE, col=colorvar, alpha = transparency, box = FALSE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F)
  dev.off()
}

plot.gene.2d.7 = function(sample, gene, m1, m2, m3, m4, m5, m6, m7, s1, s2, s3, s4, s5, s6, s7, x, y, transparency, min, max, con){
  
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
    
    # row.names(mf1)
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
    
    print(i)
    
    x1 = as.numeric(genes.barcodes.1[,1])
    x1 = x1[!is.na(x1)]
    y1 = as.numeric(genes.barcodes.1[,2])
    y1 = y1[!is.na(y1)]
    z = as.numeric(genes.barcodes.1[,3])

    s =  interp(x1, y1, z, nx = x, ny = y)
    if (i == "1") {mat.1 = s$z}
    if (i == "2") {mat.2 = s$z}
    if (i == "3") {mat.3 = s$z}
    if (i == "4") {mat.4 = s$z}
    if (i == "5") {mat.5 = s$z}
    if (i == "6") {mat.6 = s$z}
    if (i == "7") {mat.7 = s$z}
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
  
}

plot.gene.2d.inf.7 = function(sample, cluster, m1, m2, m3, m4, m5, m6, m7, s1, s2, s3, s4, s5, s6, s7, x, y, transparency, min, max){
  
  
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
    library(akima)
    s =  interp(x1, y1, z, nx = x, ny = y)
    if (i == "1") {mat.1 = s$z}
    if (i == "2") {mat.2 = s$z}
    if (i == "3") {mat.3 = s$z}
    if (i == "4") {mat.4 = s$z}
    if (i == "5") {mat.5 = s$z}
    if (i == "6") {mat.6 = s$z}
    if (i == "7") {mat.7 = s$z}
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
  pdf(paste(sample,"_6_Contour_Inf",".pdf", sep=""))
  image2D(z = mat.6, contour = T, smooth = TRUE, col=colorvar, alpha = transparency, box = FALSE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F)
  dev.off()
  pdf(paste(sample,"_7_Contour_Inf",".pdf", sep=""))
  image2D(z = mat.7, contour = T, smooth = TRUE, col=colorvar, alpha = transparency, box = FALSE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F)
  dev.off()

}

plot.gene.2d.inf.4 = function(sample, cluster, m1, m2, m3, m4, s1, s2, s3, s4, x, y, transparency, min, max){
    
    
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
        library(akima)
        s =  interp(x1, y1, z, nx = x, ny = y)
        if (i == "1") {mat.1 = s$z}
        if (i == "2") {mat.2 = s$z}
        if (i == "3") {mat.3 = s$z}
        if (i == "4") {mat.4 = s$z}
    }
    
    colorvar = colorRampPalette(c("#4ebceb","#f15f48","gray80"))(3)
    
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
    
}
