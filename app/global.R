### Load data
load('data/RA_1/RA1_norm.exp.values.1')
t1.1 = m1
load('data/RA_1/RA1_norm.exp.values.2')
t1.2 = m2
load('data/RA_1/RA1_norm.exp.values.3')
t1.3 = m3
load('data/RA_1/RA1_norm.exp.values.4')
t1.4 = m4

load('data/RA_2/RA2_norm.exp.values.1')
t2.1 = m1
load('data/RA_2/RA2_norm.exp.values.2')
t2.2 = m2
load('data/RA_2/RA2_norm.exp.values.3')
t2.3 = m3
load('data/RA_2/RA2_norm.exp.values.4')
t2.4 = m4
load('data/RA_2/RA2_norm.exp.values.5')
t2.5 = m5
load('data/RA_2/RA2_norm.exp.values.6')
t2.6 = m6
load('data/RA_2/RA2_norm.exp.values.7')
t2.7 = m7

### Load data
load('data/RA_3/RA3_norm.exp.values.1')
t3.1 = m1
load('data/RA_3/RA3_norm.exp.values.2')
t3.2 = m2
load('data/RA_3/RA3_norm.exp.values.3')
t3.3 = m3
load('data/RA_3/RA3_norm.exp.values.4')
t3.4 = m4

### Load data
load('data/RA_4/RA4_norm.exp.values.1')
t4.1 = m1
load('data/RA_4/RA4_norm.exp.values.2')
t4.2 = m2
load('data/RA_4/RA4_norm.exp.values.3')
t4.3 = m3
load('data/RA_4/RA4_norm.exp.values.4')
t4.4 = m4
load('data/RA_4/RA4_norm.exp.values.5')
t4.5 = m5

### Load data
load('data/RA_5/RA5_norm.exp.values.1')
t5.1 = m1
load('data/RA_5/RA5_norm.exp.values.2')
t5.2 = m2
load('data/RA_5/RA5_norm.exp.values.3')
t5.3 = m3

# Names to be used in app
t1_names <- c(row.names(t1.1), row.names(t1.2), row.names(t1.3), rownames(t1.4))
t1_names = unique(t1_names)

# Names to be used in app
t2_names <- c(row.names(t2.1), row.names(t2.2), row.names(t2.3), rownames(t2.4),row.names(t2.5), row.names(t2.6),rownames(t2.7))
t2_names = unique(t2_names)

#Names to be used in app
t3_names <- c(row.names(t3.1), row.names(t3.2), row.names(t3.3), rownames(t3.4))
t3_names = unique(t3_names)

#Names to be used in app
t4_names <- c(row.names(t4.1), row.names(t4.2), row.names(t4.3), rownames(t4.4), rownames(t4.5))
t4_names = unique(t4_names)

#Names to be used in app
t5_names <- c(row.names(t5.1), row.names(t5.2), row.names(t5.3))
t5_names = unique(t5_names)

#Load z spots
load('data/RA_1/RA1_1_selected_adjusted_spots_3D_manual_app')
s1.1 = s1
load('data/RA_1/RA1_2_selected_adjusted_spots_3D_manual_app')
s1.2 = s2
load('data/RA_1/RA1_3_selected_adjusted_spots_3D_manual_app')
s1.3 = s3
load('data/RA_1/RA1_4_selected_adjusted_spots_3D_manual_app')
s1.4 = s4

load('data/RA_2/RA2_1_selected_adjusted_spots_3D_manual_app')
s2.1 = s1
load('data/RA_2/RA2_2_selected_adjusted_spots_3D_manual_app')
s2.2 = s2
load('data/RA_2/RA2_3_selected_adjusted_spots_3D_manual_app')
s2.3 = s3
load('data/RA_2/RA2_4_selected_adjusted_spots_3D_manual_app')
s2.4 = s4
load('data/RA_2/RA2_5_selected_adjusted_spots_3D_manual_app')
s2.5 = s5
load('data/RA_2/RA2_6_selected_adjusted_spots_3D_manual_app')
s2.6 = s6
load('data/RA_2/RA2_7_selected_adjusted_spots_3D_manual_app')
s2.7 = s7

load('data/RA_3/RA3_1_selected_adjusted_spots_3D_manual_app')
s3.1 = s1
load('data/RA_3/RA3_2_selected_adjusted_spots_3D_manual_app')
s3.2 = s2
load('data/RA_3/RA3_3_selected_adjusted_spots_3D_manual_app')
s3.3 = s3
load('data/RA_3/RA3_4_selected_adjusted_spots_3D_manual_app')
s3.4 = s4

load('data/RA_4/RA4_1_selected_adjusted_spots_3D_manual_app')
s4.1 = s1
load('data/RA_4/RA4_2_selected_adjusted_spots_3D_manual_app')
s4.2 = s2
load('data/RA_4/RA4_3_selected_adjusted_spots_3D_manual_app')
s4.3 = s3
load('data/RA_4/RA4_4_selected_adjusted_spots_3D_manual_app')
s4.4 = s4
load('data/RA_4/RA4_5_selected_adjusted_spots_3D_manual_app')
s4.5 = s5

load('data/RA_5/RA5_1_selected_adjusted_spots_3D_manual_app')
s5.1 = s1
load('data/RA_5/RA5_2_selected_adjusted_spots_3D_manual_app')
s5.2 = s2
load('data/RA_5/RA5_3_selected_adjusted_spots_3D_manual_app')
s5.3 = s3

## All names together
t_names = c(t1_names, t2_names, t3_names, t4_names, t5_names)
t_names = unique(t_names)

### Load functions
plot.gene.3d.4 = function(sample, gene, m1, m2, m3, m4, s1, s2, s3, s4, x, y, transparency, min, max){
  
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
  }
  
    image3Drgl(z = 10, colvar = mat.1, smooth = TRUE, alpha = transparency*2, box = FALSE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F)
    image3Drgl(z = 20, colvar = mat.2, smooth = TRUE, alpha = transparency*3, add = TRUE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F)
    image3Drgl(z = 30, colvar = mat.3, smooth = TRUE, alpha = transparency*4, add = TRUE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F)
    image3Drgl(z = 40, colvar = mat.4, smooth = TRUE, alpha = transparency*5, add = TRUE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F)
  #   bgplot3d({
  #   plot.new()
  #   title(main = paste(sample,"_", gene, sep=""), line = 2, cex.main = 1.7, font.main = 2, col.main = c("black"))
  # })
}

plot.gene.3d.9 = function(sample, gene, m1, m2, m3, m4, m5, m6, m7, s1, s2, s3, s4, s5, s6, s7, x, y, transparency, min, max){
  
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

    
    intersecting.spots <- intersect(rownames(spots), colnames(exp.values))
    #genes.barcodes = exp.values[rownames(exp.values) == gene,]
    exp.values = exp.values[ ,intersecting.spots]
    #exp.values = exp.values[rowSums(exp.values) != 0, ]
    genes.barcodes = exp.values[gene, ]
    spots = spots[intersecting.spots, ]
    #print(nrow(spots)); print(length(genes.barcodes))
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

  }
  
  image3Drgl(z = 10, colvar = mat.1, smooth = TRUE, alpha = transparency*2, box = FALSE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F)
  image3Drgl(z = 20, colvar = mat.2, smooth = TRUE, alpha = transparency*2, add = TRUE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F)
  image3Drgl(z = 30, colvar = mat.3, smooth = TRUE, alpha = transparency*2, add = TRUE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F)
  image3Drgl(z = 40, colvar = mat.4, smooth = TRUE, alpha = transparency*3, add = TRUE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F)
  image3Drgl(z = 50, colvar = mat.5, smooth = TRUE, alpha = transparency*3, add = TRUE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F)
  image3Drgl(z = 60, colvar = mat.6, smooth = TRUE, alpha = transparency*3, add = TRUE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F)
  image3Drgl(z = 70, colvar = mat.7, smooth = TRUE, alpha = transparency*5, add = TRUE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F)
  # bgplot3d({
  #   plot.new()
  #   title(main = paste(sample,"_", gene, sep=""), line = 2, cex.main = 1.7, font.main = 2, col.main = c("black"))
  # })
}

plot.gene.3d.5 = function(sample, gene, m1, m2, m3, m4, m5, s1, s2, s3, s4, s5, x, y, transparency, min, max){
  
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
  }
  
  image3Drgl(z = 10, colvar = mat.1, smooth = TRUE, alpha = transparency*1, box = FALSE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F)
  image3Drgl(z = 20, colvar = mat.2, smooth = TRUE, alpha = transparency*2, add = TRUE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F)
  image3Drgl(z = 30, colvar = mat.3, smooth = TRUE, alpha = transparency*3, add = TRUE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F)
  image3Drgl(z = 40, colvar = mat.4, smooth = TRUE, alpha = transparency*4, add = TRUE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F)
  image3Drgl(z = 50, colvar = mat.5, smooth = TRUE, alpha = transparency*5, add = TRUE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F)
   #   bgplot3d({
  #   plot.new()
  #   title(main = paste(sample,"_", gene, sep=""), line = 2, cex.main = 1.7, font.main = 2, col.main = c("black"))
  # })
}

plot.gene.3d.3 = function(sample, gene, m1, m2, m3, s1, s2, s3, x, y, transparency, min, max){
  
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

    genes.barcodes = exp.values[rownames(exp.values) == gene,]
    exp.values = exp.values[,colnames(exp.values) %in% rownames(spots)]
    exp.values= exp.values[rowSums(exp.values) != 0,]
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

  }
  
  image3Drgl(z = 10, colvar = mat.1, smooth = TRUE, alpha = transparency*1, box = FALSE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F)
  image3Drgl(z = 20, colvar = mat.2, smooth = TRUE, alpha = transparency*2, add = TRUE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F)
  image3Drgl(z = 30, colvar = mat.3, smooth = TRUE, alpha = transparency*3, add = TRUE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F)
  #image3Drgl(z = 40, colvar = mat.4, smooth = TRUE, alpha = transparency*4, add = TRUE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F)
  #image3Drgl(z = 50, colvar = mat.5, smooth = TRUE, alpha = transparency*5, add = TRUE, inttype = 1, clim = c(min, max), NAcol = "transparent", x = c(1:x), y = c(1:y), scale = F)
  #   bgplot3d({
  #   plot.new()
  #   title(main = paste(sample,"_", gene, sep=""), line = 2, cex.main = 1.7, font.main = 2, col.main = c("black"))
  # })
}