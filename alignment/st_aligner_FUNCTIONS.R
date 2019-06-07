new.coordinates = function(spot.image, spots, sample, well){

  require(data.table)
  test.y = attributes(spot.image)$dim[1] / 34
  test.x = attributes(spot.image)$dim[2] / 32
  rm(spot.image)

  V1 = (spots$X)/ (test.x) + 1
  V2 = (spots$Y)/ (test.y)  + 1
  img = cbind(V1, V2)

  bc.x = matrix(nrow=nrow(img))
  bc.x = bc.x[complete.cases(bc.x),]

  for (bac in img[,1]) {
    if (bac%%1 <= 0.5)
      bc.x = rbind(bc.x,round(bac, digits = 0))
    else bc.x = rbind(bc.x,ceiling(bac))}
  bc.x = bc.x[!is.na(bc.x)]


  bc.y = matrix(nrow=nrow(img))
  bc.y = bc.y[complete.cases(bc.y),]

  for (bac in img[,2]) {
    if (bac%%1 <= 0.5)
      bc.y = rbind(bc.y,round(bac, digits = 0))
    else bc.y = rbind(bc.y,ceiling(bac))}
  bc.y = bc.y[!is.na(bc.y)]

  # r.coords.x = bc.x$V1
  # r.coords.y = bc.y$V2
  
  r.coords.x = bc.x
  r.coords.y = bc.y

  rt.coords.x = round(img[,1],digits = 3)
  rt.coords.y = round(img[,2],digits = 3)

  new.barcode.names = paste(r.coords.x,r.coords.y, sep="_")
  selected.barcodes = cbind(new.barcode.names, spots$X, spots$Y)

  return(selected.barcodes)

}

write.xy.matrix = function(exp.values) {
  write.table(exp.values, file=paste(sample, "_", well,   "_xy.tsv",sep=""), sep = "\t", quote=F, col.names=T, row.names=T)
}



