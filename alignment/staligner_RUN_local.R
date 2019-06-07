# Load libraries 
require(jpeg)
require(raster)
require(scales)
require(gdata)

# Expect command line args at the end. 
args = commandArgs(trailingOnly = TRUE)

# Source the functions file (make sure path is valid)
source('st_aligner_FUNCTIONS.R')

# You need to input correct name/files here
all_barcodes <- read.delim(args[1], sep="\t",quote="",row.names = 1) # this is to run from cmd line

# Set names
sample = sapply(strsplit(args[1],split="_"),"[[",1)
print(sample) # qc step

well = sapply(strsplit(args[1],split="_"),"[[",2)
print(well) # qc step

# Reads in image files and a Fiji tsv file
spot.image = readJPEG(paste(sample, "_", well,'.jpg', sep="")) # Can be HE or spots file
spots = read.table(paste(sample,"_", well,'.jpg_spots.tsv', sep="")) # file created in Fiji with centroid pixel coordinates of each spot under the tissue

# Get new coordinates and check alignment
colnames(all_barcodes) = sapply(strsplit(colnames(all_barcodes),split="X"),"[[",2)
selected.barcodes = data.frame(new.coordinates(spot.image, spots, sample, well))
cn = data.frame(colnames(all_barcodes))
xt = merge(selected.barcodes, cn, by = 1)
xt = xt[!duplicated(xt$new.barcode.names),]
row.names(xt) = xt[,1]
xt = xt[,-1]
colnames(xt) = c("pix_x", "pix_y")

# Test that the transformation is correct
png(paste(sample, well, "QC_alignment_check.png", sep="_"))
xa = sapply(strsplit(as.character(row.names(xt)), split = "_"), "[[",1)
ya = sapply(strsplit(as.character(row.names(xt)), split = "_"), "[[",2)
plot(as.numeric(xa), as.numeric(ya), cex = 2, ylim=c(35,1), xlim=c(1,33),ylab = '', xlab = '', xaxt='n',   yaxt='n', bty = 'n', col = c("black"), pch = 16)
xa = sapply(strsplit(colnames(all_barcodes), split = "_"), "[[",1)
ya = sapply(strsplit(colnames(all_barcodes), split = "_"), "[[",2)
par(new=T)
plot(as.numeric(xa), as.numeric(ya), cex = 2, ylim=c(35,1), xlim=c(1,33),ylab = '', xlab = '', xaxt='n',   yaxt='n', bty = 'n', col = c("red"), pch = 1)
dev.off()

### Write out aligned coordinates with (x,y) and (pix_x, pix_y) file
write.xy.matrix(xt)
