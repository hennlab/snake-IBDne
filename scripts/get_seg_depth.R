# Get segment depth
# Author: SSG
#------------
# This script accepts a file listing segments (i.e. PLINK-generated RoH file, GERMLINE2 generated .match file, IBDSeq/Beagle generated .ibd file) and a .bim file
# It returns the depth of RoH at every SNP position in the .bim file
# The required arguments are: Input (full file names)
# Ex.Rscript get_RoH_depth.R CHABU.Ethiopian_MEGA_RoH_segsJoined.hom CHABU.Ethiopian_MEGA_HWE-filtered.bim

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# Check for correct number of arguments
if (length(args)<2) {
  stop("At least two arguments must be supplied", call.=FALSE)
} else if (length(args)>2) {
  stop("Too many arguments", call.=FALSE)
}

# Read in data
options(stringsAsFactors=F)

fileName <- unlist(strsplit(args[1], split=".", fixed=TRUE))
ext <- fileName[length(fileName)]
if (ext == "hom") {
  seg.file <- read.table(args[1], header=T)
} else if (ext == "ibd") {
  seg.file <- read.table(args[1], header=F)
  colnames(seg.file)[5:7] <- c("CHR", "POS1", "POS2")
} else if (ext == "match") {
  seg.file <- read.table(args[1], header=F)
  if (ncol(seg.file==13)) {
    colnames(seg.file)[3:5] <- c("CHR", "POS1", "POS2")
  } else if (ncol(seg.file==15)) {
    colnames(seg.file)[5:7] <- c("CHR", "POS1", "POS2")
  }
}

bim <- read.table(args[2])

depth.counts <- as.data.frame(matrix(ncol=3, nrow=nrow(bim)))
depth.counts <- cbind(bim[,c(1,4)],0)
  colnames(depth.counts) <- c("CHR","BP","COUNT")

n.segs <- nrow(seg.file) 
for (n in 1:n.segs) {
  segment <- seg.file[n,]
  to.update <- which(bim$V4 >= segment$POS1 & bim$V4 <= segment$POS2 & bim$V1==segment$CHR)
  depth.counts[to.update,3] <- depth.counts[to.update,3]+1
  #if (n %% 10000 == 0) {
  #  cat(paste("Depth calculated for ", round(signif(n/n.segs*(100), 3)), "% of segments", "\n", sep=""))
  #}
}

write.table(depth.counts, "", quote=F, row.names=F, col.names=F, sep=" ")
