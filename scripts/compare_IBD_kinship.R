# Compare IBD segments with overall PLINK-inferred kinship values
# Author: SSG
#------------
# This script accepts two files (i.e. an IBD segment file and a PLINK-generated .genome file)
# It returns test statistics that reflect the concordance between IBD segments and genome-wide kinship values
# The required arguments are: Input (full file names)
# The optional arguments are: Genome length (number in bp, default 2.8e9), produce plot (T/F - default F), PI_HAT threshold (default 1/16)
# Ex.Rscript compare_IBD_kinship.R 

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

options(stringsAsFactors=F)
fileNameSplit <- strsplit(args[1], split=".", fixed=T)[[1]]
ext.ind <- which(fileNameSplit %in% c("match", "hom", "ibd"))
stem <- paste(fileNameSplit[-ext.ind], collapse=".")
ext <- fileNameSplit[ext.ind]

genome.length = 2.8e9
Plot = F
PI_HAT.thresh <- 1/16

# Update defaults, check for correct number of arguments
if (length(args)<2) {
  stop("At least two arguments must be supplied", call.=FALSE)
} else if (length(args)==3) {
  genome.length <- as.numeric(args[3])
} else if (length(args)==4) {
  genome.length <- as.numeric(args[3])
  Plot = args[4]
} else if (length(args)==5) {
  genome.length <- as.numeric(args[3])
  Plot = args[4]
  PI_HAT.thresh <- args[5]
} else if (length(args)>5) {
  stop("Too many arguments", call.=FALSE)
}

if (ext=="hom") {
  seg.file <- read.table(args[1], header=T)
  colnames(seg.file)[7:8] <- c("start", "end")
} else if (ext=="match") {
  seg.file <- read.table(args[1])
  if (ncol(seg.file)==13) {
    colnames(seg.file)[1:5] <- c("IID1","IID2","CHR", "start", "end")
  } else if (ncol(seg.file)==15) {
    colnames(seg.file)[c(2,4,5:7)] <- c("IID1","IID2","CHR", "start", "end")
  }
} else if (ext=="ibd") {
  seg.file <- read.table(args[1])
  colnames(seg.file)[c(1,3,5:7)] <- c("IID1","IID2","CHR", "start", "end")
}

genome.file <- read.table(args[2], header=T)
genome.file <- cbind(genome.file, as.data.frame(matrix(ncol=2, nrow=nrow(genome.file))))
  colnames(genome.file)[15:16] <- c("seg.kin.prop", "genome.kin.prop")

ID1 <- unique(genome.file$IID1)
ID2 <- unique(genome.file$IID2)

for (i in ID1) {
  for (j in ID2) {
    pair.segs <- seg.file[which(seg.file$IID1==i & seg.file$IID2==j),]
    total.IBD <- sum(as.numeric(pair.segs$end - pair.segs$start))
    # Calculate proportion of genome covered by IBD segments
    seg.kin.prop <- total.IBD/genome.length
    pair.genome <- genome.file[which(genome.file$IID1==i & genome.file$IID2==j),]
    # Calculate Z1 + 2 x Z2 from PLINK genome file
    genome.kin.prop <- pair.genome$Z1 + 2*pair.genome$Z2
    genome.file[which(genome.file$IID1==i & genome.file$IID2==j),c("seg.kin.prop", "genome.kin.prop")] <- c(seg.kin.prop, genome.kin.prop)
  }
}

if (Plot==T) {
  max.val <- max(genome.file[,c("seg.kin.prop", "genome.kin.prop")])
  plot(genome.file$genome.kin.prop, genome.file$seg.kin.prop, xlim=c(0, max.val), ylim=c(0, max.val), xlab="Z1 + (2 x Z2) from PLINK", ylab="Proportion of genome covered by IBD", main=paste(args[1], "-", args[2]))
    abline(a=0, b=1)
    abline(v=PI_HAT.thresh, lty=2, col="red")
}

# Calculate root mean squared error (RMSE) for pairs with a PI_HAT value over threshold
resids <- genome.file[genome.file$PI_HAT > PI_HAT.thresh, "genome.kin.prop"] - genome.file[genome.file$PI_HAT > PI_HAT.thresh, "seg.kin.prop"]
RMSE <- sqrt(mean(resids^2))
# Normalize RMSE by observed kinship values
NRMSE <- RMSE/(max(genome.file$genome.kin.prop)-min(genome.file$genome.kin.prop))

cat(NRMSE, "\n")
