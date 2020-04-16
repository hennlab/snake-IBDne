# Get segment depth v3
# Author: SSG
# v2-v3 edits by AWR
#------------
# This script accepts a file listing segments (i.e. PLINK-generated RoH file, GERMLINE2 generated .match file, IBDSeq/Beagle generated .ibd file) and a .bim file
# It returns the depth of RoH at every SNP position in the .bim file
# The required arguments are: Input (full file names)
# Ex.Rscript get_RoH_depth.R CHABU.Ethiopian_MEGA_RoH_segsJoined.hom CHABU.Ethiopian_MEGA_HWE-filtered.bim
#------------
# v2 uses functions from data.table to make calculating counts more efficient
# v2.1 splits the seg.file and bim file by chromosome to get around an overflow issue
# when calculating counts on large datasets, as well as
# reading in files using fread to make this step more efficient
#------------
# v3 takes a new approach to conserve memory
# using the GenomicRanges library to count the overlap segments
#------------

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
# Check for correct number of arguments
if (length(args)<2) {
  stop("At least two arguments must be supplied", call.=FALSE)
} else if (length(args)>2) {
  stop("Too many arguments", call.=FALSE)
}

# Define functions
pkgTest <- function(x){
  if (!require(x,character.only = TRUE, quietly = T)){
    install.packages(x,dep=TRUE)
    if(!require(x,character.only = TRUE, quietly = T)) stop("Package not found")
  }
}

# Call packages
pkgTest("data.table")
pkgTest("GenomicRanges")

# Read in data
options(stringsAsFactors=F)

fileName <- unlist(strsplit(args[1], split=".", fixed=TRUE))
ext <- fileName[length(fileName)]
if (ext == "hom") {
  seg.file <- fread(args[1], header=T)
  colnames(seg.file)[4]<-"chr"
  colnames(seg.file)[7:8]<-c("segment.start", "segment.end")
} else if (ext == "ibd") {
  seg.file <- fread(args[1], header=F)
  colnames(seg.file)[5:7] <- c("chr", "segment.start", "segment.end")
} else if (ext == "match") {
  seg.file <- fread(args[1], header=F)
  if (ncol(seg.file==13)) {
    colnames(seg.file)[3:5] <- c("chr", "segment.start", "segment.end")
  } else if (ncol(seg.file==15)) {
    colnames(seg.file)[5:7] <- c("chr", "segment.start", "segment.end")
  }
}

bim <- fread(args[2])
colnames(bim)[1:4] <- c("chr","rsid","dist","start")
bim[, end := start ]

#set key columns for data.tables
setkey(seg.file, chr, segment.start, segment.end)
setkey(bim, chr, start, end)

#calculate depth counts
#assign data.tables to GRanges objects
GRseg <- makeGRangesFromDataFrame(seg.file) 
GRbim <- makeGRangesFromDataFrame(bim)
#make new dataframe combining bim site info with counts
newdepthcounts <- as.data.frame(bim[,c("chr","start")])
#count number of segments that overlap sites in bim file assign to new column
newdepthcounts$count <- countOverlaps(GRbim, GRseg, type = "any") 
colnames(newdepthcounts)<-c("CHR","BP","COUNT") #assign final column names


write.table(newdepthcounts, "", quote=F, row.names=F, col.names=F, sep=" ") #write to out
