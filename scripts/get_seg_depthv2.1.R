# Get segment depth v2.1
# Author: SSG
# v2 and v2.1 edits by AWR
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
  if (!require(x,character.only = TRUE)){
    install.packages(x,dep=TRUE)
    if(!require(x,character.only = TRUE)) stop("Package not found")
  }
}

# Call packages
pkgTest("data.table")

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
colnames(bim)[1:4] <- c("chr","rsid","dist","pos")
bim[, pos2 := pos ]

#set key columns for data.tables
setkey(seg.file, chr, segment.start, segment.end)
setkey(bim, chr, pos, pos2)

#calculate depth counts
match_list <-list() #create empty list to fill with chr data.frames
for (i in 1:22){
  seg.chr <- seg.file[which(seg.file$chr == i),] #subset seg.file by chromosome
  bim.chr <- bim[which(bim$chr == i),] #subset bim file by chromosome
  match_list[[i]] <- foverlaps(bim.chr, seg.chr) [, .(count = sum(!is.na(segment.start))), by = .(chr,pos, pos2) ][, pos2 := NULL ] #find overlaps between segs and bim markers and count number of times the bim markers are coverged by the segments
  invisible(gc()) #this is a memory intensive process so clean up mem after each overlap count
}

newdepthcounts <- do.call(rbind,match_list) #combine chr counts into single dataframe
colnames(newdepthcounts)<-c("CHR","BP","COUNT") #assign col names

write.table(newdepthcounts, "", quote=F, row.names=F, col.names=F, sep=" ") #write to out
