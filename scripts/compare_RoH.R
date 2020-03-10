# Compare PLINK, GERMLINE2, IBDSeq RoH distributions
# Author: SSG
#------------
# This script accepts two RoH files (i.e. GERMLINE-generated .match file, IBDSeq-generated .hbd file or PLINK-generated .hom file)
# It returns test statistics that reflect the similarity of the second RoH distribution to the first
# The required arguments are: Input (full file names)
# The optional arguments are: Minimum length (in kb) of segments that should be considered (default 1000)
#                           : Plot - TRUE or FALSE (default FALSE)
# Ex.Rscript compare_RoH.R Chabu_Ethiopian_MEGA_RoH.hom CHABU_RoH_test_results.match 3000 F

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# Defaults
min.length = 1000
Plot = F

# Update defaults, check for correct number of arguments
if (length(args)<2) {
  stop("At least two arguments must be supplied", call.=FALSE)
} else if (length(args)==3) {
  min.length = as.numeric(args[3])
} else if (length(args)==4) {
  min.length = as.numeric(args[3])
  Plot <- as.logical(args[4])
} else if (length(args)>4) {
  stop("Too many arguments", call.=FALSE)
}

# Figure out what type of files are being given
arg1.split <- unlist(strsplit(args[1], split=".", fixed=T))
    arg1.ext <- arg1.split[length(arg1.split)]
arg2.split <- unlist(strsplit(args[2], split=".", fixed=T))
    arg2.ext <- arg2.split[length(arg2.split)]

# Read in these files
if (arg1.ext=="match") {
  Dist.1 <- read.table(args[1], stringsAsFactors=F)
    if (ncol(Dist.1)==13) {
      colnames(Dist.1)[2:5] <- c("IID", "CHR", "POS1", "POS2")
      Dist.1 <- Dist.1[Dist.1$V1==Dist.1$IID,]
    } else if (ncol(Dist.1)==15) {
      colnames(Dist.1)[c(2,5:7)] <- c("IID", "CHR", "POS1", "POS2")
      Dist.1 <- Dist.1[Dist.1$V3==Dist.1$IID,]
    }
  } else if (arg1.ext=="hom") {
    Dist.1 <- read.table(args[1], stringsAsFactors=F, header=T)
  } else if (arg1.ext=="hbd" | arg1.ext=="ibd"){
    Dist.1 <- read.table(args[1], stringsAsFactors=F, header=F)
    colnames(Dist.1)[c(1,5:7)] <- c("IID", "CHR", "POS1", "POS2")
    Dist.1 <- Dist.1[Dist.1$V3==Dist.1$IID,]
  }

if (arg2.ext=="match") {
  Dist.2 <- read.table(args[2], stringsAsFactors=F)
  if (ncol(Dist.2)==13) {
    colnames(Dist.2)[2:5] <- c("IID", "CHR", "POS1", "POS2")
    Dist.2 <- Dist.2[Dist.2$V1==Dist.2$IID,]
  } else if (ncol(Dist.2)==15) {
    colnames(Dist.2)[c(2,5:7)] <- c("IID", "CHR", "POS1", "POS2")
    Dist.2 <- Dist.2[Dist.2$V3==Dist.2$IID,]
  }
} else if (arg2.ext=="hom") {
    Dist.2 <- read.table(args[2], stringsAsFactors=F, header=T)
} else if (arg2.ext=="hbd" | arg2.ext=="ibd"){
    Dist.2 <- read.table(args[2], stringsAsFactors=F, header=F)
    colnames(Dist.2)[c(1,5:7)] <- c("IID", "CHR", "POS1", "POS2")
    Dist.2 <- Dist.2[Dist.2$V3==Dist.2$IID,]
}

for (POS in c("POS1", "POS2")) {
  Dist.1[,POS] <- as.numeric(Dist.1[,POS])
  Dist.2[,POS] <- as.numeric(Dist.2[,POS])
}

# Restrict analysis to only segments exceeding the defined threshold
Dist.1 <- Dist.1[which((Dist.1$POS2-Dist.1$POS1)>(min.length*1000)),]
Dist.2 <- Dist.2[which((Dist.2$POS2-Dist.2$POS1)>(min.length*1000)),]

# Match IDs in both distributions
if (TRUE %in% is.na(match(Dist.1$IID, Dist.2$IID))) {
  if (sum(na.omit(pmatch(Dist.1$IID, Dist.2$IID)))!=0) {
    for (i in unique(Dist.1$IID)) {
      Dist.2[grep(i, Dist.2$IID),"IID"] <- i
    }
  } else if (sum(na.omit(pmatch(Dist.2$IID, Dist.1$IID)))!=0) {
    for (i in unique(Dist.2$IID)) {
      Dist.1[grep(i, Dist.1$IID),"IID"] <- i
    }
  }
}
IDs <- unique(Dist.1$IID)

total.RoH <- sum(as.numeric(Dist.1$POS2 - Dist.1$POS1)) # calculate total RoH in dataset
total.diff <- 0 # start counters
total.match <- 0

if (Plot==T) {
  Dist.2.forPlotting <- Dist.2
}

suppressMessages(library("GenomicRanges"))
for (ind in IDs) {
  Dist.1.ind <- Dist.1[which(Dist.1$IID==ind),]
  Dist.1.ind.GRange <- makeGRangesFromDataFrame(Dist.1.ind[,c("IID","CHR","POS1","POS2")], seqnames.field="CHR", start.field="POS1", end.field="POS2")
  Dist.2.ind <- Dist.2[which(Dist.2$IID==ind),]
  if (nrow(Dist.2.ind)>0) {
    Dist.2.ind.GRange <- makeGRangesFromDataFrame(Dist.2.ind[,c("IID","CHR","POS1","POS2")], seqnames.field="CHR", start.field="POS1", end.field="POS2")  # Find overlaps between Dist.1 and Dist.2
    Hits <- suppressWarnings(findOverlaps(Dist.1.ind.GRange, Dist.2.ind.GRange)) # find all overlaps
    if (length(Hits)>0) {
      total.diff <- total.diff + sum(abs(as.numeric(Dist.1.ind[queryHits(Hits),"POS1"]-Dist.2.ind[subjectHits(Hits),"POS1"]))) + sum(abs(as.numeric(Dist.1.ind[queryHits(Hits),"POS2"]-Dist.2.ind[subjectHits(Hits),"POS2"])))
      Overlaps <- suppressWarnings(pintersect(Dist.1.ind.GRange[queryHits(Hits)], Dist.2.ind.GRange[subjectHits(Hits)]))
      percentOverlap <- width(Overlaps) / width(Dist.1.ind.GRange[queryHits(Hits)])
      total.match <- total.match + sum(as.numeric(apply(cbind(percentOverlap, Dist.1.ind[queryHits(Hits),"POS2"]-Dist.1.ind[queryHits(Hits),"POS1"]), MARGIN=1, prod)))
      Dist.2 <- Dist.2[-match(rownames(Dist.2.ind[unique(subjectHits(Hits)),]), rownames(Dist.2)),]
    } else {
      # add unmatched segments from first distribution to the total difference
      total.diff <- total.diff + sum(as.numeric(abs(Dist.1.ind$POS2-Dist.1.ind$POS1)))
    }
  } else {
    total.diff <- total.diff + sum(as.numeric(abs(Dist.1.ind$POS2-Dist.1.ind$POS1)))
  }
}

total.diff <- total.diff + sum(as.numeric(Dist.2[,"POS2"]-Dist.2[,"POS1"])) # add all unmatched segments from second distribution
perc.diff <- total.diff/total.RoH
perc.match <- total.match/total.RoH

cat(perc.diff, perc.match, "\n")

if (Plot==T) {
  Dist.2 <- Dist.2.forPlotting[which((Dist.2.forPlotting$POS2-Dist.2.forPlotting$POS1)>(min.length*1000)),]
  par(mfrow=c(2,4))
  # par(mfrow=c(4,6))
  chrs <- sort(unique(c(Dist.1$CHR, Dist.2$CHR)))
  IDs <- unique(Dist.1$IID)
  Ylim <- c(1,length(IDs)*2+1)
  for (chr in chrs) {
    Dist.1.chr <- Dist.1[which(Dist.1$CHR==chr),]
    Dist.2.chr <- Dist.2[which(Dist.2$CHR==chr),]
    Xlim <- range(c(Dist.1.chr$POS1, Dist.1.chr$POS2, Dist.2.chr$POS1, Dist.2.chr$POS2))
    plot(seq(Xlim[1],Xlim[2],1e7), type="n", xlim=Xlim, ylim=Ylim, ann=F, yaxt="n")
    title(main=paste("chr",chr))
    axis(2, at=seq(Ylim[1]+0.5, Ylim[2]-0.5, 2), labels=IDs, las=1, cex.axis=0.7)
    for (ind in IDs) {
      Dist.1.chr.ind <- Dist.1.chr[which(Dist.1.chr$IID==ind),]
      Dist.2.chr.ind <- Dist.2.chr[which(Dist.2.chr$IID==ind),]
      indPos <- match(ind, IDs)
      if (nrow(Dist.1.chr.ind)>0) {
        for (n in 1:nrow(Dist.1.chr.ind)) {
          segments(x0=Dist.1.chr.ind[n,"POS1"], x1=Dist.1.chr.ind[n,"POS2"], y0=indPos*2, col="red")
        }
      }
      if (nrow(Dist.2.chr.ind)>0) {
        for (n in 1:nrow(Dist.2.chr.ind)) {
          segments(x0=Dist.2.chr.ind[n,"POS1"], x1=Dist.2.chr.ind[n,"POS2"], y0=indPos*2-1, col="blue")
        }
      }
    }
    #chr.regs <- exclude.regions[exclude.regions$V1==chr,]
    #for (i in 1:nrow(chr.regs)) {
    #  if (chr.regs[i,4]=="depth_outlier") {
    #    polygon(x=as.numeric(chr.regs[i,c(2:3,3:2)]), y=c(0,0,90,90), density=NA, col=rgb(.25,.25,.85,.2))
    #  } else {
    #    polygon(x=as.numeric(chr.regs[i,c(2:3,3:2)]), y=c(0,0,90,90), density=NA, col=rgb(.75,.75,.75,.2))
    #  }
    #}
  }
}
