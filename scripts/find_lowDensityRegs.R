# AUTHOR: SSG
# This script identifies and outputs regions of low SNP density
# It requires a .bim file from which to estimate SNP density
# The required arguments are: Bim file (the full name the file)
# 		   	    : Output prefix (minus file extension)
#                           : The window size (1 Mb recommended) 
# The optional arguments are: Percentage threshold - what cutoff to use to define low density regions? Default is the top 5% lowest density regions on the chromsome or 1/50000 * window size (whichever is higher)
# USAGE:
## Rscript find_lowDensityRegs.R {input.bim} {output.prefix} 1e6 0.05

##!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)


# Defaults
perc.thresh <- 0.05

# Error handling and updating default values
if (length(args)<3) {
  stop("At least three arguments must be supplied", call.=FALSE)
} else if (length(args)==4) {
  perc.overlap <- as.numeric(args[4])
} else if (length(args)>4) {
  stop("Too many arguments", call.=FALSE)
}

options(stringsAsFactors=F)
Bim <- read.table(args[1])
WinSize <- as.numeric(args[3])

fileNameSplit <- strsplit(args[1], split=".", fixed=T)[[1]]
ext.ind <- match("bim", fileNameSplit)
stem <- paste(fileNameSplit[-ext.ind], collapse=".")

genome.LDRs <- as.data.frame(matrix(ncol=3, nrow=0))
for (chr in unique(Bim$V1)) {
  # Calculate SNP density across the chromosome in windows
  chr.data <- Bim[Bim$V1==chr,]
  start <- chr.data[1, "V4"]
  end <- chr.data[nrow(chr.data), "V4"]
  Hist <- hist(chr.data$V4, breaks=seq(start, end+WinSize, WinSize), plot=F)
  Counts <- Hist$counts
  Breaks <- Hist$breaks
  # Store count information per region
  chrWinInfo <- as.data.frame(matrix(ncol=4, nrow=length(Counts)))
  colnames(chrWinInfo) <- c("chr", "start", "end", "nSNPs")
  chrWinInfo$chr <- chr
  chrWinInfo$start <- Breaks[1:length(Counts)]
  chrWinInfo$end <- Breaks[2:length(Breaks)]
  chrWinInfo$nSNPs <- Counts
  # Find low density regions
  min.SNPs <- max(sort(Counts)[length(Counts)*perc.thresh], 1/5e4*WinSize)
  low.dens.wins <- chrWinInfo[which(chrWinInfo$nSNPs < min.SNPs), 1:3]
  # If low density regions are found join continguous windows into regions
  if (nrow(low.dens.wins) > 0) {
    Low.Density.Regions <- low.dens.wins[1, ]
    if (nrow(low.dens.wins) > 1) {
      for (x in 2:nrow(low.dens.wins)) {
        if (low.dens.wins[x, "start"] == Low.Density.Regions[nrow(Low.Density.Regions), "end"]) {
          Low.Density.Regions[nrow(Low.Density.Regions), "end"] <- low.dens.wins[x, "end"]
        } else {
          Low.Density.Regions[nrow(Low.Density.Regions)+1, ] <- low.dens.wins[x, ]
        }
      }
    }
  }
  genome.LDRs <- rbind(genome.LDRs, Low.Density.Regions)
}

genome.LDRs <- cbind(genome.LDRs, "low_density_region")

out_pre=args[2]

write.table(genome.LDRs, paste(out_pre, "_lowDensityRegions.txt", sep=""), quote=F, row.names=F, col.names=F, sep="\t")
