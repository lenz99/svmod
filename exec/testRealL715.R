# cf realPatients.R

library(svmod)

#BASEDIR <- "/Users/kuhnmat/SVdata/"
chrom <- "chr13" #"chr5"

library(dplyr)
library(foreach)

bamFile_MUT <- file.path(BASEDIR, "biotec", "L709.tr.s.dedup.bam")
bamFile_WT  <- file.path(BASEDIR, "biotec", "L708.tr.s.dedup.bam")

stopifnot( file.exists(bamFile_WT), file.exists(bamFile_MUT) )
myCovRegP27 <- findCoveredRegionsInTarget(bamFile_WT = bamFile_WT, bamFile_MUT = bamFile_MUT, chrom = chrom,
                                              minCov_WT = 30, minAvgCov_MUT=10)


# covering regions of patients
# alternative: use subsetByOverlaps..
# regNPM1 <- GenomicRanges::GRanges(seqnames = "chr5",
#                                   ranges = IRanges(start = 170812708, end=170839888))
# stopifnot( chrom == 'chr5' )

# FLT3 in hg19
# ITDseek focusses on chr13:28607161-28609590
# Pindel finds ITD at enrichment target  chr13:28607945-28608622
regFLT3 <- GenomicRanges::GRanges(seqnames = "chr13",
                                  ranges = IRanges::IRanges(start = 28575411, end=28676729))
stopifnot( chrom == 'chr13' )


coveredRegions <- myCovRegP27[GenomicAlignments::findOverlaps(query = regFLT3,
                                               subject = myCovRegP27)@subjectHits]


margin.bp=ceiling(READL/4); maxNbrCoveredRegions=100L; propPos.train=.25


patMappingInfo <- getPatientMappingInfo(bamFile_WT = bamFile_WT, bamFile_MUT = bamFile_MUT, chrom = chrom)

# from realpatients..
clippedRes_MUT <- foreach::foreach(i=seq_along(coveredRegions), .combine='rbind') %do% {
  #cat("\n ** i=",i, "**\n")
  # assess clipped-base peaks in covered region: either 0 or 1 or 2 entries (3p-5p peak pair) that correspond to highest potential SV
  regionStart <- start(coveredRegions)[i]
  regionEnd <- end(coveredRegions)[i]
  # findClippedPeaks_G_OLD
  peaks_MUT <- assessClippedBasesInRegion(bamFile = bamFile_MUT, chrom = chrom, minThreshold = 0.05, minThresholdAbs=3L, do.plot=TRUE,
                                          .targetStartPos = regionStart, .targetEndPos = regionEnd)
  cInfo <- peaks_MUT[["clippedInfo"]]
  #pInfo <- peaks_MUT[["peakInfo"]]
  #negCandPos <- peaks_MUT[["nonPeakPos"]]
  cInfo
}#foreach clippedRes


clippedRes_MUT


targetBedFile = file.path(BASEDIR, "refGen", "TruSeq_exome_targeted_regions.hg19.bed")


# feature extraction for a positive case (test data)
# example region w/ index i
i <- 3L
regStart <- clippedRes_MUT[i, "regionStart"]
regEnd <- clippedRes_MUT[i, "regionEnd"]
mySVCandStart <- clippedRes_MUT[i, "pos.left"]
mySVCandEnd   <- clippedRes_MUT[i, "pos.right"]
mySVClipInfo <- list(peakInfo=clippedRes_MUT[i, c("nbrPeaks", "posDiff")])

cbind(
  data.frame(patId=patId, chrom=chrom), SVmargin=margin.bp,
  t(getFeaturesFromMapping(patId=patId, bamFile_WT = bamFile_WT, bamFile_MUT = bamFile_MUT, patMappingInfo = patMappingInfo, 
                           chrom = chrom, startPos = mySVCandStart - margin.bp, endPos = mySVCandEnd + margin.bp,
                           targetStartPos = regStart, targetEndPos = regEnd, clipInfo = mySVClipInfo)[-1L])
)


# training!
# negative case
j <- 6
mySV.types <- sample(SV_TYPES, size = length(negCandidateRegion.ind), replace = TRUE)
mySV.lengths <- abs(ceiling(rnorm(length(negCandidateRegion.ind), mean=SV_LENGTH, sd=9L)))+1L
