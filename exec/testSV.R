# mkuhn, 2016-01-07
# test SVs
# DUPLICATION shows differences in train and test data for coverage


library(svmod)
library(GenomicRanges)

#setupParallel(cpus = 3L)


#
# patient 45
# pat00045  chr5:1277755-1279911        germline  (not well covered)
# pat00045  chr5:3599722-3601795        germline  (well covered 3599752-3601622)
# pat00045  chr5:35703655-35705764      DSV duplication

ngsProp <- ngs.prop()
bamPat45_WT <- file.path(getVirtualPatientPath(what = "mapping", .ngs.prop = ngsProp), "pat00045_N.s.bam")
bamPat45_MUT <- file.path(getVirtualPatientPath(what = "mapping", .ngs.prop = ngsProp), "pat00045_T.s.bam")
stopifnot( checkFile( bamPat45_MUT), checkFile(bamPat45_WT) )


# find covered regions in Pat45
covReg <- findCoveredRegionsInTarget(bamFile_WT = bamPat45_WT, bamFile_MUT = bamPat45_MUT, chrom = CHROM, extra.margin = 25L)
targetReg <- covReg[3L]


SVdup <- SV(type="dup", chrom=CHROM, posL = 3600000L, seqLen = 3) 
SVdup
getLocus(SVdup)

#grDup <- GenomicRanges::GRanges(seqnames = CHROM, ranges = IRanges::IRanges(start = 1277995, end = 1278005))
#grDup <- GenomicRanges::GRanges(seqnames = CHROM, ranges = IRanges::IRanges(start = 1277995, end = 1278055))
grDup <- GenomicRanges::GRanges(seqnames = CHROM, ranges = IRanges::IRanges(start = 3599995, end = 3600005))
width(grDup)
# duplication is added at breakpoint, ie. BP lies between the two copies
getContextSeq(sv = SVdup, gr = grDup)


# deletion
injectedBAMPat45 <- injectSVInMapping(bamFile = bamPat45_MUT, sv = SVdup, localMappingMargin = 1e5, ngsProp=ngsProp)

GenomicRanges::subsetByOverlaps(query = covReg, subject = grDup)
GenomicRanges::findOverlaps(query = covReg, subject = grDup)
assessClippedBasesInRegion(bamFile = bamPat45_MUT, chrom = CHROM, .targetStartPos = )




# big duplication
SVdup2 <- SV(type="dup", chrom=CHROM, posL = 3600000L, seqLen = 75) 
SVdup2
getLocus(SVdup2)

injectedBAMPat45_dup2 <- injectSVInMapping(bamFile = bamPat45_MUT, sv = SVdup2, localMappingMargin = 1e5, ngsProp = ngsProp)

clipInfo <- assessClippedBasesInRegion(bamFile = injectedBAMPat45_dup2[["mergedBam"]], chrom = chrom, do.plot=T,
                           minThreshold = 0.05,  windowSize=31L, minThresholdAbs=2L, thresholdQuantile=0.95,
                           .targetStartPos = start(targetReg), .targetEndPos = end(targetReg))

clipInfo.pos <- clipInfo[["peakInfo"]]
clipInfo.pos[c("pos.5p", "pos.3p")]
# continue here with checking coverage, insert size,  etc..
