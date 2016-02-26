
library(svmod)
# Pindel finds ITD at enrichment target  chr13:28607945-28608622
regFLT3 <- GenomicRanges::GRanges(seqnames = "chr13",
                                  #ranges = IRanges::IRanges(start = 28575411, end=28676729)) #FLT3 region
                                  ranges = IRanges::IRanges(start = 28607945, end=28608622))

REAL_DIR <- file.path(BASEDIR, "biotec")
INTERMEDIATE_DIR <- file.path(REAL_DIR, "intermediate")

patIds <- c("P27", "P39", "P40", "P136")
for (patId in patIds){
  
  bamFile_MUT <- list.files(path = REAL_DIR, pattern = paste0("^", patId, "_T.+[.]s[.]dedup.bam$"), full.names = TRUE)
  stopifnot( length(bamFile_MUT) == 1L )
  coveredRegions <- readRDS(file.path(INTERMEDIATE_DIR, paste0(patId, "_chr13_covReg.rds")))
  mapInfo <- readRDS(file.path(INTERMEDIATE_DIR, paste0(patId, "_mappingInfo.rds")))
  
  myMinThreshold <- mapInfo[["clippedBases"]]$mean + 3 * mapInfo[["clippedBases"]]$sd
  cat(sprintf("Patient %s: clipped mean %f with SD %f and peak minThreshold %f and myMinThreshold %f.\n",
              patId, mapInfo[["clippedBases"]]$mean, mapInfo[["clippedBases"]]$sd, mapInfo[["clippedBases"]]$peakMinThreshold, myMinThreshold))
  
  coveredRegionsFLT3 <- coveredRegions[GenomicAlignments::findOverlaps(query = regFLT3, subject = coveredRegions)@subjectHits]
  
  cat(sprintf("Patient %s: Covered regions at FLT3-ITD %d\n", patId, length(coveredRegionsFLT3)))
  if (length(coveredRegionsFLT3) == 1L){
    
    clipPeaks <- assessClippedBasesInRegion(bamFile = bamFile_MUT, chrom = as.character(GenomicAlignments::seqnames(regFLT3)),
                                            minThreshold = mapInfo[["clippedBases"]]$peakMinThreshold, do.plot=TRUE,
                                            .targetStartPos = IRanges::start(coveredRegionsFLT3), .targetEndPos = IRanges::end(coveredRegionsFLT3)) #[["clippedInfo"]]
    title(sprintf("Patient %s, Region %s\nPeaks %d", patId, toString(regFLT3), clipPeaks$peakInfo[1]))
  }
  
}



# regionStarts <- IRanges::start(coveredRegionsFLT3)
# regionEnds <- IRanges::end(coveredRegionsFLT3)
# 
# 
# for (i in 1:length(coveredRegionsFLT3)){
#   regionStart <- regionStarts[i]; regionEnd <- regionEnds[i]
#   assessClippedBasesInRegion(bamFile = bamFile_MUT, chrom = chrom,
#                              minThreshold = clippedPropThreshold, do.plot=TRUE,
#                              .targetStartPos = regionStart, .targetEndPos = regionEnd) #[["clippedInfo"]]
#   title(paste("Region ",i))
# }
