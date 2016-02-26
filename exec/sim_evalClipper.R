#!/usr/bin/env Rscript
# 
# mkuhn, 20150128
# A script to evaluate the clipped read finder on simulated data
######

DEV_MODE <- TRUE

library(logging); basicConfig()



if (isTRUE(DEV_MODE)){
  if (require(devtools))
    dev_mode(on = TRUE)
  else
    logwarn("Package devtools is missing. Hence, devmode is not used.")
}

## LOAD Package: SVMOD
svmod_pck <- require(svmod)
if (! svmod_pck){
  logerror("Package 'svmod' is missing.")
  q(save="no", status = 1)
}
stopifnot( svmod_pck )



###
### SETTINGS
###

RSEED <- 123L
set.seed(RSEED)


baseDir <- svmod::BASEDIR
myNGSProp <- ngs.prop()
mySampleProp <- sample.prop()




# patient information: which regions were simulated, where is the DSV
patDataFile <- paste0(baseDir, "patient/patData.rds")

patData <- readRDS(patDataFile)
stopifnot( NROW(patData) >= 3 )
#head(patData, n=2)
loginfo("Read in data from %d patients from file %s.", NROW(patData), patDataFile)

# mkuhn, 20141217: focus on new set of patients with rather small SVs
patData <- dplyr::filter_(patData, ~patId >= 1000)
patData <- dplyr::mutate(patData, SVstartPos=targetStartPos + SVstart, SVendPos=SVstartPos + SVlength)

# MAPPING directory
mappingDir <- getVirtualPatientPath(baseDir, what = "mapping", .sample.prop=mySampleProp, .ngs.prop=myNGSProp)
loginfo("Using mapping directory: %s.", mappingDir)




SV_TYPE <- "insertion"
TOLERANCE <- 30 # bp

clipperEval <- dplyr::select(patData, patId:SVstart, SVstartPos, SVendPos)
clipperEval$hit <- FALSE
clipperEval$nCandRegions <- 0

#ZZZ this is a quick and dirty evaluation of clipper
for (i in 1:100){
  tumorBamFile <- paste0(mappingDir, patData$patIdStr[i], "_T.bam")
  mySVtype <- patData[i, "SVtype"]
  mySV <- patData[i, c("SVlength", "SVstart")]
  clipperOut <- exploreClipping_G(bamFile = tumorBamFile, startPos = patData$startPos[i], endPos = patData$endPos[i], SVtype = mySVtype, SVlength = mySV[1], SVstart = mySV[2],
                                  windowSize = 20L, showPeak = TRUE) #, all.Sbases=FALSE, focus.SV = TRUE)
  
  # 
  leftMost <- pmin(clipperOut$pos.5p, clipperOut$pos.3p) - TOLERANCE
  rightMost <- pmax(clipperOut$pos.5p, clipperOut$pos.3p) + TOLERANCE
  
  clipperEval[i, "nCandRegions"] <- length(leftMost)
  # hit?
  clipperEval[i, "hit"] <- any(leftMost <= patData[i, "SVstartPos"] & rightMost >= patData[i, "SVstartPos"])
  
}
