#!/usr/bin/env Rscript
# 
# mkuhn, 20140718
# A script to explore patient mappings from a simulation run.
# It produces graphical illustration of clipped read pattern at simulated SV-positions


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

SV_TYPE <- "inversion" #  "insertion" # "deletion" #"duplication" #"no" #
NBR_REGIONS <- 5L
RSEED <- 123L

set.seed(RSEED)
stopifnot( NBR_REGIONS %% 5 == 0 ) #for easy plotting


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

# # MAPPING directory
mappingDir <- getVirtualPatientPath(baseDir, what = "mapping", .sample.prop=mySampleProp, .ngs.prop=myNGSProp)
loginfo("Using mapping directory: %s.", mappingDir)

# featureFilename <- paste0(mappingDir, "features.rds")
# featureDat0 <- NULL
# 
# if ( file.exists(featureFilename) ){
#   featureDat0 <- readRDS(featureFilename)
#   loginfo("Found %d existing features with %d columns.", NROW(featureDat0), NCOL(featureDat0))
# }


# sample some patients regions with the given SV_TYPE
patData.sv <- patData[patData$SVtype == SV_TYPE,]
patData.sv <- patData.sv[sample(NROW(patData.sv), size = NBR_REGIONS),]
patData.sv <- patData.sv[order(patData.sv$SVlength),]



opar <- par(mfrow=c(min(5, NBR_REGIONS), 1))
for (i in 1:NBR_REGIONS){

  mySVtype <- patData.sv[i, "SVtype"]
  mySV <- patData.sv[i, c("SVlength", "SVstart")]
  tumorBamFile <- paste0(mappingDir, patData.sv$patIdStr[i], "_T.bam")
  sv_cand.Regions <- exploreClipping_G(bamFile = tumorBamFile, startPos = patData.sv$startPos[i], endPos = patData.sv$endPos[i], SVtype = mySVtype, SVlength = mySV[1], SVstart = mySV[2],
                      windowSize = 20L) #, all.Sbases=FALSE, focus.SV = TRUE)
#   exploreClipping_BAM(bamFile = paste0(mappingDir, patData.sv$patIdStr[i], "_T.bam"), startPos = patData.sv$startPos[i], endPos = patData.sv$endPos[i], SVtype = mySVtype, SVlength = mySV[1], SVstart = mySV[2],
#                                   windowSize = 10L, all.Sbases=FALSE, add.Q1Pos=FALSE, focus.SV = TRUE)

}

par(opar)

