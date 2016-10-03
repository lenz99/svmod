#!/usr/bin/env Rscript
# mkuhn, 20140624
# mkuhn, 20140716: added optparse
# mkuhn, 20150320: added coverage and tumorLoad as parameters for the sampling
# mkuhn, 20151128: added flag useSimInfo if we rely on screening or not (i.e. use directly the knowledge what region has an SV and where exactly the simulated SV is)
#
# Example script how to extract features from simulated patient mappings.
# It uses the installed version of svmod package
#
# INput: simulated patient data (mapping (BAM file) for [N]- and [T]-samples corresponding to coverage-tumorload simulation)
# OUTput: feature data as RDS-object (while keeping old feature data). Learners can be trained on the features of the training data and evaluated on the test data (see benchmarkLearners.R)
#
# 


# for quick setup
if (FALSE){
  nbrCores <- 4L
  myNGSProp <- svmod::ngs.prop()
  mySampleProp <- svmod::sample.prop(cov = 60, tumorLoad = 90)
  
  USE_SIM_SV_INFO <- FALSE
  startID <- 1
  endID <- 1
  
  myMargin <- 50L
  myMaxNbrCovReg <- 6L
  rs <- 1234L
}

baseDir <- svmod::BASEDIR
myCHROM <- "chr5"


library(logging, quietly = TRUE); basicConfig()
logging::setLevel(level='INFO')  #set level of root logger
# logging::basicConfig()
# logging::setLevel('DEBUG', getHandler('basic.stdout')) # set level of handler

library(dplyr, quietly = TRUE)


TIMESTAMP <- format(Sys.time(), "%Y%m%d_%H%M") #format(Sys.Date(), "%Y%m%d")
addHandler(writeToFile, file=paste0(baseDir, "log/extractFeatures_", myCHROM, "_", TIMESTAMP, ".log"), level="WARN")
#level="DEBUG")


library(parallel, quietly = TRUE)
library(foreach, quietly = TRUE)
library(doParallel)
library(IRanges, quietly = TRUE)



# options and parameters -----

library(optparse, quietly = TRUE)


myOptions <- list(
  make_option(c("-v", "--devmode"), action="store_true", help="svmod is installed in dev-mode!", default=FALSE),
  make_option(c("-t", "--testrun"), action="store_true", help="Do only a test run. Do not save features to disk.", default=FALSE),
  make_option(c("-c", "--coverage"), type="integer", help="Coverage in sampling. [%default]", default=60L),
  make_option(c("-l", "--tumorLoad"), type="integer", help="Tumorload [%] in sampling. [%default]", default=90L),
  make_option(c("-i", "--useSimInfo"), action="store_true", help="Use coordinates of simulated SVs? If FALSE, rely on clipped base screening. Default is [%default]", default=FALSE),
  make_option(c("-m", "--margin"), type="integer", help="Use what margin left and right of simulated SV region for feature extraction [%default bp]", default=25L),
  make_option(c("-p", "--parallel"), type="integer", help="Number of Cores to use. 0: automatic, 1: no parallel. [%default]", default=2L),
  make_option(c("-s", "--startPatID"), type="integer", help="First Patient ID to do feature extraction [%default]", default=-1L),
  make_option(c("-e", "--endPatID"), type="integer", help="Last Patient ID to do feature extraction [%default]", default=-1L),
  make_option(c("-x", "--randomSeed"), type="integer", help="Set random seed [%default]", default=1234L),
  make_option(c("-r", "--maxCovRegion"), type="integer", help="Maximal number of covered regions per patient. 0 means no restriction. For now only in use when not using simulation. [%default]", default=0L)
)


myParser <- OptionParser(usage="%prog [options]", description="Extract features from simulated patients.", option_list=myOptions)
parsed_args <- parse_args(myParser, positional_arguments=TRUE)

given_opts <- parsed_args$options
given_args <- parsed_args$args


if (length(given_args) > 0L || any(nzchar(given_args)) ){
  logwarn("This script does not have any arguments, only options. Sorry.")
  print_help(myParser)
  q("no", status=1L, runLast=FALSE)
}


startID <- given_opts$startPatID
endID <- given_opts$endPatID

myMargin <- as.integer(ceiling(given_opts$margin))
stopifnot( is.integer(myMargin), myMargin > 0L )

nbrCores <- if ( given_opts$parallel <= 0L )
  max(1L, ceiling(sqrt(detectCores())), na.rm = TRUE) else given_opts$parallel

# how many covered regions per patient to consider at max for feature extraction?
myMaxNbrCovReg <- if (is.numeric(given_opts$maxCovRegion) && given_opts$maxCovRegion > 0L) given_opts$maxCovRegion else NULL






# setup -----

## LOAD Package: svmod
if (isTRUE(given_opts$devmode)){
  library(devtools, quietly = TRUE)
  dev_mode(on = TRUE)
}

library(svmod, quietly = FALSE)

#' Do we want to directly use the information on simulated SVs in regions?
#' TRUE to use info on SV presence and position. Allows to solely evaluate the learning algorithm (and not the filtering of covered/clipped region)
#' FALSE will use feature extraction from "realPatients.R" 
USE_SIM_SV_INFO <- given_opts$useSimInfo #TRUE


myNGSProp <- ngs.prop()
mySampleProp <- sample.prop(cov = given_opts$coverage, tumorLoad = given_opts$tumorLoad)



svmod::setupParallel(cpus = nbrCores)



# patient information: which regions were simulated, where is the DSV
patDataFile <- file.path(getVirtualPatientPath(baseDir, what = "patient"), "patData.rds")

patData <- readRDS(patDataFile)
stopifnot( NROW(patData) >= 3L )
#head(patData, n=2L)

patIds <- unique(patData$patId)
logging::loginfo("Found data in %s for %d patients.", patDataFile, length(patIds))




# first check endID, then startID
if ( is.numeric(endID) && endID >= 1L )  patIds <- patIds[which(patIds <= endID)]
if ( is.numeric(startID) && startID > -1L) patIds <- patIds[which(patIds >= startID)]


if ( length(patIds) == 0L ){
  logwarn("No patients for feature extraction after patID-filtering. Quitting.")
  q(save = "no", status = 1L, runLast = FALSE)
} else logging::loginfo("Feature extraction requested for %d patients %s", length(patIds), paste(patIds, collapse="+"))



# features for all candidate regions in all patients
logging::loginfo("Start to extract features from patient's mapping patterns with %d bp margin %s.",
        myMargin, if (USE_SIM_SV_INFO) "around simulated SV region" else "after screening round")


# mapping directory for given setting of coverage and TL
mappingDir <- getVirtualPatientPath(baseDir, what = "mapping", .sample.prop=mySampleProp, .ngs.prop=myNGSProp)


# mkuhn, 2016-02-05
rs <- given_opts$randomSeed
if ( is.numeric(rs) && rs > 0L ){
  logging::loginfo("Set random seed to %d.", ceiling(rs))
  set.seed(seed = ceiling(rs))
}


# mkuhn, 2016-02-08: SV-length dist for SV-injection
# should best match the setting of the simulation..
mySVmean <- 175
mySVsd <- 15



# build feature dataframe (train and test data) -----

featureDat <- if ( isTRUE(USE_SIM_SV_INFO) ){
  # getFeaturesFromMappings was  parallel in the regions, now I swapped that:
  # parallel at whole patient level. This seems more thread-safe.
  # I use simulation data of patient, i.e. the known coordinates of simulated differential SVs
  #featureDat <- foreach(patId=patIds, .combine = 'rbind', .multicombine = FALSE, .errorhandling = "stop") %dopar% { ##} %do% {  #.maxcombine = 4L
  foreach::foreach(myPatId=patIds, .combine = 'rbind', .multicombine = FALSE, .errorhandling = "stop") %dopar% {  
    
    logging::loginfo("Extract features for patient %s using simulated SV information.", myPatId)
    getVirtualPatientFeatures(patId = myPatId, patData = patData, baseDir=baseDir, sample.prop=mySampleProp, ngs.prop=myNGSProp, margin.bp=myMargin)
    
  }
} else {
  
  #featureList <- foreach::foreach(myPatId=patIds, .combine = 'list', .multicombine = TRUE, .errorhandling = "stop") %dopar% {  
  #featureList <- foreach::foreach(myPatId=patIds, .combine = 'list', .multicombine = TRUE, .errorhandling = "pass") %dopar% {    
  # TESTING: single patient
  myPatId <- patIds[1]
  logging::logwarn("Use only first patient %s", myPatId)
    
    logging::loginfo("Extract features for patient %s **without** using the simulated SV information.", myPatId)
    
    pat.n.bamFile <- list.files(path=mappingDir, pattern=paste0("^", getPatStr(myPatId),"[_]N.s.bam$"), full.names = TRUE)
    pat.t.bamFile <- list.files(path=mappingDir, pattern=paste0("^", getPatStr(myPatId),"[_]T.s.bam$"), full.names = TRUE)
    
    stopifnot( length(pat.n.bamFile) == 1L, length(pat.t.bamFile) == 1L )
    stopifnot( checkFile(pat.n.bamFile), checkFile(pat.t.bamFile) )
    featDat.pat <- getFeatureDataFromPatient(patId = myPatId, bamFile_WT = pat.n.bamFile, bamFile_MUT = pat.t.bamFile,
                                             chrom = myCHROM, localMappingMargin=1e5, propPos.train=.5, margin.bp = myMargin,
                                             SVlength_mean = mySVmean, SVlength_sd = mySVsd,
                                             maxNbrCoveredRegions = myMaxNbrCovReg, useIntermediate = TRUE, saveIntermediate = TRUE)
    
    if (is.null(featDat.pat)){
      logging::logwarn("No feature data for patient %d", myPatId)
    }
    
    featureList <- featDat.pat
    
    #featDat.pat
  ##}#foreach featureList 
  
  ## DEbugging
  saveRDS(featureList, file = paste0("~/featureList_P", paste0(patIds, collapse = "_"), toString(mySampleProp), ".rds"))
  ## DEBUGGING
  
  # feature data frame of training and test data
  #features.train <- do.call(rbind, lapply(featureList, function(l) l[["train"]]))
  #features.test0 <- do.call(rbind, lapply(featureList, function(l) l[["test"]]))
  # DEBUGGING for single patient!!
  features.train <- featureList[["train"]]
  features.test0 <- featureList[["test"]]
    
  logging::loginfo("Features of training data has first columns: %s", paste(head(names(features.train)), collapse=" - "))
  
  # myMargin is re-used here for estimating the status on test data
  features.test <- addStatusToTestFromSim(features.test0, patData, myCHROM = myCHROM, evalMargin = myMargin, regionIsMatch=FALSE, featNames = names(features.train))
  
  
  rbind(features.train, features.test)
  
  
}#esle



logging::loginfo("featureDat is a %s of length %d with %d rows.", class(featureDat), length(featureDat), NROW(featureDat))
# combine to a dataframe
#featureDat <- do.call(rbind, featureDat)


if (isTRUE(given_opts$testrun)){
  
  head(featureDat, n=2L)
  loginfo("TEST RUN: NOT writing out %d features entries.", NROW(featureDat))
  q(save = "no", status = 0L)
}

## read in existing FEATURE data from disk
stopifnot( ! isTRUE(given_opts$testrun) )




# write out -----

# FEATURE directory
featureDir <- getVirtualPatientPath(baseDir, what = "feature", .sample.prop=mySampleProp, .ngs.prop=myNGSProp)
logdebug("Using feature directory: %s.", featureDir)
try(dir.create(featureDir, showWarnings = FALSE, recursive = TRUE), silent = TRUE)


featBaseName <- if (isTRUE(USE_SIM_SV_INFO)) "features_sim_noScreen_margin" else "features_sim_screen_margin"
featureFilename <- file.path(featureDir, paste0(featBaseName, myMargin, ".rds"))

# #####
# # mkuhn, 2015-12-24: TESTING
# saveRDS(featureDat, file=featureFilename)
# 
# cat("\nExiting.. TESTING\n")
# q(save = "no", status=0L)




newPatients <- unique(featureDat[, "patId"])
loginfo("Created in total %d new feature rows with %d feature columns from %d new patients.",
        NROW(featureDat), NCOL(featureDat), length(newPatients))






# existing data?
if ( file.exists(featureFilename) ){
  featureDat0 <- readRDS(featureFilename)
  logging::loginfo("Found %d existing features with %d columns.", NROW(featureDat0), NCOL(featureDat0))
  
  
  # check if old feature dataframe has same format than new features
  if ( (! is.null(featureDat0)) && NROW(featureDat0) > 0L && identical(names(featureDat0), names(featureDat)) ){
    oldFeatures <- featureDat0[! featureDat0[, "patId"] %in% newPatients, ]
    
    if ( NROW(oldFeatures)>0L ){
      loginfo("Keeping %d old feature rows from previous feature extraction runs.", NROW(oldFeatures))
      featureDat <- rbind(oldFeatures, featureDat)
    }#fi
  }
}#fi featureDat0


# both data sets should fit together
if ( NROW(featureDat) != NROW(patData) && startID == -1L && endID == -1L ){
  logging::logwarn("Resulting feature data set does not have same number of rows as patient data set. Proceeding anyway. Dropped loci or patID-filter?")
}

saveRDS(featureDat, file=featureFilename)
loginfo("Feature data saved to disk in ::%s::", featureFilename)

