#!/usr/bin/env Rscript
# mkuhn, 20160127: based on sim_extractFeatures.R script
#
# Example script how to extract features from real patient mappings.
# It uses the installed version of svmod package
#
# INput: mapping (BAM file) for [N]- and [T]-samples of real patients
# OUTput: feature data as RDS-object


# for quick setup
if (FALSE){
  patIdentifier <- "P27"
  nbrCores <- 3L
  myMargin <- 25L
  myMaxNbrCovReg <- 3L
}


baseDir <- svmod::BASEDIR
myCHROM <- "chr13"


library(logging, quietly = TRUE); basicConfig()
logging::setLevel(level='DEBUG')  #set level of root logger
# logging::basicConfig()
# logging::setLevel('DEBUG', getHandler('basic.stdout')) # set level of handler

library(dplyr, quietly = TRUE)



library(parallel, quietly = TRUE)
library(foreach, quietly = TRUE)
library(doParallel)


# options and parameters -----

library(optparse, quietly = TRUE)
library(IRanges, quietly = TRUE)


myOptions <- list(
  make_option(c("-v", "--devmode"), action="store_true", help="svmod is installed in dev-mode!", default=FALSE),
  make_option(c("-t", "--testrun"), action="store_true", help="Do only a test run. Do not save features to disk.", default=FALSE),
  make_option(c("-m", "--margin"), type="integer", help="Use what margin left and right of candidate SV region for feature extraction [%default bp]", default=25L),
  make_option(c("-p", "--parallel"), type="integer", help="Number of Cores to use. 0: automatic, 1: no parallel. [%default]", default=2L),
  make_option(c("-x", "--randomSeed"), type="integer", help="Set random seed [%default]", default = 1234L),
  make_option(c("-r", "--maxCovRegion"), type="integer", help="Maximal number of covered regions per patient. 0 means no restriction. For now only in use when not using simulation. [%default]", default=0L)
)


myParser <- OptionParser(usage="%prog [options] patId", description="Extract features from real patients.", option_list=myOptions)
parsed_args <- parse_args(myParser, positional_arguments=TRUE)

given_opts <- parsed_args$options
given_args <- parsed_args$args


if ( length(given_args) != 1L ) {
  logging::logwarn("Please provide a single patient identifier!")
  print_help(myParser)
  q("no", status=1L, runLast=FALSE)
}

# numeric patient identifier because it is used within matrix-data (and would trigger everything being character otherwise)
patIdentifier <- as.numeric(sub(pattern = "P", replacement = "", given_args[1], fixed=TRUE))
stopifnot(is.numeric(patIdentifier), patIdentifier > 0L)

myMargin <- as.integer(ceiling(given_opts$margin))
stopifnot( is.integer(myMargin), myMargin > 0L )

nbrCores <- if ( given_opts$parallel <= 0L )
  max(1L, ceiling(sqrt(parallel::detectCores())), na.rm = TRUE) else given_opts$parallel

# how many covered regions per patient to consider at max for feature extraction?
myMaxNbrCovReg <- if (is.numeric(given_opts$maxCovRegion) && given_opts$maxCovRegion > 0L) given_opts$maxCovRegion else NULL


TIMESTAMP <- format(Sys.time(), "%Y%m%d_%H%M") #format(Sys.Date(), "%Y%m%d")
addHandler(writeToFile, file=paste0(baseDir, "log/extractFeatures_real_P", patIdentifier, "_", myCHROM, "_", TIMESTAMP, ".log"), level="DEBUG") # level="INFO")




# setup -----


if (isTRUE(given_opts$devmode)){
  library(devtools, quietly = TRUE)
  dev_mode(on = TRUE)
}

## LOAD Package: svmod
library(svmod, quietly = FALSE)
svmod::setupParallel(cpus = nbrCores)


realPatDir <- file.path(baseDir, "biotec")
stopifnot( file.exists(realPatDir) )

patBams <- list.files(path = realPatDir, pattern = paste0(patIdentifier, ".*dedup.bam$"), full.names = TRUE)
if ( length(patBams) != 2L ){
  logging::logerror("I expect precisely two BAM files for patient %s", patIdentifier)
  q("no", status=1L, runLast=FALSE)
}

rs <- given_opts$randomSeed
if (is.numeric(rs) && rs > 0L){
  rs <- ceiling(rs)
  logging::loginfo("Set random seed to %d.", rs)
  set.seed(rs)
}


# mkuhn, 2016-02-08: SV-length dist for SV-injection
# what values do I put here best? 
mySVmean <- 150
mySVsd <- 15


# build feature dataframe (train and test data) -----


logging::loginfo("Start to extract features for real patient %s.", patIdentifier)

pat.n.bamFile <- grep(pattern = "_CTRL_", patBams, value=TRUE)
pat.t.bamFile <- grep(pattern = "_T_", patBams, value=TRUE)

stopifnot( length(pat.n.bamFile) == 1L, length(pat.t.bamFile) == 1L )
stopifnot( checkFile(pat.n.bamFile), checkFile(pat.t.bamFile) )

featDat.pat <- getFeatureDataFromPatient(patId = patIdentifier, bamFile_WT = pat.n.bamFile, bamFile_MUT = pat.t.bamFile,
                                         chrom = myCHROM, localMappingMargin=1e5, propPos.train=.5, margin.bp = myMargin,
                                         SVlength_mean = mySVmean, SVlength_sd = mySVsd,
                                         maxNbrCoveredRegions = myMaxNbrCovReg, saveIntermediate=TRUE, useIntermediate=TRUE)

if ( is.null(featDat.pat) ){
  logging::logwarn("No feature data for patient %s", patIdentifier)
}

logging::loginfo("Gathered feature data featDat.pat as %s with lengths %s.", class(featDat.pat), paste(lengths(featDat.pat), collapse="-"))
# ## TESTING mkuhn, 2016-02-02
# featureDir <- file.path(realPatDir, "feature")
# try(dir.create(featureDir, showWarnings = FALSE, recursive = TRUE), silent = TRUE)
# 
# featBaseName <- paste0("featuresList_screen_testing_P", patIdentifier, "_margin")
# featureFilename <- file.path(featureDir, paste0(featBaseName, myMargin, ".rds"))
# saveRDS(featDat.pat, file=featureFilename)
# ## TESTING


# feature data frame of training and test data
features.train <- featDat.pat[["train"]]
features.test0 <- featDat.pat[["test"]]


logging::loginfo("Features of training data has first columns: %s", paste(head(names(features.train)), collapse=" - "))


#### myMargin is re-used here for estimating the status on test data
features.test <- addStatusToBiotecTest(featDat.test = features.test0, evalMargin = myMargin, regionIsMatch = FALSE, featNames = names(features.train))
logging::loginfo("Added status columns to test data for %d test data of patient %s.", NROW(features.test0), patIdentifier)

featureDat <- rbind(features.train, features.test)
logging::loginfo("Combined train/test featureDat is a %s of length %d with %d rows (%d test cases).", class(featureDat), length(featureDat), NROW(featureDat), sum(featureDat$type == 'test', na.rm = TRUE))


if (isTRUE(given_opts$testrun)){
  
  head(featureDat, n=4L)
  logging::loginfo("TEST RUN: NOT writing out %d features entries.", NROW(featureDat))
  q(save = "no", status = 0L)
}

## read in existing FEATURE data from disk
stopifnot( ! isTRUE(given_opts$testrun) )




# write out -----

# FEATURE directory
featureDir <- file.path(realPatDir, "feature")
logging::logdebug("Using feature directory: %s.", featureDir)
try(dir.create(featureDir, showWarnings = FALSE, recursive = TRUE), silent = TRUE)

featBaseName <- paste0("features_real_screen_P", patIdentifier, "_margin")
featureFilename <- file.path(featureDir, paste0(featBaseName, myMargin, ".rds"))



saveRDS(featureDat, file=featureFilename)
logging::loginfo("Feature data for real patient %s saved to disk in ::%s::", patIdentifier, featureFilename)

