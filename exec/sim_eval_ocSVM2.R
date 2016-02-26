#!/usr/bin/env Rscript
# mkuhn, 2015-03-10
# Quick evaluation of oneclass SVM as twoclass classifier
# Using the mlr framework for this.
#
# INPuT: patient data and their feature data
# OUTPUT: evaluation stuff.. to be done

library(logging); basicConfig()
library(dplyr)



library(optparse)

myOptions <- list(
  make_option(c("-v", "--devmode"), action="store_true", help="svmod is installed in dev-mode!", default=FALSE),
  make_option(c("-p", "--parallel"), type="integer", help="Number of Cores to use. 0: automatic, 1: no parallel. [%default]", default=2L),
  make_option(c("-r", "--random_seed"), type="integer", help="Set random seed. 0: do not set random seed. [%default]", default=0L),
  #make_option(c("-s", "--start"), type="integer", help="Patient ID to start with [%default]", default=1L),
  make_option(c("-c", "--coverage"), type="integer", help="Coverage in sampling. [%defalt]", default=60L),
  make_option(c("-l", "--tumorLoad"), type="integer", help="Tumorload (in %) in sampling. [%default]", default=90L)
)

myParser <- OptionParser(usage="%prog [options]", option_list=myOptions)
parsed_args <- parse_args(myParser, positional_arguments=TRUE)

given_opts <- parsed_args$options
given_args <- parsed_args$args


###
### Parse ARGUMENTS
###

if (length(given_args) >= 1L || any(nzchar(given_args)) ){
  logwarn("This script does not have arguments, only options. Sorry.")
  print_help(myParser)
  q("no", status=1, runLast=FALSE)
}



if (isTRUE(given_opts$devmode)){
  library(devtools, quietly = TRUE)
  dev_mode(on = TRUE)
}

## LOAD Package: SVMOD
library(svmod)


nbrCores <- 1L; do.parallel <- FALSE
nbrCores <- given_opts$parallel

if (is.numeric(nbrCores) && nbrCores >= 2){
  do.parallel <- TRUE
  library(parallelMap)

  parallelStartMulticore(cpus = floor(nbrCores))
}


library(mlr)


###
### SETTINGS
###

if (is.numeric(given_opts$random_seed) && given_opts$random_seed > 0) set.seed(given_opts$random_seed)


baseDir <- svmod::BASEDIR
myNGSProp <- ngs.prop()
mySampleProp <- sample.prop(cov = given_opts$coverage, tumorLoad = given_opts$tumorLoad)  #sample.prop() ##default values



# patient information: which regions were simulated, where is the DSV
patDataFile <- file.path(getVirtualPatientPath(baseDir, "patient"), "patData.rds")
stopifnot( file.exists(patDataFile))

patData <- readRDS(patDataFile)
stopifnot( NROW(patData) >= 3L )
#head(patData, n=2)
loginfo("Read in data from %d entries from file ::%s::", NROW(patData), patDataFile)

# mkuhn, 20150420: focus on recently simulated patients
patData <- dplyr::filter_(patData, ~patId <= 50)
patData <- dplyr::mutate_(patData, SVstartPos=~targetStartPos + SVstart, SVendPos=~SVstartPos + SVlength)


# FEATURE directory
featureDir <- getVirtualPatientPath(baseDir, what = "feature", .sample.prop=mySampleProp, .ngs.prop=myNGSProp)
loginfo("Using feature directory ::%s::", featureDir)



# features dataframe
featureFilename <- file.path(featureDir, "features.rds")
featureDat <- NULL
if (! file.exists(featureFilename) ){
  logwarn("Did not find feature matrix at %s.", featureFilename)
  q(status = 1L)
}

featureDat <- readRDS(featureFilename)
stopifnot( is.data.frame(featureDat) )
loginfo("Found feature dataframe with %d existing features with %d columns.", NROW(featureDat), NCOL(featureDat))

#featDat0 <- data.frame(featureDat) %>% arrange_(~patId, ~targetStartPos, ~targetEndPos)
#loginfo("Found %d feature entries in ::%s::. Every simulated patient position that is not filtered out because of too low read coverage corresponds to a feature entry.", NROW(featDat0), featureDir)
#
## add simulation info: status and SVtype. #ZZZ chrom should be in featureDat0 as well
#featDat <- dplyr::inner_join(x=patData %>% dplyr::select_(~patId, ~chrom, ~targetStartPos, ~targetEndPos, ~status, ~SVtype),
#                      y = featDat0,
#                      by = c("patId", "targetStartPos", "targetEndPos"))



# subsampling to get the right proportions of outliers. Returns indices of selected observations in train/test data.
ssi <- subsampleInd(dat.status = (featDat$status == "germline"), prop.outl = 0.05)

featDat <- featDat[ssi$train.ind,]
# result variable
featDat.status <- ifelse(featDat$status == 'germline', 0L, 1L) # 0='germline', 1='differential SV'
featDat.SVtype <- featDat$SVtype
# remove patient info from feature dataframe
featDat <- featDat %>% 
  dplyr::select_(~-c(patId:endPos))

loginfo("Feature selection and subsampling yield %d features and %d entries.", NCOL(featDat), NROW(featDat))

#ZZZ continue here

# status variable with levels IN and OUT
featDatY2 <- cbind(status=factor(featDat.status, levels = c(0,1), labels = c("IN", "OUT")), featDat)
# keep only status and difference features  
featDatY2d <- dplyr::select_(featDatY2, ~status, ~ends_with(".d"), ~ends_with(".t"))
# mkuhn, 2015-03-16: drop eventual NA-values
featDatY2d <- featDatY2d[complete.cases(featDatY2d),]


diffSV_Task <- makeClassifTask(id = "diffSV", data = featDatY2d, target = "status") %>%
  removeConstantFeatures(perc = 0.002)

ocsvm2_Learner <- makeLearner("classif.ocsvm2")
svm2_Learner <- makeLearner("classif.svm")

# mlr::train(ocsvm2_Learner, diffSV_Task)
# mlr::train(svm2_Learner, diffSV_Task)

# tuning parameters
tuning.ps <-  makeParamSet(
  makeDiscreteParam("nu", values = c(0.01, 0.05, 0.1, 0.25, 0.5)),
  makeDiscreteParam("gamma", values = c(0.0625, 0.125, 0.5, 1, 2)) #-3:1, trafo = function(x) 2^x)
)

# grid search
tuneCtrl <-  makeTuneControlGrid()
resampDesc.inner <-  makeResampleDesc("CV", iters = 5L, stratify = TRUE)

# # tune the hyperparameters (by hand)
# res <- tuneParams(ocsvm2_Learner, task = diffSV_Task, resampling = resampDesc.inner, par.set = tuning.ps,
#                   control = tuneCtrl, measures = list(acc, setAggregation(acc, test.sd)), show.info = FALSE)

# tuned learners
ocsvm2_Learner.t <- ocsvm2_Learner %>% makeTuneWrapper(resampDesc.inner, measures = list(acc),
                                                       par.set = tuning.ps, control = tuneCtrl)
svm2_Learner.t <- svm2_Learner %>% makeTuneWrapper(resampDesc.inner, measures = list(acc),
                                                   par.set = tuning.ps, control = tuneCtrl)

resampDesc.outer <-  makeResampleDesc("CV", iters=2L, stratify = TRUE)

# resample (by hand)
# r <- resample(ocsvm2_Learner.t, diffSV_Task, resampling = resampDesc.outer, show.info = FALSE, measures = list(acc, mmce))
# r$measures.test


## model multiplexer allows to tune across different learners. Useful here? Or use benchmark instead?

## Benchmark
bmark <- benchmark(learners = list(ocsvm2_Learner.t, svm2_Learner.t),
                   tasks = diffSV_Task, resamplings = resampDesc.outer, measures = list(acc, mmce))
getBMRTuneResults(bmark, as.df = TRUE)
#getBMRPerformances(bmark, as.df = TRUE)
bmark.res <- getBMRPerformances(bmark, as.df=TRUE) %>% mutate_(cov=mySampleProp$cov, tl=mySampleProp$tumorLoad,
                                                               ngsporp=~toString(myNGSProp))


if (isTRUE(do.parallel)){
  # clean up
  parallelStop()
}

# write out result data
# RESULT directory
resultDir <- getVirtualPatientPath(baseDir, what = "result")
resultFilename <- file.path(resultDir, "result.rds")

result0 <-  if (file.exists(resultFilename)) readRDS(resultFilename) else NULL
if (NROW(result0) > 0L){
  # keep only old entries that have no match
  result0 <- dplyr::anti_join(result0, bmark.res)
}

result <- rbind(result0, bmark.res)

saveRDS(result, file=resultFilename)
#ZZZ continue here: write out performance and collect it
