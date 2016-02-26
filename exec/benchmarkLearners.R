#!/usr/bin/env Rscript
# mkuhn, 2016-01-06
# benchmark different learners using the installed svmod-package
# 
# mkuhn, 20160111: added optparse
# mkuhn, 20160202: added options for biotec patient data


library(logging, quietly = TRUE); basicConfig()
logging::setLevel(level='DEBUG')  #set level of root logger
library(dplyr, quietly = TRUE)

# chromosome not needed: we use training and test data that are available
##myCHROM <- "chr5"


# ### TESTING on mackie
# devtools::dev_mode(on=TRUE)

# logging::basicConfig()
# logging::setLevel('DEBUG', getHandler('basic.stdout')) # set level of handler

# for quick setup
if (FALSE){
  baseDir <- svmod::BASEDIR
  nbrCores <- 2L
  ngsProp <- ngs.prop()
  sampleProp <- sample.prop()
  
  myMargin <- 25L
  shortEval <- TRUE
  myMaxNbrCovReg <- 3L
}


library(parallel, quietly = TRUE)
library(foreach, quietly = TRUE)
library(doParallel)



library(optparse, quietly = TRUE)
library(IRanges, quietly = TRUE)


myOptions <- list(
  make_option(c("--biotec"), action = "store_true", help = "Flag if given patient is patient ID from Biotec"),
  make_option(c("--patId"), type = "integer", help = "Take features from this specified single patient ID. In particular useful for single patient analysis (for instance for real data)"),
  make_option(c("-c", "--coverage"), type="integer", help="Coverage in sampling. [%default]", default=60L),
  make_option(c("-l", "--tumorLoad"), type="integer", help="Tumorload [%] in sampling. [%default]", default=90L),
  make_option(c("-m", "--margin"), type="integer", help="Use what margin left and right of identified regions for feature extraction [%default bp]", default=25L),
  make_option(c("-p", "--parallel"), type="integer", help="Number of Cores to use. 0: automatic, 1: no parallel. [%default]", default=3L),
  make_option(c("--tuneMode"), type="character", help='Tune mode (grid, random, irace). With grid option tuneIter is not used.', default='grid'),
  make_option(c("-t", "--tuneIter"), type="integer", help="Number of tuning iterations. Not used for grid tuning.", default = 60L),
  make_option(c("-x", "--randomSeed"), type="integer", help="Set random seed [%default]", default = 1234L),
  make_option(c("-i", "--useSimInfo"), action = "store_true", help = "Flag if you want to use simulation information about somatic SVs for learner evaluation, i.e. avoid screening step. [%default]"),
  make_option(c("--short"), action = "store_true", help = "Flag if you want only a short evaluation [%default]"),
  make_option(c("-s", "--suffix"), type="character", help="Suffix string for benchmark filename", default="")
)


myParser <- OptionParser(usage="%prog [options]", description="Run benchmark of machine learning algorithms on simulated or biotec data", option_list=myOptions)
parsed_args <- parse_args(myParser, positional_arguments=TRUE)

given_opts <- parsed_args$options
given_args <- parsed_args$args


if ( length(given_args) > 0L || any(nzchar(given_args)) ){
  logwarn("This script does not have any arguments, only options. Sorry.")
  print_help(myParser)
  q("no", status=1L, runLast=FALSE)
}


shortEval <- given_opts$short

myMargin <- as.integer(ceiling(given_opts$margin))
stopifnot( is.integer(myMargin), myMargin > 0L )

nbrCores <- if ( given_opts$parallel <= 0L )
  max(1L, ceiling(sqrt(detectCores())), na.rm = TRUE) else given_opts$parallel


tuneMode <- as.character(given_opts$tuneMode)
tuneIter <- as.integer(ceiling(given_opts$tuneIter))
suffix <- as.character(given_opts$suffix)

patId <- as.numeric(given_opts$patId)



TIMESTAMP <- format(Sys.time(), "%Y%m%d_%H%M")   
logFilename <- sprintf("benchmark_%s_%s.log", if (isTRUE(given_opts$biotec)) paste0("biotec", given_opts$patId) else "sim", TIMESTAMP)
addHandler(writeToFile, file=file.path(svmod::BASEDIR, "log", logFilename), level="INFO")  #myCHROM, "_",



library(mlr)
library(svmod)

rs <- given_opts$randomSeed
if ( is.numeric(rs) && rs > 0L ){
  rs <- ceiling(rs)
  logging::loginfo("Set random seed to %d", rs)
  set.seed(rs)
}
 

svmod::setupParallel(cpus = nbrCores, myLevel = "mlr.tuneParams")

# build list input
sources <- if ( isTRUE(given_opts$biotec) ){
  stopifnot( ! is.null(given_opts$patId) )
  logging::loginfo("~Start Benchmarking for Biotec data of patient %s in tune-mode %s with %d tuning iterations and with suffix [%s]~", patId, tuneMode, tuneIter, suffix)
  list(isSimulation=FALSE, useSimSVInfo=FALSE, patIds=patId, margin.bp=myMargin)
} else {
  ngsProp <- ngs.prop() # is fix
  sampleProp <- sample.prop(cov = given_opts$coverage, tumorLoad = given_opts$tumorLoad)
  
  logging::loginfo("~ Start Benchmarking for simulated data in tune-mode %s with %d tuning iterations and with suffix [%s] on sample prop %s ~", tuneMode, tuneIter, suffix, toString(sampleProp))
  list(isSimulation=TRUE, sampleProp=sampleProp, ngsProp=ngsProp, patIds=patId, useSimSVInfo=given_opts$useSimInfo, margin.bp=myMargin)
}

# bmrRes <- evalLearners2(sample.prop=sampleProp, ngs.prop=ngsProp,
#                         tuneMode = tuneMode, tune.niter = tuneIter, k=10L,
#                         saveEvaluation = TRUE, suffix = suffix, verbose = TRUE)
bmrRes <- evalLearners2(sources = sources, tuneMode = tuneMode, tune.niter = tuneIter, shortEval = shortEval,
                        k = 10L, saveEvaluation = TRUE, suffix = suffix, verbose = TRUE)



if (is.null(bmrRes))
  logging::logwarn("Evaluation failed.") else {
    logging::loginfo("Finished evaluation with sources %s", paste(names(sources), sources, sep=":", collapse = " - "))
    bmrRes
  }

logging::loginfo("~Fini~")
