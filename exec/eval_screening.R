#!/usr/bin/env Rscript
# mkuhn, 20151201: run clipped-base peak screening on simulated covered regions and evaluate (assess) its performance (SEN/SPEC).
# Influencing factors on parameters are parameters
# + SV type
# + SV length
# + coverage
# + tumor load TL
#
# Because of simulation we know if a mutation is present and where it precisely was simulated.
# We use the "realpatient" screening method to define the positive candidate regions and can hence decide if it is a correct hit or not (using a myEvalMargin)
# INput: simulated patient data (mapping (BAM file) [T]-samples corresponding to different coverage-tumorload simulation)
# OUTput: performance measures for screening step


library(logging, quietly = TRUE); basicConfig()
library(dplyr, quietly = TRUE)

library(optparse, quietly = TRUE)
library(IRanges, quietly = TRUE)


myOptions <- list(
  make_option(c("-v", "--devmode"), action="store_true", help="svmod is installed in dev-mode!", default=FALSE),
  make_option(c("-c", "--coverage"), type="integer", help="Coverage in sampling. [%default]", default=60L),
  make_option(c("-r", "--chrom"), type="character", help="Chromosome. [%default]", default="chr5"),
  make_option(c("-l", "--tumorLoad"), type="integer", help="Tumorload [%] in sampling. [%default]", default=90L),
  make_option(c("-m", "--evalmargin"), type="integer", help="Use what margin left and right of predicted SV-interval to decide if screening produced a hit or not [%default bp]", default=25L),
  make_option(c("-p", "--parallel"), type="integer", help="Number of Cores to use. 0: automatic, 1: no parallel. [%default]", default=2L),
  make_option(c("-s", "--startPatID"), type="integer", help="First Patient ID to do feature extraction [%default]", default=-1L),
  make_option(c("-e", "--endPatID"), type="integer", help="Last Patient ID to do feature extraction [%default]", default=-1L)
)

myParser <- OptionParser(usage="%prog [options]", option_list=myOptions)
parsed_args <- parse_args(myParser, positional_arguments=TRUE)

given_opts <- parsed_args$options
given_args <- parsed_args$args


if (length(given_args) > 0L || any(nzchar(given_args)) ){
  logwarn("This script does not have any arguments, only options. Sorry.")
  print_help(myParser)
  q("no", status=1L, runLast=FALSE)
}


myChrom <- given_opts$chrom


startID <- given_opts$startPatID
endID <- given_opts$endPatID


# first check endID, then startID
myPatIds <- NULL 
if ( is.numeric(endID) && endID >= 1L ){
  myPatIds <- seq_len(endID)
  if ( is.numeric(startID) && startID > -1L) myPatIds <- myPatIds[which(myPatIds >= startID)]
  loginfo("Filtered %d patients for screen evaluation to %s.", length(myPatIds), paste(myPatIds, collapse="-"))
}


if (isTRUE(given_opts$devmode)){
  library(devtools, quietly = TRUE)
  dev_mode(on = TRUE)
}


library(parallel, quietly = TRUE)
library(foreach, quietly = TRUE)


## LOAD Package:
## svmod

library(svmod, quietly = FALSE)



TIMESTAMP <- format(Sys.time(), "%Y%m%d_%H%M") #format(Sys.Date(), "%Y%m%d")
addHandler(writeToFile, file=paste0(svmod::BASEDIR, "log/evalScreening_", TIMESTAMP, ".log"), level="DEBUG")


nbrCores <- if ( given_opts$parallel <= 0L )
  max(1L, ceiling(sqrt(detectCores())), na.rm = TRUE) else given_opts$parallel
svmod::setupParallel(cpus = nbrCores)


##
## SETUP
##
baseDir <- svmod::BASEDIR
myNGSProp <- ngs.prop()
mySampleProp <- sample.prop(cov = given_opts$coverage, tumorLoad = given_opts$tumorLoad)

myTargetBedFile <- file.path(baseDir, "refGen", "TruSeq_exome_targeted_regions.hg19.bed")


myEvalMargin <- as.integer(ceiling(given_opts$evalmargin))
stopifnot( is.integer(myEvalMargin), myEvalMargin >= 0L )

logging::loginfo("Start with evaluation of screening for patients with NGS %s and sampling %s and with evaluation margin of %d.",
                 toString(myNGSProp), toString(mySampleProp), myEvalMargin)


evalScreenDat <- evalScreen(baseDir = baseDir, patIds = myPatIds, myTargetBedFile=myTargetBedFile, 
                            myChrom = myChrom, sample.prop = mySampleProp, ngs.prop = myNGSProp,
                            margin.bp = myEvalMargin, saveEvaluation = TRUE)



# margin 
evalScreenDat %>% 
  xtabs(~SVtruth + screenCall, data=.) %>% 
  prop.table(1L) %>% 
  knitr::kable(caption="Sensitivity and Specificity for margin") %>% 
  print


# region
evalScreenDat %>% 
  xtabs(~SVtruth + screenCallReg, data=.) %>% 
  prop.table(1L) %>% 
  knitr::kable(caption="Sensitivity and Specificity for region") %>% 
  print

cat("\n ~ Fine ~\n")
