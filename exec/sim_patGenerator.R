#!/usr/bin/env Rscript
# mkuhn, 20140226: Example script how to set a simulation run.
# mkuhn, 20140715: Adding optparse
# 
# Patient Simulator:
# A patient has an individual genome and an individual tumor which has its own tumor genome that stems from the patient's healthy genome.
# Idea is to be able to easily simulate new virtual "patients" which are sequenced on certain positions.
#
# Read Store (or pool): Reads are generated from pure normal and pure tumor with high coverage
# and with certain sequencing characteristics like PE insert mean and sd.
# Starting from a read store you can easily mix together what represents a sequenced normal and tumor sample (with coverage value and tumor load for tumor sample).
# Then a mapping for the sequenced samples will be done and BAM-files are created.
#######


library(logging, quietly = TRUE);
basicConfig()
setLevel(level='DEBUG')  # level of root logger
setLevel('DEBUG', getHandler('basic.stdout')) # level of handler


library(IRanges, quietly = TRUE)
library(optparse, quietly = TRUE)

myOptions <- list(
  make_option(c("-v", "--devmode"), action="store_true", help="Use svmod from dev-mode!", default=FALSE),
  make_option(c("-p", "--parallel"), type="integer", help="Number of Cores to use. Currently used only for mapping. 0: automatic, 1: no parallel. [%default]", default=2L),
  make_option(c("-x", "--randomSeed"), type="integer", help="Set random seed. 0: do not set random seed. [%default]", default=1234L),
  make_option(c("-s", "--start"), type="integer", help="Patient ID to start with [%default]", default=1L),
  make_option(c("-l", "--nbrLoci"), type="integer", help="Number of Loci per Patient [%default]", default=100L),
  make_option(c("-d", "--nbrDSV"), type="integer", help="Number of differential SVs per Patient [%default]", default=10L), 
  make_option("--SVLength_mean", type="integer", help="Mean of length of SVs [%default]", default=175L),
  make_option("--SVLength_sd", type="double", help="Standard deviation of length of SVs [%default]", default=15L),
  make_option("--preserve", action="store_true", help="Use existing patients, i.e. do *not* create new patients (patient FASTA and her readstore)! This is useful when resampling data for existing patients with new coverage or tumorload.", default=FALSE),
  make_option(c("-c", "--coverage"), type="integer", help="Coverage in sampling. [%default]", default=60L),
  make_option(c("--tumorLoad"), type="integer", help="Tumorload [%] in sampling. [%default]", default=90L)
)

myParser <- OptionParser(usage="%prog [options] nbrPatients", option_list=myOptions)
parsed_args <- parse_args(myParser, positional_arguments=TRUE)

given_opts <- parsed_args$options
given_args <- parsed_args$args


###
### Parse ARGUMENTS
###

if (length(given_args) != 1 || ! isTRUE(all(nzchar(given_args))) ){
  logwarn("Incorrect number of arguments. Please specify the number of patients to simulate.")
  print_help(myParser)
  q("no", status=1, runLast=FALSE)
}

# Read in the given positional arguments: project name and target file
nbrPats <- as.integer(given_args[1])
loginfo("%d Patients to simulate are requested.", nbrPats)

if ( is.null(nbrPats) || nbrPats <= 0){
  logwarn("Could not parse number of patients requested, %s, (arg #1).", given_args[1])
  print_help(myParser)
  q("no", status=1, runLast=FALSE)
}



patIdStart <- given_opts$start
if ( is.null(patIdStart) || patIdStart <= 0){
  logwarn("Strange start ID %d for patients.", patIdStart)
  print_help(myParser)
  q("no", status=1, runLast=FALSE)
}

patIdVect <- seq(from=patIdStart, length.out = nbrPats)

nbrLoci <- given_opts$nbrLoci #40L
nbrDSV <- given_opts$nbrDSV #2L #differential SV
SVLength_mean <- given_opts$SVLength_mean #250L
SVLength_sd <- given_opts$SVLength_sd


if (isTRUE(given_opts$devmode)){
  library(devtools, quietly = TRUE)
  dev_mode(on = TRUE)
}


# load package. Fail on error.
library(svmod, quietly = TRUE)

myBaseDir <- svmod::BASEDIR

# additional file handler
TIMESTAMP <- format(Sys.time(), "%Y%m%d_%H%M") # .. to the minute
addHandler(writeToFile, file=paste0(myBaseDir, "log/patSimulator_", TIMESTAMP, ".log"), level="INFO")


library(parallel, quietly = TRUE)
library(foreach, quietly = TRUE)


nbrCores <- if ( given_opts$parallel <= 0 ) max(1L, ceiling(sqrt(detectCores())), na.rm = TRUE) else given_opts$parallel
svmod::setupParallel(cpus = nbrCores)


rs <- given_opts$randomSeed
if (is.numeric(rs) && rs > 0L){
  rs <- ceiling(rs)
  logging::loginfo("Setting random seed to %s.", rs)
  set.seed(rs) # and below in call to processVirtualPatients
}


patData0 <- NULL
patDataFilename <- paste0(svmod::BASEDIR, "patient/patData.rds")

if ( isTRUE(file.exists(patDataFilename)) ){
  patData0 <- readRDS(patDataFilename)
  logging::loginfo("Found existing patData file with %d entries.", NROW(patData0))
}


# using default settings for NGS properties
patSim_sample.prop <- sample.prop(cov = given_opts$coverage, tumorLoad = given_opts$tumorLoad)
patSim_NGS.prop <- ngs.prop()

myMinMargin <- patSim_NGS.prop$pe.ins.mean + 6L * patSim_NGS.prop$pe.ins.sd
patData <- NULL

patData <- processVirtualPatients(patIds = patIdVect, basedir = myBaseDir, preserve = given_opts$preserve, do.mapping = TRUE,
                                  nbrCores = nbrCores, randomSeed = rs,
                                  sample.prop = patSim_sample.prop, ngs.prop = patSim_NGS.prop, 
                                  nbrLoci=nbrLoci, nbrDSV=nbrDSV, minMargin=myMinMargin, SVlength_mean = SVLength_mean, SVlength_sd=SVLength_sd)


# mkuhn, 2015-04-17: add patient info when new patients have been created, no matter if preserve or not
if ( ! is.null(patData) && NROW(patData) >= 1L ){  #! isTRUE(given_opts$preserve) && 
  newPatients <- unique(patData$patId)
  loginfo("Created in total %d new rows from %d new patients.", NROW(patData), length(newPatients))
  if ( ! is.null(patData0) && NROW(patData0) > 0L ){
    oldPatients <- patData0[! patData0$patId %in% newPatients, ] # keep those existing patients that are not in the set
    if (NROW(oldPatients) > 0L){
      loginfo("Keeping %d old rows from previous patient simulator runs.", NROW(oldPatients))
      patData <- rbind(oldPatients, patData)
    }
  }
  loginfo("Writing out patient data to file %s.", patDataFilename)
  saveRDS(patData, patDataFilename)
}#fi patData



# continue with
# extract features from mapping
