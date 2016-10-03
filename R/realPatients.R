# mkuhn, 2015-03-21
# Functions for analysis of real patient data. This includes the SV-injection.




#' Extract the target region coordinates from a target BED file.
#' 
#' This is a little helper function.
#' @param targetBedFile the target regions as BED file
#' @param chrom chromosome to work on
#' @param extra.margin integer. number of basepairs to widen the target intervals
#' @param reduce logical flag if we want to reduce the GRanges.
#' @return Genomic ranges of target region as \code{GRanges} object
#' @export
getTargetRegions <- function(targetBedFile=file.path(BASEDIR, "refGen", "TruSeq_exome_targeted_regions.hg19.bed"), chrom=CHROM, extra.margin=EXTRA_MARGIN_TARGET_REGIONS,
                             reduce=FALSE) {
  
  stopifnot( checkFile(targetBedFile) )
  
  # form enrichment GRanges object (from BED file)
  targetBed <- read.delim(file = targetBedFile, header = FALSE)
  stopifnot( NCOL(targetBed) >= 3L )
  colnames(targetBed)[1L:3L] <- c("chr", "start", "end") #, "id", "width", "strand")
  
  # mkuhn, 2016-10-03: drop NA-rows
  targetBed <- droplevels(targetBed[!is.na(targetBed$start) & ! is.na(targetBed$end), ])
  # BED-file is 0-based, GRanges are 1-based
  targetBed$start <- targetBed$start + 1L
  targetRegions <- GenomicRanges::makeGRangesFromDataFrame(targetBed, ignore.strand = TRUE, keep.extra.columns = FALSE)
  stopifnot( length(targetRegions) >= 1L )
  try(rm(targetBed))
  
  # use only a single fixed chromosome (chr5)
  targetRegions.chr <- targetRegions[seqnames(targetRegions) == chrom]
  # widen the intervals
  GenomicRanges::start(targetRegions.chr) <- pmax(1L, GenomicRanges::start(targetRegions.chr) - extra.margin)
  GenomicRanges::end(targetRegions.chr) <- GenomicRanges::end(targetRegions.chr) + extra.margin  # ZZZ mind the chromosome end?!
  
  if (isTRUE(reduce)) targetRegions.chr <- GenomicRanges::reduce(targetRegions.chr, drop.empty.ranges=TRUE, ignore.strand=TRUE, min.gapwidth=1L)
  
  return( targetRegions.chr )
}







#' #' Extract positions from peak info
#' #' @return start and end coordinate of SV candidate
#' #' @note Obsolete
#' getPeakRange <- function(pInfo){
#'   
#'   peakRange <- NULL
#'   
#'   if (is.data.frame(pInfo) && all( c("pos.3p", "pos.5p") %in% names(pInfo)) && NROW(pInfo) >= 1L )
#'     peakRange <- range(pInfo$pos.5p, pInfo$pos.3p)
#'   
#'   return(peakRange)
#' }




#' Select features from a patient's tumor and normal mapping.
#' 
#' First, identifies sufficiently covered loci within the target region and then apply a (clipped-reads based) screening to identify SV-candidates.
#' Taking non-hit regions as basis for training data: using SV-injection we come to positive training data. TL is given or taken from simulation setting.
#' Screening-positive candidate regions are basis for test data. If we have a simulated patient we know the truth at these test regions.
#' 
#' For now, I use for training only regions where there is no clipped-base peak high enough. I should check that we have enough training data.
#' If not enough training data one could think to use also positive screening regions, but at a non-peak spot in the region..
#' 
#' @param patId Patient ID for identification of the feature data (for simulation as numeric e.g. 17). Real patients will typically have a string identifier (e.g. L761)
#' @param bamFile_WT filename for BAM-file of normal sample
#' @param bamFile_MUT filename for BAM-file of tumor sample
#' @param TL tumor load used for SV-injection. If tumor probe is from simulation, then the value is taken from there and this option here ignored.
#' @param chrom restriction to specific chromosome
#' @param targetBedFile coordinates of target regions in BED format
#' @param margin.bp window width (in nucleotides) to the left and right of the SV used for feature extraction. Defaults to 1/4 of read length
#' @param localMappingMargin integer. radius defining local mapping neighbourhood around SV. 
#' @param maxNbrCoveredRegions maximal number of covered regions where from where to extract features. Defaults to \code{NULL}, i.e. take all available.
#' @param propPos.train proportion of positive cases (through injected SVs) within training data. Defaults to 50\%.
#' @param saveIntermediate flag if we want to save intermediate output. Might be useful for real patient analysis.
#' @param useIntermediate flag if we want to quickly load intermediate results. Might be useful for real patient analysis.
#' @return list of feature data, for training (based on screening negative regions) and for testing (the positive candidate regions) or \code{NULL} if no covered regions available
#' @export
getFeatureDataFromPatient <- function(patId, bamFile_WT, bamFile_MUT, TL=0.9 * .5, #isSimulatedPatient=FALSE, #moved to sim_extractFeatures
                                      chrom=CHROM, #seqBaseFile=paste0(BASEDIR, "refGen/hg19_chr5_truseq.b1000.fa"),
                                      targetBedFile = file.path(BASEDIR, "refGen", "TruSeq_exome_targeted_regions.hg19.bed"),
                                      margin.bp=ceiling(READL/4L), localMappingMargin=1e5,
                                      maxNbrCoveredRegions=NULL, propPos.train=.5,
                                      SVlength_mean, SVlength_sd, # for SV-injection
                                      #file.path(BASEDIR, "biotec", "intermediate")
                                      saveIntermediate=FALSE, saveIntermediateDir=dirname(bamFile_MUT), useIntermediate=FALSE) {  
  
  
  if (FALSE){
    patId <- 48 
    
    # paste0("mapping/cov20_TL50/PEm200sd25_R100/", getPatStr(patId), "_N.s.bam"))
    # file.path(BASEDIR, paste0("mapping/cov20_TL50/PEm200sd25_R100/", getPatStr(patId), "_T.s.bam"))
    bamFile_WT  <- file.path(BASEDIR, paste0("mapping/cov60_TL90/PEm200sd20_R100/", getPatStr(patId), "_N.s.bam"))
    bamFile_MUT <- file.path(BASEDIR, paste0("mapping/cov60_TL90/PEm200sd20_R100/", getPatStr(patId), "_T.s.bam"))
    chrom <- "chr5" #"chr13" #
    targetBedFile = file.path(BASEDIR, "refGen", "TruSeq_exome_targeted_regions.hg19.bed")
    margin.bp=ceiling(READL/4); maxNbrCoveredRegions=50L; propPos.train=.5; localMappingMargin=1e5
    TL <- 0.9 #0.9 * .5 #heterozygous
    SVlength_mean <- 175; SVlength_sd <- 15
    useIntermediate <- TRUE
    saveIntermediate <- TRUE
    saveIntermediateDir <- dirname(bamFile_MUT) #file.path(BASEDIR, "biotec", "intermediate")
  }
  
  
  logging::loginfo("Looking at patient %s with mappings %s and %s in %s.", patId, basename(bamFile_WT), basename(bamFile_MUT), dirname(bamFile_MUT))
  
  
  # mapping info for patient -------
  
  # mkuhn, 2016-01-30: intermediate stuff
  if (isTRUE(saveIntermediate)){
    try(dir.create(saveIntermediateDir, showWarnings = FALSE, recursive = TRUE))
  }
  
  
  mappingInfoFile <- file.path(saveIntermediateDir, paste0(getPatStr(patId), "_mappingInfo.rds"))
  if ( isTRUE(useIntermediate) && file.exists(mappingInfoFile) ){
    patMappingInfo <- readRDS(file = mappingInfoFile)
    patMappingNGS <- patMappingInfo[["ngsProp"]]
    
  } else {
    patMappingInfo <- getPatientMappingInfo(bamFile_WT = bamFile_WT, bamFile_MUT = bamFile_MUT, wBioconductor=TRUE)
    patMappingNGS <- estimateNGSPropFromMapping(bamFile = bamFile_WT, wBioconductor=TRUE)
    patMappingInfo[["ngsProp"]] <- patMappingNGS
    
    # mkuhn, 2016-01-30: intermediate stuff
    if (isTRUE(saveIntermediate)){
      saveRDS(patMappingInfo, file = mappingInfoFile)
    }
  }
  
  
  clippedBaseInfoWT <- patMappingInfo[["clippedBases"]]
  clippedPropThreshold <- clippedBaseInfoWT[["peakMinThreshold"]]
  clippedNormalMaxThreshold <- clippedBaseInfoWT[["normalMaxThreshold"]]
  
  
  
  
  # find covered regions ------
  
  # covered regions should be wide enough to have at least enough space for margin left and right of SV
  stopifnot( 2L * (margin.bp+1L) < MIN_COV_REGION_WIDTH )
  
  
  # get covered regions in WT and MUT BAM file (GRanges object)
  
  coverageRegFile <- file.path(saveIntermediateDir, paste0(getPatStr(patId), "_", chrom,"_covReg.rds"))
  if ( isTRUE(useIntermediate) && file.exists(coverageRegFile) ){
    coveredRegions <- readRDS(file = coverageRegFile)
    logging::loginfo("Have loaded %d covered regions from file %s.", length(coveredRegions), coverageRegFile)
  } else {
    coveredRegions  <- findCoveredRegionsInTarget(bamFile_WT = bamFile_WT, bamFile_MUT = bamFile_MUT, 
                                                  chrom = chrom, targetBedFile = targetBedFile, 
                                                  k = 23L, minWidth = MIN_COV_REGION_WIDTH, minGapWidth = READL, minCov_WT = MINCOV_WT, minAvgCov_MUT = 5L)
    
    logging::loginfo("Found %d covered regions for patient %s. Focus on the %d regions that lie in %s",
                     length(coveredRegions), patId,
                     length(which(chrom == as.character(GenomicRanges::seqnames(coveredRegions)))), chrom)
    coveredRegions <- coveredRegions[which(chrom == as.character(GenomicRanges::seqnames(coveredRegions)))]
    
    if (isTRUE(saveIntermediate)){
      saveRDS(coveredRegions, file = coverageRegFile)
      logging::loginfo("Saved %d covered regions to file %s.", length(coveredRegions), coverageRegFile)
    }
  }
  
  # every covered region is on target chromosome
  stopifnot( all(chrom == as.character(GenomicRanges::seqnames(coveredRegions))) )
  
  
  # restrict number of covered regions? (after shuffling)
  if ( is.numeric(maxNbrCoveredRegions) && maxNbrCoveredRegions > 0L ) coveredRegions <- head(coveredRegions[sample(NROW(coveredRegions)),], n=round(maxNbrCoveredRegions))
  
  
  # mkuhn, 2016-01-05: check if we have coverage
  if ( NROW(coveredRegions) == 0L ){
    logging::logwarn("Found no covered regions for patient %s..", patId)
    return(invisible(NULL))
  }
  
  
  
  
  # assess covered regions (w.r.t clipped base peaks in tumor probe) -----------
  
  clipRegFile <- file.path(saveIntermediateDir, paste0(getPatStr(patId), "_", chrom, "_clipReg.rds"))
  if ( isTRUE(useIntermediate) && file.exists(clipRegFile) && is.null(maxNbrCoveredRegions) ){
    clippedRes_MUT <- readRDS(file = clipRegFile)
    logging::loginfo("Have loaded %d clipped peak assessments for covered regions from file %s.", NROW(clippedRes_MUT), clipRegFile)
  } else {
  
    clippedRes_MUT <- foreach::foreach(i=seq_along(coveredRegions), .combine='rbind') %do% {
      
      # assess clipped-base peaks in covered region: either 0 or 1 or 2 entries (3p-5p peak pair) that correspond to highest potential SV
      regionStart <- IRanges::start(coveredRegions)[i]
      regionEnd <- IRanges::end(coveredRegions)[i]
      
      logging::logdebug("Assess covered region %d for patient %s w.r.t. clipped base peaks", i, patId)  
      
      # findClippedPeaks_G_OLD
      # windowSize=31L, minThresholdAbs=2L, minThreshold=0.05, thresholdQuantile=0.95){ #, minTargetBorderDist=1L){
      assessClippedBasesInRegion(bamFile = bamFile_MUT, chrom = chrom,
                                 minThreshold = clippedPropThreshold, do.plot=FALSE,
                                 .targetStartPos = regionStart, .targetEndPos = regionEnd)[["clippedInfo"]]
      #cInfo <- peaks_MUT[["clippedInfo"]]
      #pInfo <- peaks_MUT[["peakInfo"]]
      #negCandPos <- peaks_MUT[["nonPeakPos"]]
      
      #cInfo
    }#foreach clippedRes
    
    if (isTRUE(saveIntermediate) && is.null(maxNbrCoveredRegions)){
      saveRDS(clippedRes_MUT, file = clipRegFile)
      logging::loginfo("Saved %d results of clipped-peak assessment of all covered regions to file %s!", NROW(clippedRes_MUT), clipRegFile)
    }
  }

  # just in case: column vector to matrix (when there is only a single covered region considered)
  if (length(coveredRegions) == 1L) clippedRes_MUT <- t(clippedRes_MUT)
  
  
  stopifnot( NROW(coveredRegions) == NROW(clippedRes_MUT), all(c("nbrPeaks", "relS_5p.99") %in% colnames(clippedRes_MUT)) )
  # mkuhn, 2016-02-06: check that we are really taking about the same regions!
  stopifnot( all( IRanges::start(coveredRegions) == clippedRes_MUT[ , "regionStart"]) ) 
  
  
  # number of negative cases for a positive case in the training data
  # preliminary counts and numbers (in case of duplicates, see below)
  # for 1 positive case I want nbrNeg0 negative cases
  nbrNeg0 <- max(1L, ceiling(1L / propPos.train - 1L)) # calculate nbr of neg cases per 1 positive case
  
  # minimal covered length necessary to generate nbrNeg0 different negative cases and one positive case
  # +25 to have same extra safety margin
  myMinExtendedSV_width <- 25L + 2L * margin.bp + SV_MIN_LENGTH + nbrNeg0 + 1L 
  
  # mkuhn, 2016-01-09: filter out too short covered regions as they may give weird mapping pattern features
  lengthFilterInd <- GenomicRanges::width(coveredRegions) >= myMinExtendedSV_width  #2L * (margin.bp + READL + SIM_SV_LENGTH_MEAN + 1L)
  coveredRegions <- coveredRegions[lengthFilterInd, ]
  clippedRes_MUT <- clippedRes_MUT[lengthFilterInd, ]
  
  logging::loginfo("%d covered regions remain after filtering out covered regions shorter than minimum width.", length(coveredRegions))
  stopifnot( NROW(coveredRegions) == NROW(clippedRes_MUT) )
  
  # mkuhn, 2016-01-05: check if we have coverage
  if ( NROW(coveredRegions) == 0L ){
    logging::logwarn("Found no covered regions of sufficient length for patient %s..", patId)
    return(invisible(NULL))
  }
  
  
  posCandidateRegion.ind <- which(clippedRes_MUT[, "nbrPeaks"]  > 0L)
  negCandidateRegion.ind <- which(clippedRes_MUT[, "nbrPeaks"] == 0L)
  
  logging::loginfo("Patient %s has %d clipped-base peak positive regions from %d considered covered regions.",
          patId, length(posCandidateRegion.ind), length(coveredRegions))
  
  # we want to have some regions with no clipped base peaks!
  # needs to be longer than 1 because otherwise sample(negCandidateRegion.ind) function behaves differently [mkuhn, 2016-01-06: not sure where this is really needed]
  if ( length(negCandidateRegion.ind) <= 1L ){
    logging::logwarn("Found too few negative candidate regions within covered regions for patient %s..", patId)
    return(invisible(NULL))
  }
  

  
  
  
  
  # build training data ------------------------  
  
  stopifnot( length(negCandidateRegion.ind) > 1L )
  
  
  # preparation for SV-injection: choose a SV for each negative covered region
  mySV.types <- sample(SV_TYPES, size = length(negCandidateRegion.ind), replace = TRUE)
  # mkuhn, 2016-02-04: I do not use a fixed distribution of SV-lengths but randomly draw a start position and sample SV-length from what is available to respect the margins.
  # The covered regions are required to allow to harbor a minimal SV-length
  #mySV.lengths <- pmax(SV_MIN_LENGTH, abs(ceiling(rnorm(length(negCandidateRegion.ind), mean=SIM_SV_LENGTH_MEAN, sd=SIM_SV_LENGTH_SD)))+1L)
  
  # mkuhn, 2016-01-10: extract TL if from simulation run
  # # assume a heterozygous pure tumor
  parsedTumorSampProp <- parseSampleProp(bamFile_MUT)
  if (! is.null(parsedTumorSampProp)){
    TL <- parsedTumorSampProp$tumorLoad / 100L # keep TL as given in simulation. (I used to half it because of heterozygosity but it should really stay as given)
  }
  logging::loginfo("Set tumor-load for SV-injection to %.2f.", TL)
  stopifnot( exists("TL"), TL > 0L, TL <= 1L )
  
  
  
  # run across negative regions (i.e. regions that showed no clipped base peaks)
  # Each screen-negative region gives rise to
  # + some negative examples that differ in the locus within the region where features are extracted
  # + one positive example with feature extraction at a locus where an SV has been injected
  #featureDat.train <- foreach::foreach(cr.ind=negCandidateRegion.ind, j=seq_along(negCandidateRegion.ind),  .combine='rbind') %do% {
  ##.combine=list, .multicombine = TRUE) %do% {
  #### DEBUGGING
  featureDat.train <- vector(mode = "list", length = length(negCandidateRegion.ind))
  for (j in seq_along(negCandidateRegion.ind)){
    cr.ind <- negCandidateRegion.ind[j]
    
    # covered region start and end
    covRegStart <- IRanges::start(coveredRegions[cr.ind])
    covRegEnd   <- IRanges::end(coveredRegions[cr.ind])
    
    
    # select random positions within region where to extract features
    # (1st position for positive case after SV-injection, rest for negative cases)
    
    
    # relative position of SV within the covered region
    # mkuhn, 2015-12-07: I wanted to make sure that the injected SV lies completely inside the covered region..
    # I relax this assumption because the SV can be long, even longer than the covered region..
    #mySVStart0_rel <- sample(1L+abs(width(coveredRegions[cr.ind]) - 2L * (mySV.lengths[j] + margin.bp)), size = nbrNeg0+1L)
    # mkuhn, 2016-01-09: check for also minimal SV length as well
    myCovRegWidth <- GenomicRanges::width(coveredRegions[cr.ind])
    
    logging::logdebug("Pat %s, neg covered region no %d with region-index %d with width %d", patId, j,  cr.ind, myCovRegWidth)
    
    # Check if covered region is wide enough for at least minimal SV length
    #+need enough space (and also enough choice for next sample()-call )
    # should always be the case after the covReg filter step above
    if ( myCovRegWidth - myMinExtendedSV_width >= nbrNeg0 + 1L ){
      
      # relative SV-start positions within covered region
      mySVStart0_rel <- margin.bp + sample(myCovRegWidth - myMinExtendedSV_width, size = nbrNeg0+1L)
      # absolute positions of first SV-base
      mySVStart0 <- covRegStart + mySVStart0_rel
      
      # SV interval end (both for positive and the negative cases).
      # we guarantee a minimum length of the SV
      # mkuhn, 2016-02-04: mySVEnd0: position of last SV-position (not including margin) [previously mySVEnd0 contained also the margin to the right]
      stopifnot( all(myCovRegWidth - mySVStart0_rel - SV_MIN_LENGTH  - margin.bp >= 1L) )
      mySVEnd0 <- numeric(nbrNeg0 + 1L)
      # propose SV-lengths according to wished distribution
      mySVlengths <- pmin(SV_MAX_LENGTH, pmax(SV_MIN_LENGTH, round(rnorm(n = length(mySVEnd0), mean = SVlength_mean, sd = SVlength_sd))))
      # old code was: sample(mySV.lengths, size = length(mySVStart0), replace = TRUE) + margin.bp
      
      for (i in seq_along(mySVEnd0)){
        # if SV too long than sample from what is available
        if ( mySVlengths[i] >= myCovRegWidth - mySVStart0_rel[i] - margin.bp ){
          mySVlengths[i] <- min(SV_MAX_LENGTH, SV_MIN_LENGTH + sample(myCovRegWidth - mySVStart0_rel[i] - SV_MIN_LENGTH - margin.bp, size=1L))
        }
        logging::loginfo("Chose SV-length of %d for SV-injection.", mySVlengths[i])
        
        #mySVStart0[i] + SV_MIN_LENGTH + sample(myCovRegWidth - mySVStart0_rel[i] - SV_MIN_LENGTH  - margin.bp, size = 1L)
        mySVEnd0[i] <- mySVStart0[i] + mySVlengths[i] - 1L
      }#rof
      
      if (any(mySVEnd0 > covRegEnd)){
        logging::logwarn("SHOULD NOT HAPPEN any more! %d regions are cut off at the region end by %s, with remaining %s bp in SV.", sum(mySVEnd0 > covRegEnd),
                         paste(mySVEnd0[mySVEnd0 > covRegEnd] - covRegEnd, collapse="::"), paste(covRegEnd - mySVStart0[mySVEnd0 > covRegEnd] + 1L, collapse="::"))
      }
      stopifnot( length(mySVStart0) == length(mySVEnd0) )
      
      mySVIntervals0 <- IRanges::IRanges(start = mySVStart0, end = mySVEnd0)
      myExtSVIntervals0 <- IRanges::IRanges(start = pmax(covRegStart, mySVStart0 - margin.bp),
                                            end   = pmin(covRegEnd,  mySVEnd0 + margin.bp))
      
      # mkuhn, 2016-01-09: filter SVs below minimal length of SV
      #myExtSVIntervals0 <- myExtSVIntervals0[IRanges::width(myExtSVIntervals0) >= SV_MIN_LENGTH + 2L*margin.bp,]
      stopifnot( all(IRanges::width(mySVIntervals0) >= SV_MIN_LENGTH), all(IRanges::width(myExtSVIntervals0) >= SV_MIN_LENGTH + 2L * margin.bp) )
      
      
      # drop duplicates in the negative examples
      # and update vectors accordingly
      duplicateMask <- duplicated(myExtSVIntervals0) | duplicated(mySVIntervals0)
      mySVIntervals <- mySVIntervals0[! duplicateMask]
      myExtSVIntervals <- myExtSVIntervals0[! duplicateMask]
      try(rm(mySVStart0_rel, mySVStart0, mySVEnd0))
      #mySVStart_rel <- mySVStart0_rel[! duplicateMask]
      #mySVStart <- mySVStart0[! duplicateMask]
      #mySVEnd <- mySVEnd0[! duplicateMask]
      nbrNeg <- length(myExtSVIntervals) - 1L # -1 because of the one positive case
      # at least one negative case (for the coming indices)
      stopifnot( nbrNeg >= 1L )  
      
      if (any(IRanges::start(myExtSVIntervals) == covRegStart, IRanges::end(myExtSVIntervals) == covRegEnd)){
        logging::loginfo("At least one of the extended %d injection SV-intervals lies on covRegStart or covRegEnd. Effective margins are %s and %s.",
                         length(myExtSVIntervals), paste(IRanges::start(mySVIntervals) - IRanges::start(myExtSVIntervals), collapse="-"), 
                         paste(IRanges::end(myExtSVIntervals) - IRanges::end(mySVIntervals), collapse="-"))
      }
      
      logging::loginfo("Widths of intervals (first will become positive SV-training example) are [ %s ].", paste(IRanges::width(mySVIntervals), collapse=" :: "))
      
      
      
      
      #Idea: I could use also the locus of simulated positive case (via SV-injection) as negative case (w/o SV-injection?)
      mySVIntervals.pos <- mySVIntervals[1L]
      mySVIntervals.neg <- mySVIntervals[-1L] # SV intervals for NEGATVE cases. 
      myExtSVIntervals.pos <- myExtSVIntervals[1L]
      myExtSVIntervals.neg <- myExtSVIntervals[-1L] # SV intervals for NEGATVE cases. 
      stopifnot( length(myExtSVIntervals.pos) == 1L )
      
      
      # inject SV at negative region to create a positive case!
      # Idea: use SV length that fits to the covered region at hand. and then also take only test cases with peaks nicely within covered region..
      mySV <- SV(type = mySV.types[j], chrom = chrom, posL = IRanges::start(mySVIntervals.pos), seqLen = IRanges::width(mySVIntervals.pos))
      injectedBAM <- injectSVInMapping(bamFile = bamFile_MUT, sv = mySV, localMappingMargin = localMappingMargin, ngsProp=patMappingNGS, TL = TL)
      
      if ( file.exists(injectedBAM[["mergedBam"]]) && file.exists(injectedBAM[["localBam"]]) ){
        # assess the tumor BAM file at this initially negative region again, this time **with** injected SV, locally around the inj SV
        clipInfo.pos <- assessClippedBasesInRegion(bamFile = injectedBAM[["mergedBam"]], chrom = chrom, do.plot=F,
                                                   minThreshold = clippedPropThreshold,  windowSize=31L, minThresholdAbs=2L, thresholdQuantile=0.95,
                                                   .targetStartPos = covRegStart, .targetEndPos = covRegEnd)[["peakInfo"]] #[["clippedInfo"]]
        logging::loginfo("Screening on positive training example (after SV-injection) detected %s peaks with %s posDiff.", clipInfo.pos["nbrPeaks"], clipInfo.pos["posDiff"])
        
        
        # positive case (via SV-injection) from local mappings
        # ZZZ I could reuse the WT calculations of getPatientMappingInfo
        patMappingInfo_injSV <- getPatientMappingInfo(bamFile_WT = injectedBAM[["localBam"]], bamFile_MUT = injectedBAM[["mergedBam"]], chrom = chrom)
        patMappingInfo_injSV$rc.bam.ratio <- patMappingInfo$rc.bam.ratio # keep old read count ratio
        features.pos <- getFeaturesFromMapping(patId=patId, bamFile_WT = injectedBAM[["localBam"]], bamFile_MUT = injectedBAM[["mergedBam"]],
                                               patMappingInfo = patMappingInfo_injSV,
                                               chrom = chrom, startPos = IRanges::start(myExtSVIntervals.pos), endPos = IRanges::end(myExtSVIntervals.pos),
                                               targetStartPos = covRegStart, targetEndPos = covRegEnd, clipInfo = clipInfo.pos)
        stopifnot( is.vector(features.pos), names(features.pos)[1L] == 'patId' )
        features.pos <- features.pos[-1L]
        
        logging::loginfo("Got %d feature rows for positive training case through SV-injection at covered region no. %d of patient %s. Object has class %s.",
                         NROW(features.pos), j, patId, class(features.pos))
        
        # clean up of tmp BAM files
        mySVBamFiles <- c(injectedBAM[["localBam"]], injectedBAM[["mergedBam"]])
        try(unlink(x = c(mySVBamFiles, paste0(mySVBamFiles, ".bai")), force = TRUE))
        
        
        # prepare list of clipInfos for negative cases
        
        # no peak has been detected above (cf section "assess covered regions") and we take randomly sample a clipped peak proportion (based on what is normal in the WT probe)
        # prop of clipped bases is randomly generated as I have the observed proportion not available here:
        # mkuhn, 2016-01-04: add relS fields using a random value, minimal prop of clipped bases is 1e-4 = 0.0001 (adhoc) and up to what is "normal" (based on the WT mapping, 2*sigma like rule)
        clipInfos.neg <- vector(mode="list", length=nbrNeg)
        relS_random <- 1e-4 + runif(n=2L * length(clipInfos.neg), min = 0L, max = max(0.001, clippedNormalMaxThreshold))
        for (i in seq_along(clipInfos.neg)) clipInfos.neg[[i]] <- c(nbrPeaks=0L, posDiff=0L, relS_5p.99=relS_random[2L*i-1L], relS_3p.99=relS_random[2L*i])
        
        # negative cases (as dataframe)
        logging::loginfo("Features of %d training examples from negative covered region no %d of patient %s", nbrNeg, j, patId)
        features.neg <- getFeaturesFromMappings(patId=patId, bamFile_WT = bamFile_WT, bamFile_MUT = bamFile_MUT,
                                                patMappingInfo = patMappingInfo, chrom = chrom,
                                                regions = myExtSVIntervals.neg, 
                                                targetRegions = IRanges::IRanges(start = rep(covRegStart, nbrNeg), end = covRegEnd ), 
                                                clipInfos = clipInfos.neg)
        logging::loginfo("Got %d feature rows for neg training cases. Object has class %s.", NROW(features.neg), class(features.neg))
        
        
        
        
        ## both a negative case and a positive case through SV-injection
        
        # clipInfo for feature extraction: uses NROW and posDiff
        res.train <- cbind(patId=patId, chrom=chrom, 
                           status=c(rep("IN", nbrNeg), "OUT"),
                           SVtype=c(rep("no", nbrNeg), mySV.types[j]), SVlength = c(rep(0L, nbrNeg), IRanges::width(mySVIntervals.pos)),  SVmargin=margin.bp,
                           rbind(
                             # negative cases (already as dataframe)
                             features.neg[,-which(names(features.neg) %in% c("patId", "chrom"))],
                             features.pos
                           )
        )
        logging::loginfo("Training data for patient %s, region %d with object class %s and dim %s.",
                         patId, j, class(res.train), paste(dim(res.train), collapse = "-"))
        #res.train
        ### DEBUGGING
        featureDat.train[[j]] <- res.train
      } else NULL   # no BAM files
    } else NULL   # covered region no wide enough
  }#foreach
  
  ### DEBUG mkuhn, 2016-01-28
  if ( ! all(sapply(featureDat.train, class) == "data.frame") ){
    logging::logwarn("Feature data of training data for patient %s is for some neg candidate position not a dataframe", patId, class(featureDat.train), NROW(featureDat.train))
    saveRDS(featureDat.train, file = file.path(BASEDIR, "log", paste0("P", patId, "_featureDat_trainError_", format(Sys.Date(), "%Y%m%d"), ".rds")))
  }
  ### DEBUG
  
  # in case of list-combine
  featureDat.train <- do.call(rbind, featureDat.train)
  ### DEBUG
  
  # ensure column order
  stopifnot( all( FEAT_COLS %in% names(featureDat.train)) )
  logging::loginfo("Extracted %d feature rows in training data (with %d incomplete cases) for patient %s", NROW(featureDat.train), sum(! complete.cases(featureDat.train)), patId )
  featureDat.train <- featureDat.train[complete.cases(featureDat.train), c(FEAT_COLS, setdiff(names(featureDat.train), FEAT_COLS))]
  stopifnot( NROW(featureDat.train) >= 2L ) # we have some feature training data after all
  
  featureDat.train$SVtype <- factor(featureDat.train$SVtype, levels = c("no", SV_TYPES), labels = c("no", SV_TYPES))
  
  # mkuhn, 2016-01-11: provide first column 'type' for training data
  featureDat.train <- cbind(type=factor("train", levels = c("train", "test")), featureDat.train)
  
  stopifnot( "type" %in% names(featureDat.train) )
  
  
  
  # build test data -----------------
  
  # positive candidate case (region with clipped base peaks)
  featureDat.test <- foreach(cr.ind=posCandidateRegion.ind, j=seq_along(posCandidateRegion.ind), .combine='rbind') %do% {
    logging::logdebug("Test data for patient %s from positive candidate region j=%d (candReg.ind= %d)", patId, j, cr.ind)
    
    # covered region start and end
    regStart <- IRanges::start(coveredRegions[cr.ind])
    regEnd   <- IRanges::end(coveredRegions[cr.ind])
    
    mySVCandStart <- clippedRes_MUT[cr.ind, "pos.left"]
    mySVCandEnd   <- clippedRes_MUT[cr.ind, "pos.right"]
    
    stopifnot( ! is.na(mySVCandStart), ! is.na(mySVCandEnd) )
    
    mySVCandPosition <- IRanges::IRanges(start = pmax(regStart, mySVCandStart - margin.bp),
                                         end = pmin(regEnd, mySVCandEnd + margin.bp))
    
    if ( IRanges::start(mySVCandPosition) == regStart || IRanges::end(mySVCandPosition) == regEnd ){
      logging::logwarn("Positive SV-candidate position %d has extended range (w/ margin) on region border! Effective margins are %d and %d", cr.ind,
                       mySVCandStart - IRanges::start(mySVCandPosition), IRanges::end(mySVCandPosition) - mySVCandEnd)
    }
    
    mySVCandClipInfo <- clippedRes_MUT[cr.ind, c("nbrPeaks", "posDiff", "relS_5p.99", "relS_3p.99")]
    ##stopifnot( is.vector(clipInfo), all( c("nbrPeaks", "posDiff", "relS_5p.99", "relS_3p.99") %in% names(clipInfo)) )
    
    # output as dataframe
    cbind(
      data.frame(patId=patId, chrom=chrom), SVmargin=margin.bp,
      t(getFeaturesFromMapping(patId=patId, bamFile_WT = bamFile_WT, bamFile_MUT = bamFile_MUT, patMappingInfo = patMappingInfo, 
                               chrom = chrom, startPos = IRanges::start(mySVCandPosition), endPos = IRanges::end(mySVCandPosition),
                               targetStartPos = regStart, targetEndPos = regEnd, clipInfo = mySVCandClipInfo)[-1L])
    )
  }
  
  logging::loginfo("Extracted %d feature rows for test data for patient %s", NROW(featureDat.test), patId )
  # column order for base columns that describe patient and position. featureDat.test can be empty
  stopifnot( NROW(featureDat.test) == 0L || all( FEAT_COLS0 %in% names(featureDat.test)) )
  featureDat.test <- featureDat.test[, c(FEAT_COLS0, setdiff(names(featureDat.test), FEAT_COLS0))]
  
  
  list(train=featureDat.train, test=featureDat.test)
}




#' Add the outcome variable to a feature test dataset based on simulated truth.
#' @param evalMargin numeric. How to widen the simulated SVs prior to match with found test SV-intervals 
#' @param regionIsMatch logical. Do we take the simulated truth of the whole enriched target area? (i.e. do not look into eval-margin then)
#' @param featNames feature names from training data
#' @return dataframe with features of simulation test data with updated status (like the training data)
#' @export
addStatusToTestFromSim <- function(featDat.test, patData, myCHROM, evalMargin, regionIsMatch=FALSE, featNames) {
  
  if (FALSE){
    myCHROM <- "chr5"
    evalMargin <- 25L #50L #
    patData <- readRDS(file.path(getVirtualPatientPath(what = "patient"), "patData.rds"))
    regionIsMatch <- FALSE
    
    # start from existing feature file
    allFeat <- readRDS("~/SVdata/feature/cov60_TL90/PEm200sd20_R100/features_sim_screen_margin25.rds")
    featNames <- names(allFeat)
    featDat.test <- allFeat %>%
      dplyr::filter_(~type == 'test') %>% 
      dplyr::select_(~-type, ~-status, ~-SVtype, ~-SVlength)
  }
  
  # restrict to chrom
  featDat.test <- featDat.test[featDat.test$chrom == myCHROM,]
  patData <- patData[patData$chrom == myCHROM,]
  
  
  if (NROW(featDat.test) == 0L){
    logging::loginfo("No test data found to add status column!")
    return(invisible(NULL))
  }
  
  # drop stale status columns: I will recreate them
  featDat.test %<>% dplyr::select_(~ -one_of("type" ,"status", "SVtype", "SVlength"))
  
  # look into patient data: what was simulated
  # mkuhn, 2016-02-09: keep only the simulated SVs as relevant patData entries
  stopifnot( identical( c("DSV", "germline"), unique(patData$status)) )
  relPatData <- patData %>%
    dplyr::mutate_(simStartPos=~targetStartPos + SVstart, 
                   simEndPos=~simStartPos + SVlength) %>% 
    dplyr::rename_(simTargetStartPos=~targetStartPos, simTargetEndPos=~targetEndPos) %>% 
    dplyr::filter_(~status == 'DSV') %>% 
    dplyr::select_(~patId, ~chrom, ~simTargetStartPos, ~simTargetEndPos, ~simStartPos, ~simEndPos, ~SVtype, ~SVlength)
  
  
  
  # patient loop -----
  
  # for each patient set the status variable in the feature dataframe
  # if an overlap is detected with a simulated SV (as recorded in patData) than write out the corresponding status
  # if no overlap is found than it is a "no SV"
  updatedFeatDat.test <- foreach::foreach(myPatId=unique(featDat.test$patId), .combine = 'rbind', .multicombine = FALSE, .errorhandling = "stop") %do% {
    
    fDat.test <- featDat.test %>%
      dplyr::filter_(~patId == myPatId)
    myRelPatData <- relPatData %>%
      dplyr::filter_(~patId == myPatId)
    
    logging::logdebug("From patient data, select %d rows for patient %d.", NROW(myRelPatData), myPatId)
    
    # this is the simulated truth (with added evaluation margin or using whole target region)
    # no SVs are filtered out, they would have simStartPos = simEndPos within their target region.
    simSVRanges <- with(myRelPatData,
                        GenomicRanges::GRanges(seqnames = myCHROM,
                                               ranges = if (isTRUE(regionIsMatch))
                                                 IRanges::IRanges(start = simTargetStartPos,
                                                                  end = simTargetEndPos) else
                                                                    IRanges::IRanges(start = simStartPos,
                                                                                     end = simEndPos)))
    # evalMargin moved to maxgap = , see below!
    # IRanges::IRanges(start = simStartPos - evalMargin,
    #                  end = simEndPos + evalMargin)))
    
    # this is what our screening found
    testSVRanges <- with(fDat.test,
                         GenomicRanges::GRanges(seqnames = myCHROM,
                                                ranges = IRanges::IRanges(start = startPos,
                                                                          end = endPos)))
    
    #GenomicRanges::subsetByOverlaps # subset of query-GRanges
    #GenomicRanges::overlapsAny # logical vector of length equal to query indicating those that overlap any of the subject-ranges
    #GenomicRanges::findOverlaps # returns a Hit-object (for select="all") or an integer vector of length query with a matching index in subject
    # patDataInd for each query-range (test-regions from screening) gives a matching index into the subject (or NA if no match)
    patDataInd <- GenomicRanges::findOverlaps(query = testSVRanges, subject = simSVRanges, select = "arbitrary", maxgap = evalMargin)
    
    # extract the simulated SV-type and SV-length, in the order of the feature data hits
    patDataSVtype <- myRelPatData$SVtype[patDataInd]
    patDataSVtype[is.na(patDataSVtype)] <- "no"
    
    patDataSVlen <- myRelPatData$SVlength[patDataInd]
    patDataSVlen[is.na(patDataSVlen)] <- 0L
    
    # check, double-check
    stopifnot( identical(STATUS_LEVELS, c("IN", "OUT")) )
    
    fDat.test.ext <- cbind(type=factor("test", levels = c("train", "test")),
                           status=factor(patDataSVtype != 'no', levels = c(FALSE, TRUE), labels = STATUS_LEVELS),
                           SVtype=patDataSVtype,
                           SVlength=patDataSVlen,
                           fDat.test)
    
    if (! setequal(featNames, names(fDat.test.ext)) ){
      logging::logerror("Names of training and test data are not equal at patId %d.Skipping the test data!\nNames of training data:\n%s\nNames of test data (%d rows):\n%s",
                        myPatId, paste(featNames, collapse=" - "), NROW(fDat.test.ext), paste(names(fDat.test.ext), collapse=" - "))
      return(invisible(NULL))
    } 
    
    stopifnot( setequal(featNames, names(fDat.test.ext)) )
    
    
    # column ordering
    fDat.test.ext[featNames]
    
  }#hcaerof
  
  # fix factor levels for SVtype
  updatedFeatDat.test$SVtype <- factor(updatedFeatDat.test$SVtype, levels = c("no", SV_TYPES), labels = c("no", SV_TYPES))
  updatedFeatDat.test
}



#' Add missing columns for test data of a patient in the Biotec data, based on confirmed FLT3-ITD hits.
#' 
#' Only FLT3-ITD is marked here. Other hits are not known and are not recorded.
#' @param regionIsMatch flag if whole target region counts as match or if overlap of SV +/- evalMargin is required
#' @param featNames column names vector of training data
#' @export
addStatusForFLT3ITD <- function(featDat.test, evalMargin, regionIsMatch=FALSE, featNames) {
  
  flt3CHROM <- "chr13"
  # use consensus region for all four FLT3-ITD patients from manuscript
  FLT3ITD <- GenomicRanges::GRanges(seqnames = flt3CHROM,
                                    ranges = IRanges::IRanges(start = 28608242, end = 28608280 + 66L))
  
  # FLT3-ITD target region (TruSeq)
  FLT3ITD_TARGETREGION <- GenomicRanges::GRanges(seqnames = flt3CHROM,
                                                 #ranges = IRanges::IRanges(start = 28575411, end=28676729)) # whole FLT3 region
                                                 ranges = IRanges::IRanges(start = 28607945, end=28608622)) # target region with ITD
  
  
  # restrict to chrom chr13
  featDat.test <- featDat.test[featDat.test$chrom == flt3CHROM,]
  
  if (NROW(featDat.test) == 0L){
    logging::loginfo("No test data on chrom %s found!", flt3CHROM)
    return(invisible(NULL))
  }
  
  # drop stale status columns: I will recreate them
  featDat.test %<>% dplyr::select_(~ -one_of("type" ,"status", "SVtype", "SVlength"))
  
  
  # this is the simulated truth (with added evaluation margin or using whole target region)
  FLT3ITD_EXT <- if (isTRUE(regionIsMatch)) FLT3ITD_TARGETREGION else
    GenomicRanges::resize(FLT3ITD, width = IRanges::width(FLT3ITD) + 2L * evalMargin, fix = "center")
  
  
  # this is what our screening found
  testSVRanges <- with(featDat.test,
                       GenomicRanges::GRanges(seqnames = flt3CHROM,
                                              ranges = IRanges::IRanges(start = startPos,
                                                                        end = endPos)))
  
#   SV_HitsInTest <- GenomicRanges::findOverlaps(query = testSVRanges, subject = FLT3ITD_EXT)
#   # for each query-range (test-regions from screening) gives a matching index into the subject (or NA if no match)
#   patDataInd <- S4Vectors::selectHits(SV_HitsInTest, select = "arbitrary") # any hit from simulated SV-ranges (should only be one!)
  testSVFLT3ITDMatch <- GenomicRanges::overlapsAny(query = testSVRanges, subject = FLT3ITD_EXT)
  
  
  featDat.test.ext <- cbind(type=factor("test", levels = c("train", "test")),
                         status=factor(testSVFLT3ITDMatch, levels = c(FALSE, TRUE), labels = STATUS_LEVELS),
                         SVtype=factor(ifelse(testSVFLT3ITDMatch, "duplication", "no"), levels = c("no", SV_TYPES), labels = c("no", SV_TYPES)),
                         SVlength=ifelse(testSVFLT3ITDMatch, yes = IRanges::width(FLT3ITD), no = 0L),
                         featDat.test)
  
  if (! setequal(featNames, names(featDat.test.ext)) ){
    logging::logerror("Names of training and test data are not equal at patId %d.Skipping the test data!\nNames of training data:\n%s\nNames of test data (%d rows):\n%s",
                      myPatId, paste(featNames, collapse=" - "), NROW(featDat.test.ext), paste(names(featDat.test.ext), collapse=" - "))
    return(invisible(NULL))
  } 
  
  stopifnot( setequal(featNames, names(featDat.test.ext)) )
  # column ordering
  featDat.test.ext <- featDat.test.ext[featNames]
  
  featDat.test.ext
}
