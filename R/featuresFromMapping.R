

#' Counts primary & non-supplementory read ends.
#' @param minMapQ non-negative integer. Minimal mapping quality as stated in the mapping
#' @param flagsSet integer values. sum of bitmask values which must be set (-f). Default is 0, means no flag filter.
#' @param flagsNotSet integer value. sum of bitmask values which must not be set (-F). Defaults is not secondary and not supplementary reads.
#' @return number of reads in the bamFile
#' @note Alas, \code{countBam} does not handle secondary reads..
#' @export
getReadCountsFromMapping <- function(bamFile, flagsSet=0L, flagsNotSet=FLAG_SECONDARY + FLAG_SUPPL, minMapQ=0L, regionStr=NULL){
  stopifnot( ! is.null(bamFile) )
  if (is(bamFile, "BamFile")){
    bamFile <- Rsamtools::path(bamFile)
  }
  stopifnot( checkFile(bamFile) )
  stopifnot( length(minMapQ) == 1L, is.numeric(minMapQ) )
  stopifnot( is.null(regionStr) || is.character(regionStr) && length(regionStr) == 1L )
  if (is.null(flagsSet)) flagsSet <- 0L
  if (is.null(flagsNotSet)) flagsNotSet <- 0L
  
  stopifnot(is.numeric(flagsSet), flagsSet >= 0L, is.numeric(flagsNotSet), flagsNotSet >= 0L )
  stopifnot( bitops::bitAnd(flagsSet, flagsNotSet) == 0L )
  
  minMapQ <- max(0L, ceiling(minMapQ))
  
  # Rsamtools::countBam?
  bamCountCall <- system(paste(SAMTOOLS_EXE, "view -c ",
                               if (minMapQ > 0L) paste0("-q", minMapQ),
                               if (flagsSet > 0L) paste0("-f", flagsSet),
                               if (flagsNotSet > 0L) paste0("-F", flagsNotSet),
                               bamFile,
                               if (! is.null(regionStr) && nzchar(regionStr)) regionStr),
                         intern=TRUE)
  
  as.numeric(bamCountCall)
  
}

# #' Read out the BAM header
# getBAMHeader <- function(bamFile){
#   Rsamtools::scanBamHeader(bamFile)
# }




#' Estimate NGS properties from a given mapping
#' 
#' Median read length and mean and sd of insert size distribution are estimated through sampling in mapping.
#' @param bamFile character to BAM file.
#' @return NGS property object 
estimateNGSPropFromMapping <- function(bamFile, wBioconductor=TRUE){
  
  ret <- ngs.prop(read.length=100L, pe.ins.mean=185, pe.ins.sd = 28)
  
  if (isTRUE(wBioconductor)){
    #Rsamtools::quickBamFlagSummary(bamFile)
    bf <- Rsamtools::BamFile(bamFile, yieldSize=10000)
    yield <- function(x) GenomicAlignments::readGAlignmentPairs(x, param=Rsamtools::ScanBamParam(flag = Rsamtools::scanBamFlag(isPaired=TRUE, isProperPair = TRUE, isSecondaryAlignment = FALSE),
                                                                                                 simpleCigar = TRUE,
                                                                                                 what=c( "isize", "mapq", "qwidth" ))) %>% 
      suppressWarnings
    # do.call(cbind, Rsamtools::scanBam(x, param=Rsamtools::ScanBamParam(flag = scanBamFlag(isPaired=TRUE, isProperPair = TRUE, isSecondaryAlignment = FALSE),
    #                                                     simpleCigar = TRUE,
    #                                                     what=c( "flag", "isize", "mapq" )))[[1L]][c("isize", "mapq")])
    map <- function(x) { md <- S4Vectors::mcols(GenomicAlignments::first(x)); list(isize=md$isize[md$mapq > 5L & md$isize > 0L], qwidth=md$qwidth) } # function(x) x[, "isize"][x[, "mapq"]>5L]
    ## Samples records from successive chunks of aligned reads.
    readStats <- GenomicFiles::reduceByYield(bf, yield, map, REDUCE = function(l1, l2) list(mean(c(l1[["isize"]], l2[["isize"]]), na.rm = TRUE),
                                                                                            mean(c(sd(l1[["isize"]], na.rm = TRUE), sd(l2[["isize"]]), na.rm = TRUE), na.rm = TRUE),
                                                                                            mean(c(l1[["qwidth"]], l2[["qwidth"]]), na.rm = TRUE) ))
    #GenomicFiles::REDUCEsampler(1e3L, TRUE), iterate=TRUE)
    names(readStats) <- c("isize.m", "isize.sd", "qwidth.m")
    ## mkuhn, 2015-12-22: close file handle?! --apparently, he does not like it
    ##try(close(bf))
    
    # ret <- ngs.prop(read.length = median(readStats$qwidth, na.rm=TRUE),
    #                 pe.ins.mean=round(mean(readStats$isize, na.rm=TRUE)),
    #                 pe.ins.sd=round(sd(readStats$isize, na.rm=TRUE)))
    ret <- ngs.prop(read.length = ceiling(readStats$qwidth.m),
                    pe.ins.mean = ceiling(readStats$isize.m),
                    pe.ins.sd = 1L+ceiling(readStats$isize.sd))
  } 
  
  ret
  
}


#' Patient level information from her tumor-normal mapping.
#' 
#' Added an estimate of clipped base proportion based on WT-sample. It estimates a mean clipped proportion and a measure of "non-peak"-deviation by discarding the highest 2.5% of values (avoid clipped peaks).
#' @param YIELD_SIZE size of chunk used to go through the mapping files
#' @param chrom character. chromosome
#' @param flag to use bioconductor software. 
#' @return list with different bits of information from the mapping: read count ratio between MUT to WT, clipped base info for WT (yields clipped prop thresholds), secondary and suppl reads 
#' @export
getPatientMappingInfo <- function(bamFile_WT, bamFile_MUT, chrom=CHROM, YIELD_SIZE=5e5, wBioconductor=TRUE) {
  
  
  # read count ratio -------
  
  ### READ count of full BAM file (reflects the difference in overall read count betw. tumor and normal)
  rc.bam.wt <- getReadCountsFromMapping(bamFile = bamFile_WT,  flagsNotSet = FLAG_SECONDARY + FLAG_SUPPL, minMapQ = 0L, regionStr = chrom)
  rc.bam.mut <- getReadCountsFromMapping(bamFile = bamFile_MUT, flagsNotSet = FLAG_SECONDARY + FLAG_SUPPL, minMapQ = 0L, regionStr = chrom)
  
  # mkuhn, 2015-03-03: counting only q1-reads to have less artefact records in the count
  # but very similar to primary non-suppl reads
  # rc.q1.bam.wt <- getReadCountsFromMapping(bamFile = bamFile_WT, minMapQ = 1L)
  # rc.q1.bam.mut <- getReadCountsFromMapping(bamFile = bamFile_MUT, minMapQ = 1L)
  
  
  # [mkuhn, 2015-03-03] rc.bam.ratio factor for normalization of coverage numbers
  # we have taken out secondary and supplementary reads. duplicated reads should be dropped already, an issue here? Q0-reads?
  rc.bam.ratio <-  rc.bam.mut / rc.bam.wt
  
  
  
  
  # clipped bases of WT -----
  
  bf_wt <- Rsamtools::BamFile(bamFile_WT, yieldSize=YIELD_SIZE)
  clippedStats <- if (isTRUE(wBioconductor)){
    # mkuhn, 2015-12-22: proportion of clipped bases in WT
    
    yield <- function(x) GenomicAlignments::readGAlignments(x, param=Rsamtools::ScanBamParam(flag = Rsamtools::scanBamFlag(isSecondaryAlignment = FALSE))) #c("qwidth", "cigar")) )
    map <- function(x) {
      cig <- GenomicAlignments::cigar(x);
      startClipped <- cig %>% stringr::str_extract(pattern = "^[[:digit:]]+S") %>% stringr::str_replace(pattern = "S", replacement = "") %>% as.numeric 
      startClipped[is.na(startClipped)] <- 0L
      endClipped <- cig %>% stringr::str_extract(pattern = "[[:digit:]]+S$") %>% stringr::str_replace(pattern = "S", replacement = "") %>% as.numeric
      endClipped[is.na(endClipped)] <- 0L
      # sort starts with smallest values
      clippedProp <- sort((startClipped + endClipped) / GenomicAlignments::qwidth(x))
      #list(m=weighted.mean(x=clippedProp, w = GenomicAlignments::qwidth(x)), s=sd(clippedProp))
      #weighted.mean(x=clippedProp, w = GenomicAlignments::qwidth(x))
      c( clippedBases=sum(startClipped + endClipped, na.rm = TRUE),
         nbrBases=sum(GenomicAlignments::qwidth(x), na.rm = TRUE),
         sdClippedProp=sd(clippedProp[1L:(min(length(clippedProp), 2L+ceiling(length(clippedProp)) * .975))], na.rm = TRUE) )
      #quantile(x=(startClipped + endClipped) / GenomicAlignments::qwidth(x), probs = .95)
    }
    ## Samples records from successive chunks of aligned reads.
    res <- GenomicFiles::reduceByYield(bf_wt, YIELD=yield, MAP = map,
                                       REDUCE=function(v1, v2) {
                                         vec3 <- c(v1[3L], v2[3L]); vec3 <- vec3[!is.na(vec3) & ! is.null(vec3)]
                                         c(clippedBases = sum(v1[1L],v2[1L], na.rm = TRUE),
                                           nbrBases = sum(v1[2L], v2[2L], na.rm = TRUE),
                                           #sdClippedProp = mean(v1[3L], v2[3L], na.rm = TRUE)),
                                           # #weighting through iterations, early ones get diluted by taking iterated mean, hence first entry double weight
                                           sdClippedProp = if (length(vec3) == 2L) weighted.mean(x=vec3, w=c(2L,1L), na.rm = TRUE) else mean(x=vec3, na.rm=TRUE))
                                       }, iterate = TRUE)
    #`+`, iterate=TRUE) #function(x) mean(x, trim = 0, na.rm = TRUE),
    
    list(m = (res["clippedBases"] + 1L) / (res["nbrBases"] + 2L), s = res["sdClippedProp"] + 0.005)
    # REDUCEsampler(nmbr) # number must be smaller than yield size?!
    
    # mkuhn, 2015-12-22: close file handle?!
    #try(close(bf_wt))
  }
  # else { # default values for w/o bioconductor
  #   list(m = 0.03, s = 0.004)
  # }
  
  
  # threshold that we require for a (smoothed) peak to exceed
  # at max 50%, at least 5%
  clippedPropPeakMinThreshold <- max(.05, min(.5, clippedStats[["m"]] + 6L * clippedStats[["s"]]))
  clippedPropNormalMaxThreshold <- min(clippedPropPeakMinThreshold / 2L, clippedStats[["m"]] + 2.33 * clippedStats[["s"]])  #qnorm(.99)
  
  
  
  
  # Secondary mapping -----
  
  # reads with 2nd alignment anywhere in the mapping 
  myParam_2nd <- Rsamtools::ScanBamParam(flag=Rsamtools::scanBamFlag(isPaired=NA, isProperPair=NA, isSecondaryAlignment = TRUE),
                                         what=c("qname"))
  #Rsamtools::scanBam(file=bamFile_WT,  param=myParam_2nd) 
  reads_2nd_WT  <- unique(Rsamtools::scanBam(file=bamFile_WT,  param=myParam_2nd)[[1]]$qname)
  reads_2nd_MUT <- unique(Rsamtools::scanBam(file=bamFile_MUT, param=myParam_2nd)[[1]]$qname)
  
  
  
  
  # Supplementary reads -----
  
  # mkuhn, 2015-12-22: supplementary reads anywhere in the mapping. We use the windowing of BamFile because supplementary reads 
  # are not supported in Bioconductor's Rsamtools (sadly)
  bf_mut <- Rsamtools::BamFile(bamFile_MUT, yieldSize=YIELD_SIZE)
  yieldF <- function(x) GenomicAlignments::readGAlignments(x, param=Rsamtools::ScanBamParam(flag = Rsamtools::scanBamFlag(isSecondaryAlignment = NA),
                                                                                            what = c("qname", "flag")) )#c("qwidth", "cigar")) )
  mapF <- function(ga) {
    flag <- S4Vectors::mcols(ga)$flag
    
    S4Vectors::mcols(ga)$qname[bitops::bitAnd(flag, FLAG_SUPPL) > 0L]
    # list(
    #   supplReads=mcols(ga)$qname[bitops::bitAnd(flag, FLAG_SUPPL) > 0L],
    #   secondaryReads=mcols(ga)$qname[bitops::bitAnd(flag, FLAG_SECONDARY) > 0L]
    # )
    #cat("\nret:",ret)
  }
  
  ## Samples records from successive chunks of aligned reads.
  supplReads_WT  <- unique(GenomicFiles::reduceByYield(bf_wt, YIELD=yieldF,  MAP = mapF, REDUCE = c, iterate = TRUE))
  supplReads_MUT <- unique(GenomicFiles::reduceByYield(bf_mut, YIELD=yieldF, MAP = mapF, REDUCE = c, iterate = TRUE))
  
  
  list(clippedBases=list(mean=clippedStats[["m"]], sd=clippedStats[["s"]], peakMinThreshold=clippedPropPeakMinThreshold, normalMaxThreshold=clippedPropNormalMaxThreshold),
       rc.bam.ratio = rc.bam.ratio,
       wt=list(secReads=reads_2nd_WT, supplReads=supplReads_WT), mut=list(secReads=reads_2nd_MUT, supplReads=supplReads_MUT) )
}




#' Extracts features from  tumor/normal mapping at \emph{one} selected locus.
#' 
#' Method to extract features from mapping characteristics of tumor/normal sample at single selected locus.
#' It uses baseline characteristcs from the normal sample and difference values defined as tumor - normal
#' (i.e. diff is how much you need to add to come from normal to tumor).
#' Only the local mapping (at the target region) is of interest here, global mapping information is given separately.
#' 
#' Both mapping files should be based on the same part of the reference genome.
#' The mappings may stem from simulated reads at a certain location. 
#' The genomic region describes the locus from where to extract features.
#' In simulations this locus is known, in real-world application these come from a pre-screening.
#' 
#' @param patId Patient ID (for simulation as numeric, e.g. 17)
#' @param bamFile_WT filename for BAM-file of normal sample. [mkuhn, 2015-12-15: will filter around target region]
#' @param bamFile_MUT filename for BAM-file of tumor sample  [mkuhn, 2015-12-15: will filter around target region]
#' @param patMappingInfo list with information on whole mapping for tumor/normal of patient [mkuhn, 2015-12-15: new to allow to have local BAMs for the rest]
#' @param chrom chromosome of single locus
#' @param startPos single integer value. start position of selected candidate region
#' @param endPos single integer value. end position of selected candidate region
#' @param targetStartPos single integer value. start position of target region that contains the selected region
#' @param targetEndPos single integer value. end position of target region that contains the selected region
#' @param clipInfo vector with clipping info for the selected region
#' @return Returns a feature vector that describes certain mapping characteristics that might help to differentiate tumor from normal sample.
#' @export
getFeaturesFromMapping <- function (patId, bamFile_WT, bamFile_MUT, patMappingInfo,
                                    chrom=CHROM, startPos, endPos, #startPos=170836544L, endPos=170838544L,
                                    targetStartPos=startPos, targetEndPos=endPos, clipInfo) {
  if (FALSE){
    # example patient 1000 with inversion locus
    patId <- 3L;
    bamFile_WT  <- paste0(BASEDIR, "mapping/cov60_TL90/PEm200sd25_R100/", getPatStr(patId), "_N.s.bam")
    bamFile_MUT <- paste0(BASEDIR, "mapping/cov60_TL90/PEm200sd25_R100/", getPatStr(patId), "_T.s.bam")
    patMappingInfo <- getPatientMappingInfo(bamFile_WT = bamFile_WT, bamFile_MUT = bamFile_MUT)
    chrom <- "chr5"; startPos <- 36974877 + 646 - 25; endPos <- 36974877 + 646 + 259 + 25; targetStartPos <- 36974877; targetEndPos <- 36977504;
    #clipInfo <- ???
  }
  
  stopifnot( checkFile(SAMTOOLS_EXE) )
  stopifnot( checkFile(bamFile_WT), checkFile(bamFile_MUT) )
  stopifnot( targetStartPos <= startPos, targetEndPos >= endPos )
  stopifnot( is.vector(clipInfo), all( c("nbrPeaks", "posDiff", "relS_5p.99", "relS_3p.99") %in% names(clipInfo)) )
  
  # BAM coordinates are refering to left-most mapped position, i.e. this coordinate stands for a read that has its center to the right of it.
  startPos_BAM <- round(startPos - READL/2L);  endPos_BAM <- round(endPos-READL/2L)
  targetStartPos_BAM <- round(targetStartPos - READL/2L);  targetEndPos_BAM <- round(targetEndPos - READL/2L)
  
  SELECTED_REGION_GEN <- paste0(chrom, ":", startPos, "-", endPos)
  SELECTED_REGION_BAM <- paste0(chrom, ":", startPos_BAM, "-", endPos_BAM) 
  TARGET_REGION_GEN <- paste0(chrom, ":", targetStartPos, "-", targetEndPos)
  TARGET_REGION_BAM <- paste0(chrom, ":", targetStartPos_BAM, "-", targetEndPos_BAM) 
  
  targetRange <- GenomicRanges::GRanges(seqnames = chrom,
                                        ranges = IRanges::IRanges(start = targetStartPos, end = targetEndPos))
  regionRange <- GenomicRanges::GRanges(seqnames = chrom,
                                        ranges=IRanges::IRanges(start=startPos, end=endPos))
  
  # extended region: this is used to check for reads with unmapped mates
  regionExtension <- 6L * READL #adhoc
  regionExtRange <- GenomicRanges::GRanges(seqnames = chrom,
                                           ranges=IRanges::IRanges(start=startPos-regionExtension, end=endPos+regionExtension))
  
  idInfo <- c(patId=patId, startPos=startPos, endPos=endPos, targetStartPos=targetStartPos, targetEndPos=targetEndPos)
  
  bamRegStr <- paste0("_Pat", patId, "_", TARGET_REGION_GEN)
  
  
  # mkuhn, 2015-12-15: I do not use filtered BAM files because I hope that with index we can work on whole BAM file. Otherwise, I generate lots of file operations..
  # # # mapping filtered down to target region
  # bamFileTarget_WT <- Rsamtools::filterBam(file = bamFile_WT, destination = file.path(TMPDIR, paste0(basename(bamFile_WT), bamRegStr)),
  #                                       indexDestination = TRUE, param=Rsamtools::ScanBamParam(which = targetRange))
  # 
  # bamFileTarget_MUT <- Rsamtools::filterBam(file = bamFile_MUT, destination = file.path(TMPDIR, paste0(basename(bamFile_MUT), bamRegStr)),
  #                                       indexDestination = TRUE, param=Rsamtools::ScanBamParam(which = targetRegRange))
  # 
  # on.exit(try(unlink(x=c(bamFileTarget_WT, bamFileTarget_MUT))))
  
  
  
  
  # REGION features -------------
  
  # length of the selected region
  feat.region.length.r <- endPos - startPos + 1L
  
  
  # minimal distance from selected region to target region borders,
  # ZZZ Should I have a maximum value here as upper bound?
  # ZZZ feature name: let it end with .t (and it will be picked up!)
  feat.region.targetDist.5p.r <- log(abs(startPos-targetStartPos)+1L) #min( abs(startPos-targetStartPos), 5*READL )
  feat.region.targetDist.3p.r <- log(abs(targetEndPos - endPos)+1L) #min( abs(targetEndPos - endPos), 5*READL ) 
  
  
  
  # read counts features -----
  
  rc.bam.ratio <- patMappingInfo$rc.bam.ratio
  
  # read entries in the target region
  rc.target.wt  <- getReadCountsFromMapping(bamFile = bamFile_WT,  flagsNotSet = FLAG_SECONDARY + FLAG_SUPPL, regionStr = TARGET_REGION_BAM)
  rc.target.mut <- getReadCountsFromMapping(bamFile = bamFile_MUT, flagsNotSet = FLAG_SECONDARY + FLAG_SUPPL, regionStr = TARGET_REGION_BAM)
  
  # read entries in the selected region
  rc.reg.wt  <- getReadCountsFromMapping(bamFile = bamFile_WT,  flagsNotSet = FLAG_SECONDARY + FLAG_SUPPL, regionStr = SELECTED_REGION_BAM)
  rc.reg.mut <- getReadCountsFromMapping(bamFile = bamFile_MUT, flagsNotSet = FLAG_SECONDARY + FLAG_SUPPL, regionStr = SELECTED_REGION_BAM)
  
  # secondary reads (alternative mapping positions e.g. due to repeats) in the selected region
  rc.reg.2nd.wt  <- getReadCountsFromMapping(bamFile = bamFile_WT, flagsSet = FLAG_SECONDARY, flagsNotSet = 0L, regionStr = SELECTED_REGION_BAM)
  rc.reg.2nd.mut <- getReadCountsFromMapping(bamFile = bamFile_MUT, flagsSet = FLAG_SECONDARY, flagsNotSet = 0L, regionStr = SELECTED_REGION_BAM)
  
  # supplementary reads (chimerics e.g. due to SV) in the selected region
  rc.reg.suppl.wt  <- getReadCountsFromMapping(bamFile = bamFile_WT, flagsSet = FLAG_SUPPL, flagsNotSet = 0L, regionStr = SELECTED_REGION_BAM)
  rc.reg.suppl.mut <- getReadCountsFromMapping(bamFile = bamFile_MUT, flagsSet = FLAG_SUPPL, flagsNotSet = 0L, regionStr = SELECTED_REGION_BAM)
  
  # mkuhn, 2015-03-03: commented out - as target region is not really biologically meaningful?! Use only selected region.
  #   # percentage of reads on target region
  #   feat.target.readProp.1 <- rc.target.wt / (rc.bam.wt+1L) * 100L
  #   feat.target.readProp.d <- (rc.target.mut / (rc.bam.mut+1L) * 100L) - feat.target.readProp.1
  
  # percent of reads on selected region relative to target region [mkuhn, 20150419: used to be relative to all reads in BAM file before]
  feat.readProp.1 <- rc.reg.wt / (rc.target.wt+1L) * 100L  #used to be in the denominator: (rc.bam.wt+1L) * 100L
  feat.readProp.d <- (rc.reg.mut / (rc.target.mut+1L) * 100L) - feat.readProp.1
  
  # percent of *secondary* read segments on region relative to all read segments in region (multiple mappings e.g. due to repeats in ref)
  feat.readProp.2nd.1 <- rc.reg.2nd.wt / (rc.reg.wt+1L) * 100L
  feat.readProp.2nd.d <- (rc.reg.2nd.mut / (rc.reg.mut+1L) * 100L) - feat.readProp.2nd.1
  
  # percent of *supplementary* read segments on region relative to all read segments in region (chimeric alignment, parts of a read map somewhere else e.g. due to SVs)
  feat.readProp.suppl.1 <- rc.reg.suppl.wt / (rc.reg.wt+1L) * 100L
  feat.readProp.suppl.d <- (rc.reg.suppl.mut / (rc.reg.mut+1L) * 100L) - feat.readProp.suppl.1
  
  
  # all reads segments (properly paired or not) that map into the selected genomic region of interest
  bamParam.reg <- Rsamtools::ScanBamParam(flag=Rsamtools::scanBamFlag(isPaired=NA, isProperPair=NA), #, isSecondaryAlignment=FALSE),
                                          what=setdiff(Rsamtools::scanBamWhat(), c("seq", "qual", "qwidth")), 
                                          which=regionRange) # genomic coordinates are used by ScanBamParam()
  mapping.reg.N <- data.frame(Rsamtools::scanBam(file=bamFile_WT, param=bamParam.reg)[[1]])
  mapping.reg.T <- data.frame(Rsamtools::scanBam(file=bamFile_MUT, param=bamParam.reg)[[1]])
  
  
  # Only regular mapping of reads in the selected region of interest (because we want to have a single entry per read end)
  mapping.reg.regular.N <- mapping.reg.N[which(bitops::bitAnd(mapping.reg.N$flag, FLAG_SECONDARY + FLAG_SUPPL) == 0L),]
  mapping.reg.regular.T <- mapping.reg.T[which(bitops::bitAnd(mapping.reg.T$flag, FLAG_SECONDARY + FLAG_SUPPL) == 0L),]  
  
  
  # mkuhn, 2016-01-03: target region as wider window around selected region
  bamParam.regionExt <- Rsamtools::ScanBamParam(flag=Rsamtools::scanBamFlag(isPaired=NA, isProperPair=NA), #, isSecondaryAlignment=FALSE),
                                                what=setdiff(Rsamtools::scanBamWhat(), c("seq", "qual", "qwidth")), 
                                                which=regionExtRange) # genomic coordinates are used by ScanBamParam()
  mapping.regionExt.N <- data.frame(Rsamtools::scanBam(file=bamFile_WT, param=bamParam.regionExt)[[1]])
  mapping.regionExt.T <- data.frame(Rsamtools::scanBam(file=bamFile_MUT, param=bamParam.regionExt)[[1]])
  
  
  # proportion of discordant pairs (orientation either F-F or R-R) [mkuhn, 2015-03-03]
  feat.readProp.discordant.1 <- length(which(bitops::bitAnd(mapping.reg.regular.N$flag, FLAG_READ_REV + FLAG_MATE_REV) %in% c(0L, FLAG_READ_REV + FLAG_MATE_REV))) / NROW(mapping.reg.regular.N)
  feat.readProp.discordant.d <- length(which(bitops::bitAnd(mapping.reg.regular.T$flag, FLAG_READ_REV + FLAG_MATE_REV) %in% c(0L, FLAG_READ_REV + FLAG_MATE_REV))) / NROW(mapping.reg.regular.T) - feat.readProp.discordant.1
  
  
  # also consider 2nd or supplementary reads that map outside the region but have like (are paired) to the region
  # all templates anywhere (identified by QNAME) that have secondary (2nd) alignment with a link to the region via their read-pair: 
  # they either map into that region secondarily or their primary mapping is there but they have 2nd mappings.  2nd mappings can be pretty far away from primary mapping place.
  secondaryAlignmentsFromRegion.N <- intersect(patMappingInfo[["wt"]][["secReads"]], mapping.reg.N$qname)
  secondaryAlignmentsFromRegion.T <- intersect(patMappingInfo[["mut"]][["secReads"]], mapping.reg.T$qname)
  
  # prop of secondary (i.e. non-primary) mappings relative to primary reads on selected region:
  # theoretically this number can go beyond 1 [if all primary read pairs in region have a 2nd alignment and there are 2nd alignments in region which have primary alignment outside that region]
  feat.secAlignment.prop.1 <- length(secondaryAlignmentsFromRegion.N) / (NROW(mapping.reg.regular.N)+1L)
  feat.secAlignment.prop.d <- length(secondaryAlignmentsFromRegion.T) / (NROW(mapping.reg.regular.T)+1L) - feat.secAlignment.prop.1  #which(bitAnd(mapping.2$flag, FLAG_SECONDARY) == FLAG_SECONDARY)
  
  
  # all templates anywhere (identified by QNAME) that have supplementary alignments with a link to the region via their read-pair: 
  # they either map into that region as supplementary or their primary mapping is there but they have a supplementary mapping elsewhere
  supplAlignmentsFromRegion.N <- intersect(patMappingInfo[["wt"]][["supplReads"]], mapping.reg.N$qname)
  supplAlignmentsFromRegion.T <- intersect(patMappingInfo[["mut"]][["supplReads"]], mapping.reg.T$qname)
  
  # prop of supplementary read mappings relative to primary reads on selected region:
  # theoretically this number can go beyond 1 [if all primary read pairs in region have a supplementary alignment and there are supplementary alignments in region which have primary alignment outside that region]
  feat.supplAlignment.prop.1 <- length(supplAlignmentsFromRegion.N) / (NROW(mapping.reg.regular.N)+1L)
  feat.supplAlignment.prop.d <- length(supplAlignmentsFromRegion.T) / (NROW(mapping.reg.regular.T)+1L) - feat.supplAlignment.prop.1
  
  
  
  # proportion of Q1 reads in region relative to all mapping entries in region 
  feat.q1Reads.prop.1 <- length(which(mapping.reg.regular.N$mapq>0L)) / NROW(mapping.reg.N)
  feat.q1Reads.prop.d <- length(which(mapping.reg.regular.T$mapq>0L)) / NROW(mapping.reg.T) - feat.q1Reads.prop.1
  
  
  # mean and sd of mapping quality in selected region
  feat.mapq.avg.1 <- mean(mapping.reg.regular.N$mapq, na.rm=TRUE)
  feat.mapq.avg.d <- mean(mapping.reg.regular.T$mapq, na.rm=TRUE) - feat.mapq.avg.1
  
  feat.mapq.sd.1 <- sd(mapping.reg.regular.N$mapq, na.rm=TRUE)
  feat.mapq.sd.d <- sd(mapping.reg.regular.T$mapq, na.rm=TRUE) - feat.mapq.sd.1
  
  feat.pe.isize.avg.1 <- median(abs(mapping.reg.regular.N$isize), na.rm=TRUE)
  feat.pe.isize.avg.d <- median(abs(mapping.reg.regular.T$isize), na.rm=TRUE) - feat.pe.isize.avg.1
  if ( is.na(feat.pe.isize.avg.d) ) feat.pe.isize.avg.d <- 0L
  
  feat.pe.isize.sd.1 <- mad(mapping.reg.regular.N$isize[mapping.reg.regular.N$isize>0L & mapping.reg.regular.N$isize<MAX_PAIR_WIDTH], na.rm=TRUE)
  feat.pe.isize.sd.d <- mad(mapping.reg.regular.T$isize[mapping.reg.regular.T$isize>0L & mapping.reg.regular.T$isize<MAX_PAIR_WIDTH], na.rm=TRUE) - feat.pe.isize.sd.1
  if ( is.na(feat.pe.isize.sd.d) ) feat.pe.isize.sd.d <- 0L
  
  # proportion of reads that are mapped itself within the region range but have unmapped mate
  #feat.pe.singleAnchor.prop.1 <- sum(bitops::bitAnd(mapping.reg.regular.N$flag, FLAG_MATE_UNMAPPED) == FLAG_MATE_UNMAPPED & bitops::bitAnd(mapping.reg.regular.N$flag, FLAG_READ_UNMAPPED) == 0L) / NROW(mapping.reg.regular.N)
  #feat.pe.singleAnchor.prop.d <- sum(bitops::bitAnd(mapping.reg.regular.T$flag, FLAG_MATE_UNMAPPED) == FLAG_MATE_UNMAPPED & bitops::bitAnd(mapping.reg.regular.T$flag, FLAG_READ_UNMAPPED) == 0L) / NROW(mapping.reg.regular.T) - feat.pe.singleAnchor.prop.1
  
  # use extended selected region to also get read pairs where the mapped read end is next to the selected region
  feat.pe.singleAnchor.prop.1 <- sum(bitops::bitAnd(mapping.regionExt.N$flag, FLAG_MATE_UNMAPPED) == FLAG_MATE_UNMAPPED & bitops::bitAnd(mapping.regionExt.N$flag, FLAG_READ_UNMAPPED) == 0L) / NROW(mapping.regionExt.N)
  feat.pe.singleAnchor.prop.d <- sum(bitops::bitAnd(mapping.regionExt.T$flag, FLAG_MATE_UNMAPPED) == FLAG_MATE_UNMAPPED & bitops::bitAnd(mapping.regionExt.T$flag, FLAG_READ_UNMAPPED) == 0L) / NROW(mapping.regionExt.T) - feat.pe.singleAnchor.prop.1
  
  
  featureSet_reads <- c(feat.region.length.r=feat.region.length.r,
                        feat.region.targetDist.5p.r=feat.region.targetDist.5p.r, feat.region.targetDist.3p.r=feat.region.targetDist.3p.r,
                        feat.readProp.1=feat.readProp.1, feat.readProp.d=feat.readProp.d,
                        feat.readProp.discordant.1=feat.readProp.discordant.1, feat.readProp.discordant.d=feat.readProp.discordant.d,
                        feat.readProp.2nd.1=feat.readProp.2nd.1, feat.readProp.2nd.d=feat.readProp.2nd.d, #mkuhn, 2015-12-15 renamed from sec. to 2nd. for clarity
                        feat.readProp.suppl.1=feat.readProp.suppl.1, feat.readProp.suppl.d=feat.readProp.suppl.d,
                        feat.secAlignment.prop.1=feat.secAlignment.prop.1, feat.secAlignment.prop.d=feat.secAlignment.prop.d,
                        feat.supplAlignment.prop.1=feat.supplAlignment.prop.1, feat.supplAlignment.prop.d=feat.supplAlignment.prop.d,
                        feat.q1Reads.prop.1=feat.q1Reads.prop.1, feat.q1Reads.prop.d=feat.q1Reads.prop.d,
                        feat.mapq.avg.1=feat.mapq.avg.1, feat.mapq.avg.d=feat.mapq.avg.d,
                        feat.mapq.sd.1=feat.mapq.sd.1, feat.mapq.sd.d=feat.mapq.sd.d, 
                        feat.pe.isize.avg.1=feat.pe.isize.avg.1, feat.pe.isize.avg.d=feat.pe.isize.avg.d,
                        feat.pe.isize.sd.1=feat.pe.isize.sd.1, feat.pe.isize.sd.d=feat.pe.isize.sd.d,
                        feat.pe.singleAnchor.prop.1=feat.pe.singleAnchor.prop.1, feat.pe.singleAnchor.prop.d=feat.pe.singleAnchor.prop.d)
  
  
  
  
  
  
  
  # CLIPPED BASE features -----
  
  # mkuhn, 2015-03-03
  #  clipped read peaks in genomic coordinates from screening
  #  mainly based on tumor probe (hence .t ending)
  
  ## mkuhn, 2015-11-27: old version, when peakInfo was a dataframe of 0, 1 or 2 rows
  ## feat.clipped.peak.n.t <- NROW(peakInfo)
  ## # 0 has a double meaning as it is a valid posDiff (meaning resolves based on feat.clipped.peak.n.t)
  ## feat.clipped.peak.posdiff.t <- if (NROW(peakInfo) >= 1L) peakInfo$posDiff[1L] else 0L
  #peakInfo <- clipInfo[["peakInfo"]] # clipInfo is a vector now, not a list any more
  feat.clipped.peakNbr.t <- clipInfo[["nbrPeaks"]] #[[ ]] will strip name
  # distance between the 3' and 5' peaks
  # posDiff shifted by +/-1. This is to discriminate a real peak with pos-diff of 0 (becomes 1) from no peak Info (becomes 0)
  feat.clipped.peak.posdiff.t <- if ( clipInfo[["nbrPeaks"]] >= 1L) clipInfo[["posDiff"]] + sign(clipInfo[["posDiff"]]) + (clipInfo[["posDiff"]] == 0L) else 0L 
  # mkuhn, 2016-02-09: log-transformation because of outliers
  feat.clipped.peak.posdiff.t <- if (feat.clipped.peak.posdiff.t != 0L) log10(abs(feat.clipped.peak.posdiff.t) + 3L) * sign(feat.clipped.peak.posdiff.t) else 0L
  
  # proportion of clipped bases at peak
  feat.clipped.prop.t <- mean(clipInfo[["relS_5p.99"]], clipInfo[["relS_3p.99"]])
  # mkuhn, 2016-01-10: idea to add minimum clipped peak score as feature
  #feat.clipped.minScore.t <- clipInfo[["minScore"]]
  
  
  
  ## CIGAR string of all regular reads
  # (for clipped/inserted/deleted bases)
  cigar.N <- as.character(mapping.reg.regular.N$cigar)
  cigar.T <- as.character(mapping.reg.regular.T$cigar)
  
  # no hard-clipping  (see BWA MEM -Y option)
  stopifnot( ! any(grepl("H", cigar.N)) )
  stopifnot( ! any(grepl("H", cigar.T)) )
  
  
  # cigar string for Q1-mapped regular reads:
  cigar.q1.N <- cigar.N[which(mapping.reg.regular.N$mapq >= 1L)]
  cigar.q1.T <- cigar.T[which(mapping.reg.regular.T$mapq >= 1L)]
  
  # proportion of clipped reads within the Q1-reads
  feat.clipped.q1.prop.1 <- if (length(cigar.q1.N) > 0L) sum(grepl("S", x=cigar.q1.N, fixed=TRUE)) / length(cigar.q1.N) else NA
  feat.clipped.q1.prop.d <- if (length(cigar.q1.T) > 0L) sum(grepl("S", x=cigar.q1.T, fixed=TRUE)) / length(cigar.q1.T) - feat.clipped.q1.prop.1 else NA
  
  # sum of all clipped S-bases within the Q1-reads (anywhere in the read)
  nbrS.q1.cigar.N <- sapply( stringr::str_extract_all(cigar.q1.N, pattern="[[:digit:]]+S"), function(x) sum(as.numeric(stringr::str_sub(x, end=-2L))) )
  nbrS.q1.cigar.T <- sapply( stringr::str_extract_all(cigar.q1.T, pattern="[[:digit:]]+S"), function(x) sum(as.numeric(stringr::str_sub(x, end=-2L))) )
  
  # number of all S-bases relative to number of Q1-read entries
  feat.clipped.q1.bases.1 <- if (length(cigar.q1.N) > 0L)  sum(nbrS.q1.cigar.N, na.rm=TRUE) / length(cigar.q1.N) else NA
  feat.clipped.q1.bases.d <- if (length(cigar.q1.T) > 0L) (sum(nbrS.q1.cigar.T, na.rm=TRUE) / length(cigar.q1.T)) - feat.clipped.q1.bases.1 else NA
  
  #   # sum of clipped bases within the Q1-reads, from either end
  #   nbrS.q1.5p.cigar.N <- as.numeric(str_sub(str_extract(cigar.q1.N, pattern = "^[[:digit:]]+S"), end=-2L))
  #   nbrS.q1.5p.cigar.T <- as.numeric(str_sub(str_extract(cigar.q1.T, pattern = "^[[:digit:]]+S"), end=-2L))
  #   
  #   nbrS.q1.3p.cigar.N <- as.numeric(str_sub(str_extract(cigar.q1.N, pattern = "[[:digit:]]+S$"), end=-2L))
  #   nbrS.q1.3p.cigar.T <- as.numeric(str_sub(str_extract(cigar.q1.T, pattern = "[[:digit:]]+S$"), end=-2L))
  #   
  #   
  #   #ZZZ Or should I simplify and sum both clipped bases counts (from 3' and 5') in a single feature?
  #   feat.clipped.q1.bases.5p.1 <- sum(nbrS.q1.5p.cigar.N, na.rm=TRUE)
  #   feat.clipped.q1.bases.5p.d <- sum(nbrS.q1.5p.cigar.T, na.rm=TRUE) - feat.clipped.q1.bases.5p.1
  #   
  #   feat.clipped.q1.bases.3p.1 <- sum(nbrS.q1.3p.cigar.N, na.rm=TRUE)
  #   feat.clipped.q1.bases.3p.d <- sum(nbrS.q1.3p.cigar.T, na.rm=TRUE) - feat.clipped.q1.bases.3p.1
  
  
  #ZZZ spread of clipped bases across pos!
  #ZZZ use PV= mean{ ABS(xi - xj) / max(xi,xj)} over all pairs xi,xj as robust measure of variability in a series of data
  
  
  
  # INSERTED BASES in reads ------
  
  # proportion of reads with insertion within the Q1-reads
  feat.ins.q1.prop.1 <- if (length(cigar.q1.N) > 0L) sum(grepl("I", x=cigar.q1.N, fixed=TRUE)) / length(cigar.q1.N) else NA
  feat.ins.q1.prop.d <- if (length(cigar.q1.T) > 0L) sum(grepl("I", x=cigar.q1.T, fixed=TRUE)) / length(cigar.q1.T) - feat.ins.q1.prop.1 else NA
  
  # sum of all inserted bases within the Q1-reads (anywhere in the read)
  nbrI.q1.cigar.N <- sapply(stringr::str_extract_all(cigar.q1.N, pattern="[[:digit:]]+I"),
                            function(x) sum(as.numeric(stringr::str_sub(x, end=-2L))) )
  nbrI.q1.cigar.T <- sapply(stringr::str_extract_all(cigar.q1.T, pattern="[[:digit:]]+I"),
                            function(x) sum(as.numeric(stringr::str_sub(x, end=-2L))) )
  
  # relative to number of Q1-read entries
  feat.ins.q1.bases.1 <- if ( length(nbrI.q1.cigar.N) > 0L) sum(nbrI.q1.cigar.N, na.rm=TRUE) / length(cigar.q1.N) else NA
  feat.ins.q1.bases.d <- if ( length(nbrI.q1.cigar.T) > 0L) (sum(nbrI.q1.cigar.T, na.rm=TRUE) / length(cigar.q1.T)) - feat.ins.q1.bases.1 else NA
  
  
  
  # DELETED BASES in reads -----
  
  # proportion of reads with deletion within the Q1-reads
  feat.del.q1.prop.1 <- if (length(cigar.q1.N) > 0L) sum(grepl("D", x=cigar.q1.N, fixed=TRUE)) / length(cigar.q1.N) else NA
  feat.del.q1.prop.d <- if (length(cigar.q1.T) > 0L) sum(grepl("D", x=cigar.q1.T, fixed=TRUE)) / length(cigar.q1.T) - feat.del.q1.prop.1 else NA
  
  # sum all deleted bases within the Q1-reads (anywhere in the read)
  nbrD.q1.cigar.N <- sapply(stringr::str_extract_all(cigar.q1.N, pattern="[[:digit:]]+D"),
                            function(x) sum(as.numeric(stringr::str_sub(x, end=-2L))) )
  nbrD.q1.cigar.T <- sapply(stringr::str_extract_all(cigar.q1.T, pattern="[[:digit:]]+D"),
                            function(x) sum(as.numeric(stringr::str_sub(x, end=-2L))) )
  
  # relative to number of Q1-read entries
  feat.del.q1.bases.1 <- if ( length(nbrD.q1.cigar.N) > 0L) sum(nbrD.q1.cigar.N, na.rm=TRUE) / length(cigar.q1.N) else NA
  feat.del.q1.bases.d <- if ( length(nbrD.q1.cigar.T) > 0L) (sum(nbrD.q1.cigar.T, na.rm=TRUE) / length(cigar.q1.T)) - feat.del.q1.bases.1 else NA
  
  
  featureSet_cigar <- c(feat.clipped.peakNbr.t=feat.clipped.peakNbr.t, #feat.clipped.peak.score.t=feat.clipped.peak.score.t,
                        feat.clipped.peak.posdiff.t=feat.clipped.peak.posdiff.t, 
                        feat.clipped.prop.t=feat.clipped.prop.t,
                        feat.clipped.q1.prop.1=feat.clipped.q1.prop.1, feat.clipped.q1.prop.d=feat.clipped.q1.prop.d,
                        feat.clipped.q1.bases.1=feat.clipped.q1.bases.1, feat.clipped.q1.bases.d=feat.clipped.q1.bases.d,
                        #feat.clipped.q1.bases.3p.1=feat.clipped.q1.bases.3p.1, feat.clipped.q1.bases.3p.d=feat.clipped.q1.bases.3p.d,
                        #feat.clipped.q1.bases.5p.1=feat.clipped.q1.bases.5p.1, feat.clipped.q1.bases.5p.d=feat.clipped.q1.bases.5p.d,
                        feat.ins.q1.prop.1=feat.ins.q1.prop.1, feat.ins.q1.prop.d=feat.ins.q1.prop.d, 
                        feat.ins.q1.bases.1=feat.ins.q1.bases.1, feat.ins.q1.bases.d=feat.ins.q1.bases.d,
                        feat.del.q1.prop.1=feat.del.q1.prop.1, feat.del.q1.prop.d=feat.del.q1.prop.d, 
                        feat.del.q1.bases.1=feat.del.q1.bases.1, feat.del.q1.bases.d=feat.del.q1.bases.d)
  
  
  
  
  
  
  # COVERAGE (in SELECTED_REGION) -----
  
  QUAL_FILTER <- "-Q1" # "-Q1" # There were cases where Q1 sometimes filters out everything!
  # returns an unparsed string of CHROM\tPOS\tCOUNT
  cov.samtools.N <- system(paste(SAMTOOLS_EXE, "depth ", QUAL_FILTER, "-r ", SELECTED_REGION_GEN, bamFile_WT), intern=TRUE)
  cov.samtools.T <- system(paste(SAMTOOLS_EXE, "depth ", QUAL_FILTER, "-r ", SELECTED_REGION_GEN, bamFile_MUT), intern=TRUE)
  
  # no-coverage matrix for region
  zeroCov <- matrix(c(seq(startPos, endPos), rep(0L, feat.region.length.r)), ncol=2L)
  
  
  # get smoothed coverage
  if ( length(cov.samtools.N) == 0L ){
    logwarn("No coverage (%s) in genomic region %s for mapping %s.", QUAL_FILTER, SELECTED_REGION_GEN, bamFile_WT)
    #stop("No coverage in [N]!")
    #return(invisible(NULL)) 
    cov.N <- zeroCov
  } else {
    # assuming a single fixed chromosome for the region and also not empty
    chrom.N <- unique(stringr::str_replace(cov.samtools.N, "\t.*", ""))
    stopifnot(  length(chrom.N) == 1L ) #chrom.N == chrom.T 
    cov.raw.N <- matrix(as.numeric(unlist(strsplit(stringr::str_replace(cov.samtools.N, paste0("^", chrom.N, "\t"), replacement = ""),
                                                   split="\t", fixed=TRUE))),
                        ncol=2L, byrow=TRUE)
    
    # keep position and coverage, but fill zero cov positions
    cov.N <- rbind(cov.raw.N, zeroCov[! zeroCov[,1L] %in% cov.raw.N[,1L],])
    cov.N <- cov.N[order(cov.N[,1L]),]
  }
  
  
  # cov.T0 is for raw (non-normalized) coverage
  if ( length(cov.samtools.T) == 0L ){
    logwarn("No coverage (%s) in genomic region %s for mapping %s.", QUAL_FILTER, SELECTED_REGION_GEN, bamFile_MUT)
    #stop("No coverage in [T]!")
    #return(invisible(NULL))
    cov.T0 <- cov.T <- zeroCov
    
  } else {
    # assuming a single fixed chromosome for the region and also not empty
    chrom.T <- unique(stringr::str_replace(cov.samtools.T, "\t.*", ""))
    stopifnot(  length(chrom.T) == 1L ) #chrom.N == chrom.T 
    
    cov.raw.T <- matrix(as.numeric(unlist(strsplit(stringr::str_replace(cov.samtools.T, paste0("^", chrom.T, "\t"), replacement = ""),
                                                   split="\t", fixed=TRUE))),
                        ncol=2L, byrow=TRUE)
    
    cov.T <- rbind(cov.raw.T, zeroCov[! zeroCov[,1L] %in% cov.raw.T[,1L],])
    cov.T0 <- cov.T <- cov.T[order(cov.T[,1L]),]
    # mkuhn, 2015-03-03
    # normalize coverage based on overall read numbers of whole BAMs
    cov.T[,2L] <- cov.T[,2L] / rc.bam.ratio
  }
  
  
  # every position of the region should be listed precisely once
  stopifnot( NROW(cov.T) == feat.region.length.r, NROW(cov.N) == feat.region.length.r )
  
  
  # median coverage satisfying QUAL_FILTER in candidate region relative to average coverage in target region
  # using cov.raw.X because cov.T was normalized for count-differences betw T and N. Here we have "internal normalization" by read count on target
  # mkuhn, 2015-12-26: including zero-coverages for both N and T. Use raw data (not normalized for read count) for [T]
  feat.cov.ratio.1 <- median(cov.N[,2L]) / ( rc.target.wt * READL / (targetEndPos - targetStartPos) )
  feat.cov.ratio.d <- ( median(cov.T0[,2L]) / ( rc.target.mut * READL / (targetEndPos - targetStartPos) )) - feat.cov.ratio.1
  
  
  # proportion of bases in selected region that surpass certain coverage
  # fixed thresholds are OK because we have normalized cov.T based on the read count in the normal WT-BAM
  # delivers many correlated features?!
  feat.cov.base01.prop.1 <- sum(cov.N[,2L] >= 1L) / NROW(cov.N)
  feat.cov.base01.prop.d <- sum(cov.T[,2L] >= 1L) / NROW(cov.T) - feat.cov.base01.prop.1
  feat.cov.base05.prop.1 <- sum(cov.N[,2L] >= 5L) / NROW(cov.N)
  feat.cov.base05.prop.d <- sum(cov.T[,2L] >= 5L) / NROW(cov.T) - feat.cov.base05.prop.1
  #     
  #     feat.cov.base15.prop.1 <- sum(cov.N[,2] >= 15) / NROW(cov.N)
  #     feat.cov.base15.prop.d <- sum(cov.T[,2] >= 15) / NROW(cov.T) - feat.cov.base15.prop.1
  #     
  #     feat.cov.base25.prop.1 <- sum(cov.N[,2] >= 25) / NROW(cov.N)
  #     feat.cov.base25.prop.d <- sum(cov.T[,2] >= 25) / NROW(cov.T) - feat.cov.base25.prop.1
  #     
  #     feat.cov.base50.prop.1 <- sum(cov.N[,2] >= 50) / NROW(cov.N)
  #     feat.cov.base50.prop.d <- sum(cov.T[,2] >= 50) / NROW(cov.T) - feat.cov.base50.prop.1
  
  # how steep is the ecdf betw 25% and 75% of coverage values? #should I use smoothed cov-values?
  feat.cov.iqr.1 <- IQR(cov.N[,2L])
  feat.cov.iqr.d <- IQR(cov.T[,2L]) - feat.cov.iqr.1
  
  # what percentage of bases in candidate region is covered more than average coverage in candidate region
  # this would be sensitive to coverage outliers
  feat.cov.baseAvg.prop.1 <- sum(cov.N[,2L] > mean(cov.N[,2L])) / NROW(cov.N)
  feat.cov.baseAvg.prop.d <- sum(cov.T[,2L] > mean(cov.T[,2L])) / NROW(cov.T) - feat.cov.baseAvg.prop.1
  
  
  
  
  
  # Coverage smoothing
  
  # Idea: use wide and small window. Now, it is adhoc.
  # Window size: aim at 20% of region length, at max 25, at least 5.
  #
  rollmean_window <- rollsd_window <- max(5L, min(25L, floor(feat.region.length.r / 5L))) 
  hammWeight.mean <- getHammingWindowWeights(rollmean_window)
  hammWeight.sd <- getHammingWindowWeights(rollsd_window)
  
  cov.roll.N <- RcppRoll::roll_mean(cov.N, n=rollmean_window, weights = hammWeight.mean)
  cov.roll.N <- cbind(cov.roll.N,
                      RcppRoll::roll_sd(cov.N[,2L], n=rollsd_window, weights = hammWeight.sd))
  colnames(cov.roll.N) <- c("pos", "covMean.N", "covSD.N")
  cov.roll.T <- RcppRoll::roll_mean(cov.T, n=rollmean_window, weights = hammWeight.mean)
  cov.roll.T <- cbind(cov.roll.T,
                      RcppRoll::roll_sd(cov.T[,2L], n=rollsd_window, weights = hammWeight.sd))
  colnames(cov.roll.T) <- c("pos", "covMean.T", "covSD.T")
  
  # we do not need position twice
  stopifnot( isTRUE(all.equal(cov.roll.N[,1L], cov.roll.T[,1L])) )
  cov.roll <- cbind(cov.roll.N, cov.roll.T[,-1L])
  # inner_join(cov.roll.N, cov.roll.T) # if we had data-frames
  
  # get differences [T] - [N] in mean coverage and SD
  cov.roll <- cbind(cov.roll,
                    covMean.d=cov.roll[,"covMean.T"] - cov.roll[,"covMean.N"],
                    covSD.d=cov.roll[,"covSD.T"] - cov.roll[,"covSD.N"],
                    covCV.N=cov.roll[, "covSD.N"] / (cov.roll[, "covMean.N"]+1L),
                    covCV.T=cov.roll[, "covSD.T"] / (cov.roll[, "covMean.T"]+1L) )
  
  #cov.roll2 <- cbind(cov.roll,
  # overview plot:
  #     YCOLS <- c("covMean.N", "covSD.N", "covMean.T", "covSD.T")
  #     matplot(x = cov.roll[, "pos"], y=cov.roll[,YCOLS], #, "covMean.d", "covSD.d"
  #      type="l", col = c(1, 1, 2, 2), lty=c(1,2,1,2), xlab="Pos", ylab="Coverage");  #3, 3)
  #     legend("bottom", YCOLS, col=c(1, 1, 2, 2), lty=c(1,2,1,2), inset = c(0, 0.01))
  
  
  ## Coverage features:
  
  # overall difference in mean coverage at selected region
  # cov means were normalized, so we can use the mean difference of smoothed coverage as feature
  feat.cov.d <- mean(cov.roll[, "covMean.d"])
  
  # overall spread across all Q1-coverage in selected region on un-smoothed base-resolution coverage
  # or should I use cov.roll?!
  # mad() as robust version of sd()
  feat.cov.sd.1 <- mad(cov.N[,2L], na.rm = TRUE)  
  feat.cov.sd.d <- mad(cov.T[,2L], na.rm = TRUE) - feat.cov.sd.1
  
  # features from rolling windows on coverage: maximal difference in coverage betw [T] vs [N]
  # mkuhn, 2015-04-20: replace feat.cov[.sd].max.d feature for cov (resp. sd) with each two features, namely mean of cov (resp. sd) of 10% and 90% quantile region
  #
  # Idea: use information of which.max() that might give a hint where the BP is!
  #feat.cov.max.d <- as.numeric(cov.roll[which.max(abs(cov.roll[,"covMean.d"])), "covMean.d"])
  #feat.cov.sd.max.d  <- as.numeric(cov.roll[which.max(abs(cov.roll[,"covSD.d"])), "covSD.d"])
  
  covMean.diff.quants <- quantile(cov.roll[,"covMean.d"], probs = c(0.1, 0.9))
  covSD.diff.quants <- quantile(cov.roll[, "covSD.d"], probs = c(0.1, 0.9))
  # mean smoothed coverage difference on bases of 10% smallest and 10% highest coverage difference
  feat.cov.q10.d <- median(cov.roll[which(cov.roll[, "covMean.d"] <= covMean.diff.quants[1]), "covMean.d"], na.rm = TRUE)
  feat.cov.q90.d <- median(cov.roll[which(cov.roll[, "covMean.d"] >= covMean.diff.quants[2]), "covMean.d"], na.rm = TRUE)
  
  feat.cov.sd.q10.d <- median(cov.roll[which(cov.roll[, "covSD.d"] <= covSD.diff.quants[1]), "covSD.d"], na.rm = TRUE)
  feat.cov.sd.q90.d <- median(cov.roll[which(cov.roll[, "covSD.d"] >= covSD.diff.quants[2]), "covSD.d"], na.rm = TRUE)
  
  
  # mkuhn, 2015-03-03
  # CV of smoothed coverage should be constant (as spread grows with count)
  # spread of CV
  feat.cov.cv.sd.1 <- mad(cov.roll[, "covCV.N"], na.rm = TRUE)
  feat.cov.cv.sd.d <- mad(cov.roll[, "covCV.T"], na.rm = TRUE) - feat.cov.cv.sd.1
  
  
  # mkuhn, 2015-03-03
  # proportion of positions where normalized mean coverage in [T] deviates pseudo-"significantly" from coverage in [N]
  myQuantil <- 2.58 #qnorm(.995)
  ind.cov.sigdiff.lower <- which(cov.roll[, "covMean.T"] < cov.roll[, "covMean.N"] - myQuantil * (cov.roll[, "covSD.N"]+.1) / sqrt(rollmean_window) )
  ind.cov.sigdiff.upper <- which(cov.roll[, "covMean.T"] > cov.roll[, "covMean.N"] + myQuantil * (cov.roll[, "covSD.N"]+.1) / sqrt(rollmean_window) )
  feat.cov.sig.lower.d <- length(ind.cov.sigdiff.lower) / NROW(cov.roll)
  feat.cov.sig.upper.d <- length(ind.cov.sigdiff.upper) / NROW(cov.roll)
  
  
  
  
  
  
  featureSet_cov <- c(feat.cov.d=feat.cov.d,
                      feat.cov.sd.1=feat.cov.sd.1, feat.cov.sd.d=feat.cov.sd.d, 
                      #feat.cov.max.d=feat.cov.max.d, feat.cov.sd.max.d=feat.cov.sd.max.d, 
                      feat.cov.10.d=feat.cov.q10.d, feat.cov.90.d=feat.cov.q90.d, #ZZZ .q90.d iso .90.d?
                      feat.cov.sd.10.d=feat.cov.sd.q10.d, feat.cov.sd.90.d=feat.cov.sd.q90.d,
                      feat.cov.cv.sd.1=feat.cov.cv.sd.1, feat.cov.cv.sd.d=feat.cov.cv.sd.d,
                      feat.cov.sig.lower.d=feat.cov.sig.lower.d, feat.cov.sig.upper.d=feat.cov.sig.upper.d,
                      feat.cov.ratio.1=feat.cov.ratio.1, feat.cov.ratio.d=feat.cov.ratio.d,
                      feat.cov.base01.prop.1=feat.cov.base01.prop.1, feat.cov.base01.prop.d=feat.cov.base01.prop.d,
                      feat.cov.base05.prop.1=feat.cov.base05.prop.1, feat.cov.base05.prop.d=feat.cov.base05.prop.d,
                      #feat.cov.base15.prop.1=feat.cov.base15.prop.1, feat.cov.base15.prop.d=feat.cov.base15.prop.d,
                      #feat.cov.base25.prop.1=feat.cov.base25.prop.1, feat.cov.base25.prop.d=feat.cov.base25.prop.d, 
                      #feat.cov.base50.prop.1=feat.cov.base50.prop.1, feat.cov.base50.prop.d=feat.cov.base50.prop.d,
                      feat.cov.iqr.1=feat.cov.iqr.1, feat.cov.iqr.d=feat.cov.iqr.d,
                      feat.cov.baseAvg.prop.1=feat.cov.baseAvg.prop.1, feat.cov.baseAvg.prop.d=feat.cov.baseAvg.prop.d)
  
  
  
  # finally,
  # return the extracted features as named numeric vector
  return( c( idInfo, featureSet_reads, featureSet_cov, featureSet_cigar) )
}




#' Extracts the features for different candidate loci from a tumor/normal mapping of a patient.
#' 
#' This is a wrapper around \code{\link{getFeaturesFromMapping}} for a whole set of regions.
#' In a previous version I used \code{statusInfo} a vector containing status-info as character vector as argument. But took it out because it only is available for simulation patients.
#' Parallel processing of the different loci was tried but there were concurrency problems. Hence, it works now sequentially, one loci at a time.
#' 
#' @param patId patient ID as numeric (like 13)
#' @param bamFile_WT filename of mapping from normal sample
#' @param bamFile_MUT filename of mapping from tumor sample
#' @param patMappingInfo list with certain information based on whole mapping in tumor and normal sample
#' @param chrom Chromosome of the candidate loci. Currently, only single chromosome implemented. chr5 is used.
#' @param regions IRanges object with candidate loci. A candidate locus should be contained within a target region. 
#' @param targetRegions IRanges object with target regions that are enriched by the exome sequencing kit. Candidate locus lies within a target region.
#' @param clipInfos list with clipping info of the tumor probe, per region a list entry 
#' @return a dataframe of features for the two given mapping files. Patient and genetic loci information comes first.
#' @export
getFeaturesFromMappings <- function(patId, bamFile_WT, bamFile_MUT, patMappingInfo,
                                    chrom=CHROM, regions, targetRegions=regions, clipInfos){
  
  stopifnot( inherits(regions, "IRanges"), inherits(targetRegions, "IRanges") )
  
  ## for now: we have support for only a single chromosome. we expect chr5
  stopifnot( length(chrom) == 1L )# == 'chr5' )  
  
  if (length(regions) == 0L){
    logwarn("Got no regions where to extract features.")
    return(invisible(NULL))
  } else {
    logdebug("Got %d regions to extract features.", length(regions))
  }
  
  # Start and end coordinates of regions
  regStart <- start(regions)
  regEnd <- end(regions)
  
  targetStart <- start(targetRegions)
  targetEnd <- end(targetRegions)
  
  # mkuhn, 20150423: check length constraints
  stopifnot( length(regStart) == length(targetStart), length(regStart) == length(clipInfos), length(regStart) == length(regEnd) )
  stopifnot( all(regStart < regEnd) )
  stopifnot( all(regStart >= targetStart & regEnd <= targetEnd) )
  
  #bamFile_WT <- "../../SVdata//mapping//cov60_TL90//PEm200sd25_R100//pat00001_N.bam"
  #bamFile_MUT <- "../../SVdata//mapping//cov60_TL90//PEm200sd25_R100//pat00001_T.bam"
  # # DEBUGING
  # for (i in seq_along(regStart)){
  #   print(i)
  #   startR <- regStart[i]; endR <- regEnd[i]; startT <- targetStart[i]; endT <- targetEnd[i]; clipi <- clipInfos[[i]]
  featureMat  <- foreach(startR=regStart, endR=regEnd, startT=targetStart, endT=targetEnd, clipi=clipInfos, .combine = 'rbind') %do% {
    
    getFeaturesFromMapping(patId=patId, bamFile_WT=bamFile_WT, bamFile_MUT=bamFile_MUT, patMappingInfo = patMappingInfo, 
                           chrom=chrom, startPos=startR, endPos=endR, targetStartPos = startT, targetEndPos = endT, clipInfo=clipi)
  }
  
  #   # name features according to the candidate region it comes from
  #   rownames(featureMat) <- paste0(chrom,":",regStart, "-", regEnd)
  
  # mkuhn, 2016-01-04: bug fix
  #+if only a single interval given, then result is not a matrix but a vector!!
  featureDat <- data.frame(if (is.vector(featureMat)) t(featureMat) else featureMat)
  featureDat$chrom <- chrom
  #featureDat$status <- factor(statusInfo, levels = c("germline", "DSV"))
  
  # mkuhn, 2015-04-22: order columns so that patient and genomic locus information comes first
  featureDat1 <- featureDat %>% dplyr::select_( ~c(patId, chrom, startPos, endPos, targetStartPos, targetEndPos)) #, status))
  #, margin,
  #cov, tumorLoad, read.length, pe.ins.mean, pe.ins.sd) )
  featureDat.new <- cbind(featureDat1, featureDat[, ! names(featureDat) %in% names(featureDat1)])
  rownames(featureDat.new) <- NULL
  
  featureDat.new
}


#' Combines stored feature data from different simulation scenarios.
#' 
#' It merges also patient meta-data if feature data does \emph{not} stem from screening run where SVs have been introduced by SV-injection.
#' 
#' @param baseDir base directory
#' @param featureFileNamePattern pattern string which feature files to match
#' @param fromScreen logical flag: do you want to have feature data from screening run?
#' @param marginValues numeric given the margin value used for feature extraction
#' @param covValues numeric. if specified only features for the given coverages are merged.
#' @param writeOutFeatures flag if combined features should be written out to disk
#' @return dataframe that contains features from different scenarios and patient meta-data (status) or \code{NULL}
harvestFeatures <- function(baseDir=svmod::BASEDIR, featureFileNamePattern=c("features.+[.]rds$","features.*[.]rds$"), fromScreen=TRUE, marginValues=NULL, covValues=NULL, writeOutFeatures=FALSE){
  featureFileNamePattern <- match.arg(featureFileNamePattern)
  
  patDataFile <- file.path(getVirtualPatientPath(baseDir, "patient"), "patData.rds")
  stopifnot( file.exists(patDataFile))
  
  patData <- readRDS(patDataFile)
  stopifnot( NROW(patData) > 3L)
  
  featureBasePath <- file.path(baseDir, "feature")
  featureFiles <- list.files(path=featureBasePath, pattern = featureFileNamePattern, recursive = TRUE, full.names = TRUE)
  
  if (length(featureFiles) == 0L){
    logging::logwarn("No feature files found.. Sorry.")
    return(invisible(NULL))
  }
  
  # restrict to matching filenames
  featureFiles <- grep(pattern = "_screen", x = featureFiles, value=TRUE, invert = ! isTRUE(fromScreen))
  
  if (! is.null(marginValues) && is.numeric(marginValues)){
    featureFiles <- grep(pattern = paste0("_margin", marginValues, collapse="|"), x=featureFiles, value=TRUE)
  }
  
  if (! is.null(covValues) && is.numeric(covValues)){
    featureFiles <- grep(pattern = paste0("cov", covValues, collapse="|"), x = featureFiles, value=TRUE)
  }
  
  # mkuhn, can go out!
  #   # check
  #   df.check <- data.frame(ff=featureFiles)
  #   df.check$colNames <- sapply(featureFiles, function(ff) {f <- readRDS(ff); paste(names(f), collapse="**")})
  
  combFeatures <- foreach::foreach(ff=featureFiles, .combine = 'rbind') %do% {
    f <- readRDS(ff)
    sampleProp <- parseSampleProp(s = dirname(ff))
    
    # mkuhn, 2015-12-28: adding sample characteristics (coverage and tumor-load)
    f %<>% dplyr::mutate_(coverage=~sampleProp[["cov"]], 
                          tumorLoad=~sampleProp[["tumorLoad"]])
    
    if (is.data.frame(f)) f else NULL
  }
  
  stopifnot( is.data.frame(combFeatures) )
  logging::loginfo("Found %d features from within %d feature files.", NROW(combFeatures), length(featureFiles))
  
  
  # adding further context information (e.g. simulated patient SVtype)
  if (! isTRUE(fromScreen)){
    logging::loginfo("Combine with information from simulated patient data.. (as we are in non-screening mode)")
    combFeatures <- patData %>%
      dplyr::select_(~patId, ~chrom, ~targetStartPos, ~targetEndPos, ~status, ~SVtype) %>% 
      dplyr::mutate_(status=~factor(status, levels = c("germline", "DSV"), labels = c("IN", "OUT"))) %>% 
      dplyr::inner_join(y = combFeatures, by = c("patId", "chrom", "targetStartPos", "targetEndPos", "status"))
  }
  
  if (isTRUE(writeOutFeatures)){
    featureOutFile <- paste0("allfeatures", if (isTRUE(fromScreen)) "_screen", ".rds") 
    # margin value is already part of the feature dataframe
    #if (! is.null(marginValues)) paste0("_margin", paste(marginValues, collapse = "*")), ".rds")
    saveRDS(combFeatures, file=file.path(featureBasePath, featureOutFile))
    loginfo("Written out all %d features to %s.", NROW(combFeatures), featureBasePath)
  }
  
  return(combFeatures)
}



#' Visualize feature data that are based on simulation runs.
#' 
#' Do train and test data look differently?
#' We have an option to re-evaluate the SV-outcome for test data based on a given evalMargin.
#' By default, the benchmark-script uses for evaluation margin the value of the feature margin.
#' @param patIds numeric vector of patient IDs to look at. Optional.
#' @param featureMargin read in what features, defaults to 25.
#' @param evalMargin use this margin to judge if a test entry matches the simulated SV-truth. Not used if \code{evalRegionAsMatch=TRUE}
#' @param evalRegionAsMatch flag if test data counts as SV-positive if in same region is a simulated SV. Default is FALSE.
lookAtFeatureData <- function(patIds=NULL, sampleProp = sample.prop(), ngsProp = ngs.prop(), featureMargin=25L, evalMargin=NULL, evalRegionAsMatch=FALSE){
  featPath <- getVirtualPatientPath(what = "feature", .sample.prop = sampleProp, .ngs.prop = ngsProp)
  
  stopifnot( is.numeric(featureMargin), featureMargin >= 0L )
  fDat0 <- readRDS(file = file.path(featPath, paste0("features_sim_screen_margin", featureMargin ,".rds")))
  stopifnot( NROW(fDat0) > 0L )
  
  if ( ! is.null(patIds) && is.numeric(patIds)){
    logging::loginfo("Filter features on patients %s.", paste(patIds, collapse="-"))
    fDat0 <- dplyr::filter_(fDat0, ~patId %in% patIds)
  }
  
  # mkuhn, 2016-02-09: take different evalMargin for test data!
  if (! is.null(evalMargin) && is.numeric(evalMargin) && evalMargin >= 0L || isTRUE(evalRegionAsMatch) ){
    fDat0.train <- fDat0 %>% dplyr::filter_(~type == 'train')
    fDat0.test <- fDat0 %>% dplyr::filter_(~type == 'test')
    stopifnot( NROW(fDat0.train) + NROW(fDat0.test) == NROW(fDat0) )
    
    patData <- readRDS(file.path(getVirtualPatientPath(what = "patient"), "patData.rds"))
    fDat1.test <- addStatusToTestFromSim(featDat.test=fDat0.test, patData = patData, myCHROM = CHROM,
                                         evalMargin = evalMargin, regionIsMatch = evalRegionAsMatch, featNames = names(fDat0.train))
    
    stopifnot( all.equal(dim(fDat1.test), dim(fDat0.test)) )
    fDat0 <- rbind(fDat0.train, fDat1.test)
  }
  
  fDat <- dplyr::select_(fDat0, ~type, ~status, ~SVtype, ~ends_with(".r"), ~ends_with(".t"), ~ends_with(".d"))
  
  logging::loginfo("Found %d feature rows, with %d training. %d test cases with %d 'noSV'.",
                   NROW(fDat), sum(fDat$type == 'train'), sum(fDat$type == 'test'), sum(fDat$type == 'test' & fDat$SVtype == 'no'))
  
  
  if (FALSE){
    
    fDist <- dist(dplyr::select_(fDat, ~-type, ~-status, ~-SVtype))
    mdsFit <- cmdscale(fDist, k = 2L)
    
    # mdsCoord1 <- mdsFit[,1L]
    # mdsCoord2 <- mdsFit[,2L]
    # plot(x = mdsCoord1, y = mdsCoord2, main = "Classical MDS on Feature Data", xlab = "1st coord", ylab = "2nd coord",
    #      pch = 1L + as.numeric(fDat$SVtype), col = ifelse(fDat$type == 'train', yes = "black", no = "red"))
    mdsData <- data.frame(mdsFit[, 1:2])
    mdsData <- cbind(mdsData, fDat[, c("type", "status", "SVtype")])
    ggplot2::ggplot(data = mdsData, mapping = ggplot2::aes_(x=~X1, y=~X2, colour = ~type, shape = ~SVtype)) + 
      ggplot2::geom_point() + ggplot2::labs(title="Classical MDS on Feature Data")
  }
  
  
  
  fDat2 <- tidyr::gather_(fDat, key_col = "variable", value_col = "value", 
                          gather_cols = grep(pattern = "[.][td]$", names(fDat), value = TRUE) )
  
  plotObj <- ggplot2::ggplot(data = fDat2, ggplot2::aes(x=SVtype, y=value, fill=type)) + 
    ggplot2::geom_boxplot() +
    ggplot2::facet_wrap(~variable, scales = "free_y") + 
    ggplot2::theme(axis.text.x=ggplot2::element_text(angle=45L, hjust=1L)) + 
    ggplot2::labs(title = sprintf("Features vs SVtype, evalMargin %s, %d test data with %d 'noSV'",
                                  if (isTRUE(evalRegionAsMatch)) "'region'" else if (! is.null(evalMargin)) evalMargin else "asis", sum(fDat$type == 'test'), sum(fDat$type == 'test' & fDat$SVtype == 'no')))
  
  print(plotObj)
}

