


#' Find covered regions within target regions. 
#' 
#' Look for sufficiently covered regions with homogeneous coverage within the specified target regions. It considers the tumor as well as the normal sample.
#' @param bamFile_WT file name of WT bamfile
#' @param bamFile_MUT file name of MUT bamfile
#' @param targetBedFile coordiantes of target regions in BED format
#' @param extra.margin numeric. margin (in bp) around the target regions 
#' @param k window size for smoothing coverage
#' @param minWidth minimal width of regions to be returned (if further coverage requirements are fulfilled)
#' @param minGapWidth minimal gap between two covered regions. Regions with smaller gap will be merged. This is a way to get over accidential (or SV-related) low coverage spots.
#' @param minCov_WT minimal coverage for normal WT-probe after median smoothing
#' @param minAvgCov_MUT minimal normalized number of reads that need to map into the tumor MUT-probe. The normalization is by width of region in terms of read length.
#' I.e. coverage in tumor probe at any position in the region when reads are spread equally with no overlap.
#' @return the covered regions in target on the specified chromosome as an \code{GRanges}-object 
#' @export
findCoveredRegionsInTarget <- function(bamFile_WT, bamFile_MUT, chrom=CHROM,
                                       targetBedFile=file.path(BASEDIR, "refGen", "TruSeq_exome_targeted_regions.hg19.bed"), extra.margin=EXTRA_MARGIN_TARGET_REGIONS, 
                                       k=21L, minWidth=MIN_COV_REGION_WIDTH, minGapWidth=2*READL, minCov_WT=15L, minAvgCov_MUT=5L){
  
  if (FALSE){ #for testing
    # cov60TL50/pat00001 on Toshiba external hard disk: covered region chr5:310430-312500 (found via samtools depth)
    mySampleProp <- sample.prop(cov = 60L, tumorLoad = 95L)
    bamFile_WT  <- paste0(getVirtualPatientPath(what = "mapping", .sample.prop = mySampleProp), getPatStr(1), "_N.bam")
    bamFile_MUT <- paste0(getVirtualPatientPath(what = "mapping", .sample.prop = mySampleProp), getPatStr(1), "_T.bam")
    chrom <- CHROM
    targetBedFile <- file.path(BASEDIR, "refGen", "TruSeq_exome_targeted_regions.hg19.bed")
    k <- 21L
    minWidth <- 50L
    minGapWidth <- READL
    minCov_WT <- 15L
    minAvgCov_MUT <- 5L
    extra.margin <- EXTRA_MARGIN_TARGET_REGIONS
  }
  
  if ( ! checkFile(bamFile_WT) ){
    logerror("BAM file %s of normal probe not found!", bamFile_WT)
    return(invisible(FALSE))
  }
  
  targetRegions <- getTargetRegions(targetBedFile = targetBedFile, extra.margin=extra.margin, chrom = chrom, reduce=TRUE)
  
  
  # Reads from BAM file as GAlignments object: no duplicate reads
  # Cave: Reads will be reported more than once if they cover more than a single target region.
  myParam <- Rsamtools::ScanBamParam(flag = Rsamtools::scanBamFlag(isDuplicate = FALSE),
                                     which = targetRegions)
  
  rga.wt  <- GenomicAlignments::readGAlignments(file=bamFile_WT, param = myParam)
  rga.mut <- GenomicAlignments::readGAlignments(file=bamFile_MUT, param = myParam)
  # median-based smoothing on coverage, package S4Vectors/GenomicRanges is necessary here for runmean(); exported from S4Vectors, runmed in stats
  bamCov.wt  <- runmed(x=Biostrings::coverage(rga.wt), k=k)
  #bamCov.mut <- runmed(x=Biostrings::coverage(rga.mut), k=k) #not needed!
  
  # find areas of interest that show the desired coverage in WT-probe
  hotslices.wt <- GenomicRanges::GRanges(XVector::slice(bamCov.wt, lower=minCov_WT, rangesOnly=TRUE))
  
  # eventually merge adjacent hotslices if close enough
  if (! is.null(minGapWidth) && is.numeric(minGapWidth) && minGapWidth >= 1L){
    hotslices.wt <- GenomicRanges::reduce(hotslices.wt, min.gapwidth=minGapWidth, with.revmap=FALSE)
  }
  
  # keep only regions that are sufficiently wide
  hotslices.wt <- hotslices.wt[which(width(hotslices.wt) >= minWidth)]
  # I do not care about proximity to target borders here because we have our minimal coverage in WT guaranteed
  
  # check corresponding region at MUT for minimal read coverage there in toto
  hotslices.wt <- hotslices.wt[which(IRanges::countOverlaps(query = hotslices.wt, subject = rga.mut) / width(hotslices.wt) * READL >= minAvgCov_MUT)]
  
  
  hotslices.wt
}




#' Cigar strings from a mapping BAM-file are transformed into a dataframe of positions and clipped base information (kind of pileup format for clipped reads).
#' 
#' This is a helper function used in the clipped peak filter.
#' @param pos vector of left-most position of mapped (non-clipped) base per read segment (from BAM)
#' @param cigar cigar string for each read segment. ##rLen read segment length
#' @return dataframe with genomic position and number of clipped bases, left (5') and right (3'). Can have zero rows if now clipped base is in the cigar strings.
bamReadToClipPos <- function(pos, cigar){  #, nbrS5p, nbrS3p
  
  # sum of clipped bases per given read, parsed from CIGAR string, separate for either ends (5p and 3p)
  nbrS5p <- as.integer(stringr::str_sub(stringr::str_extract(cigar, pattern = "^[[:digit:]]+S"), end=-2L)) #5' = left-most clipped bases
  nbrS3p <- as.integer(stringr::str_sub(stringr::str_extract(cigar, pattern = "[[:digit:]]+S$"), end=-2L)) #3'= right-most clipped bases
  
  nbrS5p[is.na(nbrS5p)] <- 0L
  nbrS3p[is.na(nbrS3p)] <- 0L
  
  stopifnot( length(pos) == length(nbrS5p) )
  stopifnot( length(pos) == length(nbrS3p) )
  
  # read length directly from cigar string is problematic because of I and D
  # I: read has additional bases that are not in ref. For ref-coordinates: (does not count for G-coord)
  # D/N (delete/skip on ref): read misses bases that are there in ref: + D (count up for G-count)
  # P (pad) is neutral (does not count for G-coord)
  
  # length in between the two (possibly clipped) ends of read segment
  # to get genomic coordinates for read bases, count M(/X/=)/D/N (not S or H not I nor P)
  rLenInBetween <- sapply(stringr::str_extract_all(cigar, pattern = "[[:digit:]]+[MX=DN]"),
                          function(x) sum(as.numeric(stringr::str_extract_all(x, pattern = "[[:digit:]]+"))))
  
  # prepare vector of positions for each soft-clipped base, 5p and 3p separately, to be tabulated later
  if (sum(nbrS5p) > 0L){
    S5bases <- integer(sum(nbrS5p))
    
    # 5p S: clipped bases to the left of read
    for (i in which(nbrS5p > 0L)){
      ind.S5 <- which.max(S5bases == 0L)
      S5bases[ind.S5:(ind.S5+nbrS5p[i]-1L)] <- pos[i] - (nbrS5p[i]:1L) # go to left
      # quality of 5' S bases (left-most bases in BAM) as weight?
    }
    
    S5tab <- table(S5bases)
    S5df <- as.data.frame(S5tab, responseName = "Sbases_5p")
    names(S5df)[1] <- "pos"
    S5df$pos <- as.integer(levels(S5df$pos))[S5df$pos]
  } else { 
    S5df <- data.frame(pos=integer(0), Sbases_5p=integer(0))
  }  
    
  if (sum(nbrS3p) > 0L ){  
    S3bases <- integer(sum(nbrS3p))
    
    # 3p S: clipped bases to the right of read (3' prime end)
    for (j in which(nbrS3p > 0L)){
      ind.S3 <- which.max(S3bases == 0L)
      #pos[j] - nbrS5p[j] = start of read
      S3bases[ind.S3:(ind.S3+nbrS3p[j]-1L)] <- pos[j] + rLenInBetween[j] + (1L:nbrS3p[j])  # go to right
      # quality of 3' S-bases (right most bases in BAM) as weight?
    }
    
    S3tab <- table(S3bases)
    S3df <- as.data.frame(S3tab, responseName = "Sbases_3p")
    names(S3df)[1] <- "pos"
    S3df$pos <- as.integer(levels(S3df$pos))[S3df$pos]
  } else {
    S3df <- data.frame(pos=integer(0), Sbases_3p=integer(0))
  }
    
  S53df <- dplyr::full_join(S5df, S3df, by="pos")
  S53df$Sbases_5p[is.na(S53df$Sbases_5p)] <- 0L
  S53df$Sbases_3p[is.na(S53df$Sbases_3p)] <- 0L
  
      
  # S5bases <- integer(sum(nbrS5p))
  # S3bases <- integer(sum(nbrS3p))
  # 
  # # 5p S: clipped bases to the left of read
  # for (i in which(nbrS5p > 0L)){
  #   ind.S5 <- which.max(S5bases == 0L)
  #   S5bases[ind.S5:(ind.S5+nbrS5p[i]-1L)] <- pos[i] - (nbrS5p[i]:1L) # go to left
  #   # quality of 5' S bases (left-most bases in BAM) as weight?
  # }
  # 
  # # 3p S: clipped bases to the right of read (3' prime end)
  # for (j in which(nbrS3p > 0L)){
  #   ind.S3 <- which.max(S3bases == 0L)
  #   #pos[j] - nbrS5p[j] = start of read
  #   S3bases[ind.S3:(ind.S3+nbrS3p[j]-1L)] <- pos[j] + rLenInBetween[j] + (1L:nbrS3p[j])  # go to right
  #   # quality of 3' S-bases (right most bases in BAM) as weight?
  # }
  # 
  # 
  # S5tab <- table(S5bases)
  # S5df <- as.data.frame(S5tab, responseName = "Sbases_5p")
  # names(S5df)[1] <- "pos"
  # S5df$pos <- as.integer(levels(S5df$pos))[S5df$pos]
  # 
  # S3tab <- table(S3bases)
  # S3df <- as.data.frame(S3tab, responseName = "Sbases_3p")
  # names(S3df)[1] <- "pos"
  # S3df$pos <- as.integer(levels(S3df$pos))[S3df$pos]
  # 
  # #return( union(left_join(S5df, S3df, by="pos"), right_join(S5df, S3df, by="pos")) )
  # S53df <- dplyr::full_join(S5df, S3df, by="pos")
  # S53df$Sbases_5p[is.na(S53df$Sbases_5p)] <- 0L
  # S53df$Sbases_3p[is.na(S53df$Sbases_3p)] <- 0L
  
  # mkuhn, 20141217: every row should be genomic positions with clipped bases
  stopifnot( all( S53df$Sbases_5p + S53df$Sbases_3p > 0L) )
  
  return( S53df )
}



#' Counts the number of clipped bases per genomic position.
#' 
#' This is a helper function to find clipped peaks. The clipped proportion estimate is biased downwards by using cov+1 in denumerator to avoid cov=0 problems.
#' @return dataframe with Q1-coverage and clipped base information per genomic coordinate of the target region
#' @export
getClippedBasePositions <- function(bamFile, chrom=CHROM, targetStartPos, targetEndPos){
  
  stopifnot( length(chrom) == 1L )
  logdebug("Getting clipped base information for region %s:%d-%d.", chrom, targetStartPos, targetEndPos)
  
  
  if ( ! checkFile(bamFile) ){
    logerror("BAM file %s not found!", bamFile)
    return(invisible(FALSE))
  }#fi 
  
  
  myParam <- Rsamtools::ScanBamParam(flag=Rsamtools::scanBamFlag(isPaired=TRUE, isUnmappedQuery = FALSE, isProperPair = NA, isSecondaryAlignment = FALSE),
                                     what=setdiff(Rsamtools::scanBamWhat(), c("qwidth", "seq", "qual")),  #"qual", #"seq",
                                     which=GenomicRanges::GRanges(seqnames=chrom, ranges=IRanges::IRanges(start=targetStartPos, end=targetEndPos)) )
  
  # mapping filter on Q1 reads
  mapping <- data.frame(Rsamtools::scanBam(file=bamFile, param=myParam)[[1]], stringsAsFactors = FALSE) %>% 
    dplyr::filter_(~mapq >= 1L)
  
  ## sort on position not necessary
  #mapping <- dplyr::arrange_(mapping, ~pos)
  
  # it is a sorted BAM file
  stopifnot( all( diff(mapping$pos) >= 0 ) )
  
  # no hard-clipping  (see BWA MEM -Y option)
  stopifnot( ! any(grepl("H", mapping$cigar)) )
  
  
  # # sum of clipped bases within the Q1-reads, read from CIGAR string, separate for either ends, 5p and 3p
  # mapping$Sbases_5p <- as.integer(stringr::str_sub(stringr::str_extract(mapping$cigar, pattern = "^[[:digit:]]+S"), end=-2L))
  # mapping$Sbases_3p <- as.integer(stringr::str_sub(stringr::str_extract(mapping$cigar, pattern = "[[:digit:]]+S$"), end=-2L))
  # 
  # mapping$Sbases_5p[is.na(mapping$Sbases_5p)] <- 0L
  # mapping$Sbases_3p[is.na(mapping$Sbases_3p)] <- 0L
  
  # logdebug("There are ? 5p (start) and ? 3p (end) S-bases in bamfile ::%s:: at %s at %d - %d.", 
  #         #sum(mapping$Sbases_5p, na.rm=TRUE), sum(mapping$Sbases_3p, na.rm=TRUE),
  #         basename(bamFile), chrom, targetStartPos, targetEndPos)
  
  
  
  # mkuhn, 2014-12-16: base-coordinate centric counts of clipped bases (and base quality for weight?)
  #mappingClip <- filter(mapping, Sbases_5p > 0 | Sbases_3p > 0)
  clippedBasePosition0 <- bamReadToClipPos(pos=mapping$pos, cigar=mapping$cigar) #, nbrS5p=mapping$Sbases_5p, nbrS3p=mapping$Sbases_3p)
  logdebug("In region %s:%s-%s found %d genomic positions with clipped bases.", chrom, targetStartPos, targetEndPos, NROW(clippedBasePosition0))
  
  
  ##
  ## Q1-Coverage per genomic base
  ##
  TARGET_REGION <- paste0(chrom, ":", targetStartPos, "-", targetEndPos)
  cov.samtools <- system(paste(SAMTOOLS_EXE, "depth -Q1 -r ",TARGET_REGION, bamFile), intern=TRUE)
  
  # zero-coverage matrix for whole region
  zeroCov <- matrix(data = c(seq(targetStartPos, targetEndPos), rep(0L, targetEndPos - targetStartPos + 1L)),
                    ncol=2L, dimnames = list(NULL, c("pos", "readCov")))
  
  # mkuhn, 2015-11-29 don't stop when no Q1-reads but return all 0s
  #stopifnot( length(cov.samtools) >= 1L )
  if ( length(cov.samtools) >= 1L ){
    
    # ZZZ assuming a single fixed chromosome for the region
    stopifnot(  length(chrom) == 1L && identical(chrom, unique(stringr::str_replace(cov.samtools, "\t.*", ""))) )  
    
    # keep position and coverage, but fill zero cov positions
    cov.raw <- matrix(as.numeric(unlist(strsplit(stringr::str_replace(cov.samtools, paste0("^", chrom, "\t"), replacement = ""), split="\t", fixed=TRUE))),
                      ncol=2L, byrow=TRUE)
    # fill up coverage with 0s where no coverage
    cov.mat <- rbind(cov.raw, zeroCov[! zeroCov[,1L] %in% cov.raw[,1L],])
    cov.mat <- cov.mat[order(cov.mat[,1]),]
    
    # not here, I do it below
    # cov.mat[,2L] <- cov.mat[,2L]+1L # avoid 0
    # stopifnot( all(cov.mat[,2L] > 0L) )
    colnames(cov.mat) <- c("pos", "readCov")
  } else {
    logwarn("No Q1-reads at region %s in mapping %s.", TARGET_REGION, bamFile)
    cov.mat <- zeroCov
  }
  
  # keep only all positions of the target region
  clippedBasePosition <- data.frame(cov.mat) %>% 
    dplyr::left_join(clippedBasePosition0, by="pos")  #dplyr::full_join
  clippedBasePosition$Sbases_5p[is.na(clippedBasePosition$Sbases_5p)] <- 0L
  clippedBasePosition$Sbases_3p[is.na(clippedBasePosition$Sbases_3p)] <- 0L
  #  clipped reads are not counted for the coverage `readCov`
  #  [mkuhn, 2015-11-21] I do not want this, because it dilutes the clipped peak signal
  #  [mkuhn, 2015-12-07] I turned it back on because on real data we get frequent peaks at 100% otherwise
  clippedBasePosition$readCov <- clippedBasePosition$readCov + clippedBasePosition$Sbases_5p + clippedBasePosition$Sbases_3p
  
  # relative proportion of clipped bases
  clippedBasePosition %<>% 
    dplyr::mutate_(relSbases_5p = ~pmin(1L, Sbases_5p / (readCov + 1L)),
                   relSbases_3p = ~pmin(1L, Sbases_3p / (readCov + 1L)))  %>% 
    dplyr::arrange_(~pos)
  
  return(clippedBasePosition)
}


#' Find peaks of clipped bases.
#' 
#' Smooth the proportion of clipped bases at given region.
#' Smoothing uses moving average with Hamming weights.
#' 
#' For threshold, a minimum threshold is given as parameter but also quantile values of smoothed clipped-base proportions are used.
#' An extra effort is done to rescue peaks at the borders of the region.
#' 
#' @param clippedBasePos dataframe with number of coverage and number of reads with clipped bases per position
#' @param windowSize integer size of Hamming window for smoothing
#' @note Implicitly a clipped base gets more weight on low coverage regions (because it can more easy go up to 100\%). Should we use weighted mean here?
#' @return list with peak-dataframe and specifically number of 5p and 3p peaks
findClippedPeaks <- function(clippedBasePos, windowSize=31L, minThresholdAbs=2L, minThreshold=0.05, thresholdQuantile=0.95){ #}, minTargetBorderDist=1L){
  
  # have minimal odd windowsize
  windowSize <- max(7L, windowSize)
  if (windowSize %% 2L == 0L) windowSize <- windowSize + 1L
  
  
  # normal smoothing window (for region interior) -----
  
  # start smoothing via filters
  hammWeights <- getHammingWindowWeights(windowSize)
  clippedBasePosF <- data.frame(pos=as.integer(round(RcppRoll::roll_mean(clippedBasePos$pos, n=windowSize, weights = hammWeights))),
                                readCov=RcppRoll::roll_mean(clippedBasePos$readCov, n=windowSize, weights = hammWeights), #cov can be useful to weight the relS-proportions
                                #Sbases=RcppRoll::roll_mean(clippedBasePos$Sbases_5p + clippedBasePos$Sbases_3p, n=windowSize, weights = hammWeights),
                                Sbases_5p=RcppRoll::roll_mean(clippedBasePos$Sbases_5p, n=windowSize, weights = hammWeights),
                                Sbases_3p=RcppRoll::roll_mean(clippedBasePos$Sbases_3p, n=windowSize, weights = hammWeights),
                                relSbases_5p=RcppRoll::roll_mean(clippedBasePos$relSbases_5p, n=windowSize, weights = hammWeights),
                                relSbases_3p=RcppRoll::roll_mean(clippedBasePos$relSbases_3p, n=windowSize, weights = hammWeights)
  )
  
  
  # indices (in clippedBasePosF) of local maxima that reach highest values (above given quantile)
  quant_rel5p <- quantile(clippedBasePosF$relSbases_5p, probs=thresholdQuantile, na.rm=TRUE)
  quant_rel3p <- quantile(clippedBasePosF$relSbases_3p, probs=thresholdQuantile, na.rm=TRUE)
  
  empMinThreshold_5p <- max(minThreshold, quant_rel5p, quant_rel3p * .33)
  empMinThreshold_3p <- max(minThreshold, quant_rel3p, quant_rel5p * .33)
  empMinThreshold <- max(minThreshold, quant_rel5p, quant_rel3p)
  
  localMaxInd.5p <- localMax(clippedBasePosF$relSbases_5p, threshold = empMinThreshold_5p )
  localMaxInd.3p <- localMax(clippedBasePosF$relSbases_3p, threshold = empMinThreshold_3p )
  
  # mkuhn, 2015-12-21: I do not care if clipped peaks is pretty much close to the border.. Why should I?!
  # no maximum at the borders (index-wise)
  # localMaxInd.5p <- localMaxInd.5p[! localMaxInd.5p %in% c(1L, NROW(clippedBasePosF))]
  # localMaxInd.3p <- localMaxInd.3p[! localMaxInd.3p %in% c(1L, NROW(clippedBasePosF))]
  
  # require minimal absolute number of clipped bases, 3p and 5p
  localMaxInd.5p <- intersect(x=localMaxInd.5p, y=which(clippedBasePosF$Sbases_5p >= minThresholdAbs))
  localMaxInd.3p <- intersect(x=localMaxInd.3p, y=which(clippedBasePosF$Sbases_3p >= minThresholdAbs))
  
  
  
  # fine smoothing at region borders ----
  
  # rescue peaks close to the border
  windowSizeXS <- floor(max(3L, windowSize/3L))
  hammWeightsXS <- getHammingWindowWeights(n=windowSizeXS)
  
  clippedBasePosStart <- head(clippedBasePos, n=(windowSize-1L) / 2L)
  clippedBasePosEnd   <- tail(clippedBasePos, n=(windowSize-1L) / 2L)
  
  clippedBasePosStartF <- data.frame(pos=as.integer(round(RcppRoll::roll_mean(clippedBasePosStart$pos, n=windowSizeXS, weights = hammWeightsXS))),
                                     readCov=RcppRoll::roll_mean(clippedBasePosStart$readCov, n=windowSizeXS, weights = hammWeightsXS), #cov can be useful to weight the relS-proportions
                                     #Sbases=RcppRoll::roll_mean(clippedBasePos$Sbases_5p + clippedBasePos$Sbases_3p, n=windowSize, weights = hammWeights),
                                     Sbases_5p=RcppRoll::roll_mean(clippedBasePosStart$Sbases_5p, n=windowSizeXS, weights = hammWeightsXS),
                                     Sbases_3p=RcppRoll::roll_mean(clippedBasePosStart$Sbases_3p, n=windowSizeXS, weights = hammWeightsXS),
                                     relSbases_5p=RcppRoll::roll_mean(clippedBasePosStart$relSbases_5p, n=windowSizeXS, weights = hammWeightsXS),
                                     relSbases_3p=RcppRoll::roll_mean(clippedBasePosStart$relSbases_3p, n=windowSizeXS, weights = hammWeightsXS)
  )
  
  clippedBasePosEndF <- data.frame(pos=as.integer(round(RcppRoll::roll_mean(clippedBasePosEnd$pos, n=windowSizeXS, weights = hammWeightsXS))),
                                     readCov=RcppRoll::roll_mean(clippedBasePosEnd$readCov, n=windowSizeXS, weights = hammWeightsXS), #cov can be useful to weight the relS-proportions
                                     #Sbases=RcppRoll::roll_mean(clippedBasePos$Sbases_5p + clippedBasePos$Sbases_3p, n=windowSize, weights = hammWeights),
                                     Sbases_5p=RcppRoll::roll_mean(clippedBasePosEnd$Sbases_5p, n=windowSizeXS, weights = hammWeightsXS),
                                     Sbases_3p=RcppRoll::roll_mean(clippedBasePosEnd$Sbases_3p, n=windowSizeXS, weights = hammWeightsXS),
                                     relSbases_5p=RcppRoll::roll_mean(clippedBasePosEnd$relSbases_5p, n=windowSizeXS, weights = hammWeightsXS),
                                     relSbases_3p=RcppRoll::roll_mean(clippedBasePosEnd$relSbases_3p, n=windowSizeXS, weights = hammWeightsXS)
  )
  
  
  #max(clippedBasePosStartF$relSbases_5p[which(clippedBasePosStartF$relSbases_5p >= minThreshold)])
  BORDER_FACTOR <- max(1L, windowSize / windowSizeXS) # how much is my threshold higher at the borders than normal? (protect against false positives at the borders because it is )
  localMaxIndStart.5p <- localMax(clippedBasePosStartF$relSbases_5p, threshold = BORDER_FACTOR*empMinThreshold)
  localMaxIndStart.3p <- localMax(clippedBasePosStartF$relSbases_3p, threshold = BORDER_FACTOR*empMinThreshold)
  
  localMaxIndStart.5p <- intersect(x=localMaxIndStart.5p, y=which(clippedBasePosStartF$Sbases_5p >= minThresholdAbs))
  localMaxIndStart.3p <- intersect(x=localMaxIndStart.3p, y=which(clippedBasePosStartF$Sbases_3p >= minThresholdAbs))
  
  localMaxIndEnd.5p <- localMax(clippedBasePosEndF$relSbases_5p, threshold = BORDER_FACTOR*empMinThreshold)
  localMaxIndEnd.3p <- localMax(clippedBasePosEndF$relSbases_3p, threshold = BORDER_FACTOR*empMinThreshold)
  
  localMaxIndEnd.5p <- intersect(x=localMaxIndEnd.5p, y=which(clippedBasePosEndF$Sbases_5p >= minThresholdAbs))
  localMaxIndEnd.3p <- intersect(x=localMaxIndEnd.3p, y=which(clippedBasePosEndF$Sbases_3p >= minThresholdAbs))
  

  
  
  # combined peak candidates ------
  
  # 5' and 3' max positions
  # enforce minimal distance to target border (nucleotide-position-wise)
  # rank according to the mean proportion of clipped bases: highest proportion gets smallest ranks.
  # Ranks are normalized so that highest rank number is 1 that corresponds to lowest clipped-base proportion.
  localMax.5p <- rbind(clippedBasePosF[localMaxInd.5p, ],
                       clippedBasePosStartF[localMaxIndStart.5p,],
                       clippedBasePosEndF[localMaxIndEnd.5p,])
  
  localMax.3p <- rbind(clippedBasePosF[localMaxInd.3p,],
                       clippedBasePosStartF[localMaxIndStart.3p,],
                       clippedBasePosEndF[localMaxIndEnd.3p,])
  
  
  
  # mkuhn, for clipped base peaks do not look at distance to target border
  localMax.5p %<>%
    dplyr::select_( ~ -relSbases_3p) %>% 
    #dplyr::filter_(~pos > .targetStartPos + minTargetBorderDist, ~pos < .targetEndPos - minTargetBorderDist) %>% 
    dplyr::mutate_(rank = ~rank(-relSbases_5p) / n())
  
  localMax.3p %<>% 
    dplyr::select_( ~ -relSbases_5p) %>% 
    #dplyr::filter_(~pos > .targetStartPos + minTargetBorderDist, ~pos < .targetEndPos - minTargetBorderDist) %>% 
    dplyr::mutate_(rank = ~rank(-relSbases_3p) / n())
  
  
  # Idea: every clipped base peak comes as 5p and as corresponding 3p-peak. Their distance can vary depending on SV type. Some SVs have 2 peaks, others only 1.
  # match 5' and 3' max positions based on peak height and distance. Peak height is more important?!?
  localMax.53 <- merge(localMax.5p, localMax.3p, by=NULL, suffixes = c(".5p", ".3p")) #cartesian product
  
  # score rates every pair of 3p-5p clipped bases peaks: the smaller the score the better the match
  localMax.53 %<>% 
    dplyr::mutate_(posDiff= ~pos.3p - pos.5p, ## how far (in nt) lies 5p peak ahead of 3p
                   # rank.posDiff a score that assesses the 5' and 3'-peak distance.
                   # 1st factor: Small differences get small rank, hence lower 1st factor. sqrt() for posDiff rank gives less weight to pos diff
                   # 2nd factor: give extra benefit to distances below 3 * READL.  pmin(1, pmax(0.5, sqrt(x/50))) lives betw 0.5 and 1. 
                   #+ Simpler step-function version is:  * ifelse( abs(localMax.53$posDiff) <= 50L, 0.5, 1) 
                   rank.posDiff= ~sqrt(rank(abs(posDiff))/n()) * pmin(1L, pmax(0.5, sqrt(abs(posDiff) / (1L + 3 * READL)))),
                   score= ~rank.5p * rank.3p * rank.posDiff) %>% 
    dplyr::arrange_(~score) %>% # best hits first
    dplyr::filter_(~ ! duplicated(pos.5p), ~ ! duplicated(pos.3p)) # keep only combinations with best scores, i.e. drop entries as soon as a coordinate (pos.5p or pos.3p) has been used before, with a better (=lower) score
  
  
  list(clippedBases=clippedBasePosF, peaks=localMax.53, peaks5p=localMax.5p, peaks3p=localMax.3p) #nbrPeaks5p=NROW(localMax.5p), nbrPeaks3p=NROW(localMax.3p))
}


#' Assess covered target regions with respect to clipped bases.
#' 
#' Idea is to pinpoint regions that carry many clipped bases as candidate regions for a SV.
#' False positives (i.e. \emph{no} SV in the target region) are not so much of a concern here because it is a screening step.
#' The assumption is that features extracted from mapping in the vicinity of clipped-base peaks only will be much more sensitive than the features from the whole target region.
#' 
#' The function uses the genome coordinates as reference point. The target regions typically come from an enrichment.
#' 
#' This function ignores coverage and should be run on well-covered regions only.
#' ZZZ we only find one peak or two peaks if they are very close (measured via score).
#' But there can be more (different) peaks in a covered region. Idea could be to identify peaks based on being a positive outlier for prop clipped bases.
#' Idea: learn from data to better assess clipped 3p/5p peak properties. What is "normal" property distribution depending on SV (present or not/type) and coverage
#' ZZZ weighting of clipped base proportion based on coverage not implemented! High proportion at low coverage should not count much, no? (High coverage should count more, no?!)
#' ZZZ a peak is defined for 3p *and* 5p clipped bases. At border to uncovered regions one peak should be enough..
#' 
#' @param bamFile mapping filename. Typically a tumor BAM file.
#' @param minThresholdAbs numeric. minimal number of clipped base needed
#' @param minThreshold numeric. minimal proportion of clipped base needed
#' @param thresholdQuantile numeric < 1. quantile of clipped base proportion distribution in the region that needs to be exceeded
#' @return list with elements vector clippedInfo, a named vector peakInfo (which is a subset of clippedInfo) and a vector nonPeakPos of negative candidate positions
#' @export
assessClippedBasesInRegion <- function(bamFile, chrom=CHROM, .targetStartPos, .targetEndPos, do.plot=FALSE, 
                                        windowSize=31L, minThresholdAbs=2L, minThreshold=0.05, thresholdQuantile=0.95){ #, minTargetBorderDist=1L){
  if (FALSE){
    # Example: simulated Pat1, Tumor sample w/ differential insertion of length 255bp
    bamFile <- paste0(BASEDIR, "mapping/cov60_TL50/PEm200sd25_R100/pat00001_T.bam")
    chrom <- "chr5"; .targetStartPos <- 79365843; .targetEndPos <- 79368003
    windowSize=31L; minThresholdAbs=2L; minThreshold=0.05; thresholdQuantile=0.95 #; minTargetBorderDist=1L #READL/2L
    # example case with inversion
    #windowSize <- 30L; thresholdQuantile=0.9; minTargetBorderDist=35L; chrom="chr5"
    #bamFile  <- "/Users/mkuhn/seq1/SVdata/mapping/cov60_TL90/PEm200sd25_R100/pat01020_T.bam" 
    #.targetStartPos <- 149512443; .targetEndPos <- 149514571
    #SVlen <- 197; SVoffset <- 1750
    #bamFile  <- "/Users/mkuhn/seq1/SVdata/mapping/cov60_TL90/PEm200sd25_R100/pat01000_T.bam"
    #.targetStartPos <- 10410534; .targetEndPos <- 10410737
  }
  
  stopifnot( length(chrom) == 1L )
  
  stopifnot( checkFile(bamFile) )
  
  
  # get number of clipped bases per genomic location
  clippedBasePos <- getClippedBasePositions(bamFile, chrom, targetStartPos=.targetStartPos, targetEndPos=.targetEndPos)
  stopifnot( all(clippedBasePos$pos >= .targetStartPos), all(clippedBasePos$pos <= .targetEndPos) )
  
  
  clippedPeakInfo <- findClippedPeaks(clippedBasePos, windowSize = windowSize, minThresholdAbs=minThresholdAbs, minThreshold=minThreshold, thresholdQuantile=thresholdQuantile)
  # including peaks close to the target region border
  localMax.53 <- clippedPeakInfo[["peaks"]]
  peaks5p <- clippedPeakInfo[["peaks5p"]] 
  peaks3p <- clippedPeakInfo[["peaks3p"]]
  clippedBasePosF <- clippedPeakInfo[["clippedBases"]]
  # hit consists of 5' and 3' peaks
  # inversions give two hits: but then their positions should not be very far apart..
  
  # summary information of clipped bases for target region
  clippedInfo <- c(regionStart=.targetStartPos, regionEnd=.targetEndPos,
                   minThreshold=minThreshold,
                   nbrLocalMax.5p=NROW(peaks5p), nbrLocalMax.3p=NROW(peaks3p),
                   # mkuhn, 2015-11-27: try to extract peak information from clipping information dataframe. These few numbers are not enough to represent the whole SV information
                   # posDiff, pos.5p, pos.3p only of peak with best (=lowest) score
                   # ZZZ Hence, abstraction from peaks to SV necessary
                   nbrPeaks=NROW(localMax.53),
                   minScore=if (NROW(localMax.53) >= 1L) min(localMax.53$score) else NA_real_,
                   posDiff=if (NROW(localMax.53) >= 1L) localMax.53$posDiff[which.min(localMax.53$score)] else NA_integer_,
                   pos.5p=if (NROW(localMax.53) >= 1L) localMax.53$pos.5p[which.min(localMax.53$score)] else NA_integer_,
                   pos.3p=if (NROW(localMax.53) >= 1L) localMax.53$pos.3p[which.min(localMax.53$score)] else NA_integer_,
                   pos.left=if (NROW(localMax.53) >= 1L) min(localMax.53$pos.3p, localMax.53$pos.5p) else NA_integer_,
                   pos.right=if (NROW(localMax.53) >= 1L) max(localMax.53$pos.3p, localMax.53$pos.5p) else NA_integer_,
                   relS_5p=quantile(clippedBasePosF$relSbases_5p, probs=c(.5, .75, .90, .95, .99)),
                   sd_5p=sd(clippedBasePosF$relSbases_5p),
                   relS_3p=quantile(clippedBasePosF$relSbases_3p, probs=c(.5, .75, .90, .95, .99)),
                   sd_3p=sd(clippedBasePosF$relSbases_3p)
                   #median_5p=median(clippedBasePosF$relSbases_5p), quant75_5p=quantile(clippedBasePosF$relSbases_5p, .75, names=FALSE), quant90_5p=quantile(clippedBasePosF$relSbases_5p, .90), quant95_5p=quantile(clippedBasePosF$relSbases_5p, .95), quant99_5p=quantile(clippedBasePosF$relSbases_5p, .99), mean_5p=mean(clippedBasePosF$relSbases_5p), sd_5p=sd(clippedBasePosF$relSbases_5p),
                   #median_3p=median(clippedBasePosF$relSbases_3p), quant75_3p=quantile(clippedBasePosF$relSbases_3p, .75), quant90_3p=quantile(clippedBasePosF$relSbases_3p, .90), quant95_3p=quantile(clippedBasePosF$relSbases_3p, .95), quant99_3p=quantile(clippedBasePosF$relSbases_3p, .99), mean_3p=mean(clippedBasePosF$relSbases_3p), sd_3p=sd(clippedBasePosF$relSbases_3p)
  )
  # mkuhn, 2016-01-04: remove %-sign at end of quantile-names
  names(clippedInfo) <- sub("%$", "", names(clippedInfo))
  
  
  
  
  # mkuhn, 2015-11-21: "negative candidates"
  # find loci with no peaks and not in vicinity of local maxima of either 5p or 3p
  # operate on indices of clippedBasePosF-data
  
  # mkuhn, 2015-12-22: old implementation used the index. this is problematic since I do now peak rescue at borders of target region and that can shift the indices
  #+ better to use position directly 
  # neg.cand.df <- data.frame(ind=which(clippedBasePosF$relSbases_5p <= min(minThreshold, quantile(clippedBasePosF$relSbases_5p, probs=.33)) &
  #                                       clippedBasePosF$relSbases_3p <= min(minThreshold, quantile(clippedBasePosF$relSbases_3p, probs=.33))))
  # neg.tabu.ind <- c(1L:minTargetBorderDist, localMaxInd.5p, localMaxInd.3p, (NROW(clippedBasePosF)-minTargetBorderDist):NROW(clippedBasePosF))
  # neg.cand.df$minDistToLocalMax <- apply(neg.cand.df, 1L, function(i) min(abs(i-neg.tabu.ind)))
  # # negative candidates should have at least a distance of minTargetBorderDist to the local clipped peak maxima that go beyond the quantile and the absolute threshold
  # neg.cand.pos <- clippedBasePosF$pos[neg.cand.df$ind[localMax(neg.cand.df$minDistToLocalMax, threshold = minTargetBorderDist)]]
  
  
  neg.pos.cand <- clippedBasePosF$pos[which(clippedBasePosF$relSbases_5p <= min(minThreshold, quantile(clippedBasePosF$relSbases_5p, probs=.33)) &
          clippedBasePosF$relSbases_3p <= min(minThreshold, quantile(clippedBasePosF$relSbases_3p, probs=.33)))]
  
  #neg.pos.tabu <- c(.targetStartPos:(.targetStartPos+minTargetBorderDist), peaks5p$pos, peaks3p$pos, (.targetEndPos-minTargetBorderDist):.targetEndPos)
  neg.pos.tabu <- c(.targetStartPos, peaks5p$pos, peaks3p$pos, .targetEndPos)
  neg.pos <- neg.pos.cand[localMax(sapply(neg.pos.cand, function(p) min(abs(p-neg.pos.tabu))), threshold = 3L )] #threshold=minTargetBorderDist
  
  
  
  
  
  # mkuhn, 2015-03-03: if running on small region it can be NULL
  #stopifnot( NROW(localMax.53) >= 1L )
  
  
  # mkuhn, 2015-11-27: old implementation    
  # # ZZZ assume at most one SV in target region
  # # decide if we return two peak-combinations or only one peak-combination
  # peakInfo <- NULL
  # 
  # # return two peaks if they are closer together than average (according to combined score)
  # peakInfo <- if ( NROW(localMax.53) >= 2L ) #&& (localMax.53$score[2L] - localMax.53$score[1L]) <= avgScoreDiff )
  #   localMax.53[1:2,] else 
  #     if (NROW(localMax.53) >= 1L) localMax.53[1L,]
  # # else  if (NROW(localMax.53) == 1L) # used to give NULL in case more than 1 entry that is not close together
  peakInfo <- clippedInfo[c("nbrPeaks", "posDiff", "pos.5p", "pos.3p", "relS_5p.95", "relS_3p.95", "relS_5p.99", "relS_3p.99")]
  
  # 
  # Plotting
  if (isTRUE(do.plot)){
    plot( relSbases_5p ~ pos, data=clippedBasePosF, type="l", col="darkblue", ylab="prop. clipped bases",
          ylim=c(0, 1.05*max(c(clippedBasePosF$relSbases_5p, clippedBasePosF$relSbases_3p))) )
    lines(relSbases_3p ~ pos, data=clippedBasePosF, type="l", lty=2, lwd=1.5, col="green")
    rug(x=peaks5p$pos, side=3, col="darkblue", lwd=2)
    rug(x=peaks3p$pos, side=3, col="green", ticksize = -.02, lwd=2)
    abline(h=minThreshold, lty=3, col="lightgrey")
    legend("topright", legend = c("5'", "3'"), col = c("darkblue", "green"), lwd=2, lty=1:2, inset = c(.02,.03))
  }
  
  
  #clippedBasePos=clippedBasePos, 
  return( list(clippedInfo=clippedInfo, peakInfo=peakInfo, nonPeakPos=neg.pos) )  
}




#' Assess covered (enrichment) target regions with respect to clipped bases. OLD: use now assessClippedBasesInRegion
#' 
#' The aim is to decide if a given covered target region is positive or negative for clipped base enrichment.
#' This function aims to bring focus within a bigger target region. 
#' Idea is to pinpoint regions that carry many clipped bases as candidate regions for a SV.
#' False positives (i.e. \emph{no} SV in the target region) are not a concern here because it is a screening step.
#' The hope is that by doing so the features from mapping in the vicinity only will be much more sensitive than the features from the whole target region.
#' Still, it is unclear how do I know how big I need to make the candidate region? This depends on the SV-length, no? Can I see that from clipped reads profile?
#' 
#' The function uses the genome coordinates as reference point, hence the _G for genomic. An older version used the BAM coordinates which is not precise.
#' 
#' This function ignores coverage and should be run on well-covered regions only.
#' ZZZ we only find one peak or two peaks if they are very close (measured via score).
#' But there can be more (different) peaks in a covered region. Idea could be to identify peaks based on being a positive outlier for prop clipped bases.
#' ZZZ: learn from data: assess clipped 3p/5p peak properties depending on SV (present or not/type) and coverage
#' ZZZ should I use negCandPos only for regions with no strong clipped-peak at all?
#' 
#' @param a BAM file. Typically, this will be a mapping from a tumor sample.
#' @param targetStartPos single integer. Start position (5p) of the target region.
#' @param targetEndPos single integer. End position (3p) of the target region.
#' @param do.plot logical flag if a clipped-bases plot per location should be done
#' @param windowSize width of the Hamming window for averaging.
#' @param minThreshold minimal precentage of smoothed 3p- resp 5p-clipped bases to count as positive candidate
#' @param thresholdQuantile quantile value that a peak in clipped 5p-bases should pass to be retained as local maximum
#' @param minTargetBorderDist minimal distance to border of target region
#' @return List of length 3. First, a dataframe with per genomic nucleotide clipped base info. Second, a dataframe with peak locations, either one or two entries. Third, a vector with positions of possible non-hits.
#' @note This will be important for real data where I do \emph{not} know the location of an SV. And I need to make a good guess if I do not want to simply take the whole target region.
#' @export
findClippedPeaks_G_OLD <- function(bamFile, chrom="chr5", .targetStartPos, .targetEndPos, do.plot=FALSE, 
                               windowSize=31L, minThreshold=0.05, thresholdQuantile=0.95, minTargetBorderDist=READL/2L){
  if (FALSE){
    # Example: simulated Pat1, Tumor sample w/ differential insertion of length 255bp
    bamFile <- paste0(BASEDIR, "mapping/cov60_TL50/PEm200sd25_R100/pat00001_T.bam")
    chrom <- "chr5"; .targetStartPos <- 79365843; .targetEndPos <- 79368003
    windowSize=31L; minThreshold=0.05; thresholdQuantile=0.95; minTargetBorderDist=10L
    # example case with inversion
    #windowSize <- 30L; thresholdQuantile=0.9; minTargetBorderDist=35L; chrom="chr5"
    #bamFile  <- "/Users/mkuhn/seq1/SVdata/mapping/cov60_TL90/PEm200sd25_R100/pat01020_T.bam" 
    #.targetStartPos <- 149512443; .targetEndPos <- 149514571
    #SVlen <- 197; SVoffset <- 1750
    #bamFile  <- "/Users/mkuhn/seq1/SVdata/mapping/cov60_TL90/PEm200sd25_R100/pat01000_T.bam"
    #.targetStartPos <- 10410534; .targetEndPos <- 10410737
  }
  
  stopifnot( length(chrom) == 1L )
  
  if ( ! checkFile(bamFile) ){
    logerror("BAM file %s not found!", bamFile)
    return(invisible(FALSE))
  }#fi 
  
  
  # get number of clipped bases per genomic location
  clippedBasePos <- getClippedBasePositions(bamFile, chrom, targetStartPos=.targetStartPos, targetEndPos=.targetEndPos)
  
  stopifnot( all(clippedBasePos$pos >= .targetStartPos), all(clippedBasePos$pos <= .targetEndPos))
  
  # with(clippedBasePos, matplot(x=pos, y=cbind(Sbases_5p, Sbases_3p), pch=16, cex=0.7))
  # with(clippedBasePos, matplot(x=pos, y=sqrt(cbind(Sbases_5p, Sbases_3p)), pch=16, cex=0.7))
  # # normalized
  # with(clippedBasePos, matplot(x=pos, y=cbind(Sbases_5p/(cov+1L), Sbases_3p/(cov+1L)), pch=16, cex=0.7))
  
  
  
  # do filtering of the proportion of clipped bases above genomic bases 
  # filtered cumulated mapping data (with Hamming windows)
  hammWeights <- getHammingWindowWeights(windowSize)
  clippedBasePosF <- data.frame(pos=as.integer(round(RcppRoll::roll_mean(clippedBasePos$pos, n=windowSize, weights = hammWeights))),
                                relSbases_5p=RcppRoll::roll_mean(clippedBasePos$relSbases_5p, n=windowSize, weights = hammWeights),
                                relSbases_3p=RcppRoll::roll_mean(clippedBasePos$relSbases_3p, n=windowSize, weights = hammWeights)
  )
  
  
  # indices of local maxima that reach highest values (above given quantile)
  localMaxInd.5p <- localMax(clippedBasePosF$relSbases_5p, threshold = max(minThreshold, quantile(clippedBasePosF$relSbases_5p, probs=thresholdQuantile)))
  localMaxInd.3p <- localMax(clippedBasePosF$relSbases_3p, threshold = max(minThreshold, quantile(clippedBasePosF$relSbases_3p, probs=thresholdQuantile)))
  
  # no maximum at the borders (index-wise)
  localMaxInd.5p <- localMaxInd.5p[! localMaxInd.5p %in% c(1L, NROW(clippedBasePosF))]
  localMaxInd.3p <- localMaxInd.3p[! localMaxInd.3p %in% c(1L, NROW(clippedBasePosF))]
  
  # 5' and 3' max positions
  # enforce minimal distance to target border (nucleotide-position-wise)
  # rank according to the mean proportion of clipped bases: highest proportion gets smallest ranks.
  # Ranks are normalized so that highst rank number is 1 that corresponds to lowest clipped-base proportion.
  localMax.5p <- clippedBasePosF[localMaxInd.5p,] %>%
    dplyr::select_( ~ -relSbases_3p) %>% 
    dplyr::filter_(~pos > .targetStartPos + minTargetBorderDist, ~pos < .targetEndPos - minTargetBorderDist) %>% 
    dplyr::mutate_(rank = ~rank(-relSbases_5p) / n())
  
  localMax.3p <- clippedBasePosF[localMaxInd.3p,] %>% 
    dplyr::select_( ~ -relSbases_5p) %>% 
    dplyr::filter_(~pos > .targetStartPos + minTargetBorderDist, ~pos < .targetEndPos - minTargetBorderDist) %>% 
    dplyr::mutate_(rank = ~rank(-relSbases_3p) / n())
  
  # mkuhn, 2015-11-21: "negative candidates"
  # find loci with no peaks and not in vicinity of local maxima of either 5p or 3p
  # operate on indices of clippedBasePosF-data
  neg.tabu.ind <- c(1L:minTargetBorderDist, localMaxInd.5p, localMaxInd.3p, (NROW(clippedBasePosF)-minTargetBorderDist):NROW(clippedBasePosF))
  neg.cand.df <- data.frame(ind=which(clippedBasePosF$relSbases_5p <= min(minThreshold, quantile(clippedBasePosF$relSbases_5p, probs=.33)) &
                                        clippedBasePosF$relSbases_3p <= min(minThreshold, quantile(clippedBasePosF$relSbases_3p, probs=.33))))
  neg.cand.df$minDistToLocalMax <- apply(neg.cand.df, 1, function(i) min(abs(i-neg.tabu.ind)))
  
  # negative candidates should have at least a distance of minTargetBorderDist to the local clipped peak maxima that go beyond the quantile and the absolute threshold
  neg.cand.pos <- clippedBasePosF$pos[neg.cand.df$ind[localMax(neg.cand.df$minDistToLocalMax, threshold = minTargetBorderDist)]]
  
  
  # match 5' and 3' max positions based on peak height and distance. Peak height is more important??
  localMax.53 <- merge(localMax.5p, localMax.3p, by=NULL, suffixes = c(".5p", ".3p")) #cartesian product
  
  
  # score rates every pair of 3p-5p clipped bases peaks: the smaller the score the better the match
  localMax.53 %<>% 
    dplyr::mutate_(posDiff= ~pos.3p - pos.5p, ## how far (in nt) lies 5p peak ahead of 3p
                   # rank.posDiff a score that assesses the 5' and 3'-peak distance.
                   # 1st factor: Small differences get small rank, hence lower 1st factor. sqrt() for posDiff rank gives less weight to pos diff
                   # 2nd factor: give extra benefit to distances below 50.  pmin(1, pmax(0.5, sqrt(x/50))) lives betw 0.5 and 1. 
                   #+ Simpler step-function version is:  * ifelse( abs(localMax.53$posDiff) <= 50L, 0.5, 1) 
                   rank.posDiff= ~sqrt(rank(abs(posDiff))/n()) * pmin(1, pmax(0.5, sqrt(abs(posDiff) / 50L))),
                   score= ~rank.5p * rank.3p * rank.posDiff) %>% 
    dplyr::arrange_(~score) # best hits first
  
  
  # compare scores of regions that are adjacent with respect to score
  avgScoreDiff <- median(diff(localMax.53$score))
  
  # keep only combinations with best scores
  # drop entries as soon as a coordinate (pos.5p or pos.3p) has been used before, with a better (=lower) score
  localMax.53 %<>% 
    dplyr::filter_(~ ! duplicated(pos.5p), ~ ! duplicated(pos.3p))
  
  
  # mkuhn, 2015-03-03: if running on small region it can be NULL
  #stopifnot( NROW(localMax.53) >= 1L )
  
  
  # ZZZ assume at most one SV in target region
  # decide if we return two peak-combinations or only one peak-combination
  peakInfo <- NULL
  
  # return two peaks if they are closer together than average (according to combined score)
  peakInfo <- if ( NROW(localMax.53) >= 2L ) #&& (localMax.53$score[2L] - localMax.53$score[1L]) <= avgScoreDiff )
    localMax.53[1:2,] else 
      if (NROW(localMax.53) >= 1L) localMax.53[1L,]
  # else  if (NROW(localMax.53) == 1L) # used to give NULL in case more than 1 entry that is not close together
  
  # 
  # Plotting
  if (isTRUE(do.plot)){
    plot( relSbases_5p ~ pos, data=clippedBasePosF, type="l", col="darkblue", ylab="prop. clipped bases",
          ylim=c(0, 1.05*max(c(clippedBasePosF$relSbases_5p, clippedBasePosF$relSbases_3p))) )
    lines(relSbases_3p ~ pos, data=clippedBasePosF, type="l", lty=2, lwd=1.5, col="green")
    rug(x=localMax.5p$pos, side=3, col="darkblue", lwd=2)
    rug(x=localMax.3p$pos, side=3, col="green", ticksize = -.02, lwd=2)
  }
  
  
  return( list(clippedBasePos=clippedBasePos, peakInfo=peakInfo, nonPeakPos=neg.cand.pos))
}




#' Find start and end coordinates for feature region within target/covered region
#' 
#' @param pInfo peak info dataframe
#' @return feature region (for now, as vector?)
findFeatureRegion <- function(pInfo){
  if (is.null(pInfo)){
    # find startPos and endPos
    startPos <- floor(min(pInfo$pos.5p, pInfo$pos.3p)) - 5L
    endPos <- ceiling(max(pInfo$pos.5p, pInfo$pos.3p)) + 5L
    if (NROW(pInfo) == 1L){
      startPos <- startPos - 15L
      endPos <- endPos + 15L
    }
  }
  
  return(c(startPos, endPos))
}