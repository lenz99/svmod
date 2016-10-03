

# functions via Rcpp
# f <- Rcpp::cppFunction('
# int fibonacci(const int x) {
# if (x < 2)
# return x;
# else
# return (fibonacci(x - 1)) + fibonacci(x - 2);
# }
# ')
#
# Expressions via Rcpp:
# Rcpp::evalCpp('std::numeric_limits<double>::max()')


#' Create the patient ID string in canonical form
#' @param patIdNbr the patient ID as numeric
#' @export
getPatStr <- function(patIdNbr){
  paste0("pat",formatC(patIdNbr, width=5, flag="0"))
}

#' Get the patient ID as number from patient ID string
#' @param patIdStr the patient ID as string
#' @export
getPatNbr <- function(patIdStr){
  as.numeric(substr(patIdStr, 4, 8))
}



#' Checks given files if they are valid.
#' 
#' Do the files exist and are they non-empty?
#' @param f Filename(s)
#' @export
checkFile <- Vectorize(function(f){
  if (is.null(f)) return(FALSE)
  fsize <- file.info(f)$size
  
  ! is.na(fsize) && fsize > 512L
})


#' Count number of lines in a file
#' 
#' It uses the command line tool wc.
#' @param txtFile the filename of the text file
getNumberOfLinesUsingWC <- function(txtFile){
  if (! checkFile(txtFile)){
    logwarn("File %s does not exist!", txtFile)
    return(NULL)
  }
  
  wcOutput <- system(paste0("wc -l ", txtFile), intern=TRUE)
  wcOutput.regex <- regexec("[[:digit:]]+", wcOutput)[[1]]
  
  strtoi(substr(wcOutput, start=wcOutput.regex, stop=wcOutput.regex+attr(wcOutput.regex, "match.length")-1), base=10L)
}


# Find SV-intervals within candidate regions.
#
# look into the mapping of sufficiently covered regions 
# and find intervals around SV-like points within candidate region. The candidate regions are, for instance, based on clipping proportion.
# mapping: reads that map on a covered region
# works only for single chromosome
# @return returns matrix with interval number, start and stop position (or NULL)
# @note Only used in findIntervals-function below.
###
findIntervalsInCandidateRegions <- function(mapping, chrom=CHROM){
  stopifnot( is.data.frame(mapping) && NROW(mapping) >= 3L )
  loginfo("Candidate region with %d reads.", NROW(mapping))
  
  mapping$clipped  <- factor(ifelse(grepl("S", mapping$cigar), yes=2, no=1), levels=1:2, labels=c("noS", "S"))
  
  
  #moving averages from package caTools. 
  #+zoo::rollapply alternative to caTools
  
  #plot(rollapply(as.numeric(mapping$clipped), 30, FUN=mean, na.rm=TRUE, fill=NA) ~ pos, data=mapping, main="clipped")
  #plot(rollapply(mapping$mapq, 30, FUN=mad, na.rm=T, fill=NA) ~ pos, data=mapping, type="l", main="mapq")
  
  # Take product of windowed mean clipped proportion and std. deviation of mappingQ 
  #SVscore <- rollapply((rollapply(as.numeric(mapping$clipped), 30, FUN=mean, na.rm=TRUE, fill=NA)+1) * (rollapply(mapping$mapq, 30, FUN=mad, na.rm=T, fill=NA)+1), width=80, FUN=mean, fill=NA)
  #ZZZ window-width k=: should it depend on candidate region length
  #SVscore <- runmean((runmean(as.numeric(mapping$clipped), k=51, endrule="NA")+10) * (log(runsd(mapping$mapq, k=51, endrule="NA")+1)), k=21, endrule="mean"); SVscoreType <- "clipped * mapQ"  # clipped * mapQ
  #SVscore <- runsd(mapping$mapq, k=51, endrule="NA"); SVscoreType <- "mapQ-sd"  ## only mapQ-sd
  
  #ZZZ neighborhood is defined about the ordering of reads in the mapping. It does not consider how many base-pairs the reads really lie apart
  ## For each read give the proportion of clipped reads in its neighborhood
  #require(caTools)
  SVscore <- runmean(as.numeric(mapping$clipped)-1, k=51, endrule="NA"); SVscoreType <- "clipped-prop" 
  SVscore.cutoff <- quantile(SVscore, probs=0.85, na.rm=TRUE)  # what is a "high proportion of clipped reads in a neighbourhood": only 15% of reads have that enriched neighborhood
  # look for reads where 80% of neighboring reads have "enriched neighborhood"
  SVscore.high <- caTools::runquantile(as.numeric(SVscore > SVscore.cutoff), k=51, probs=0.2, endrule="NA") 
  #SVscore.high <- rollapply(as.numeric(SVscore > SVscore.cutoff), width=11, FUN=median, na.rm=TRUE, fill=NA) ##zoo
  
  # indicating intervals of SV-candidates within covered region
  SVscore.high.ind <- which(SVscore.high == 1)
  nbr.HighSVReads <- length(SVscore.high.ind)
  loginfo("Number of reads that high SV-score: %d", nbr.HighSVReads)
  if ( nbr.HighSVReads == 0){
    loginfo("Skipping candidate region because it does not have areas with high SV reads density.")
    return(invisible(NULL))
  }
  pos.high <- mapping$pos[SVscore.high.ind]
  mapping$svreg <- -1
  is.na(mapping$svreg) <- is.na(mapping$pos)
  mapping$svreg[SVscore.high.ind] <- cut(pos.high, breaks=unique(pos.high[c(1, which(diff(SVscore.high.ind) > 1L), length(pos.high))]), include.lowest=TRUE, labels=FALSE)
  loginfo("Found ::%s:: svreg levels", unique(mapping$svreg))
  
  plot( SVscore ~ pos, data=mapping, subset=!is.na(SVscore), pch=".", main=SVscoreType)
  abline(h=SVscore.cutoff, lty=2, col="grey")
  lines(max(SVscore, na.rm=T)*SVscore.high ~ pos, data=mapping, col="red", type="l", lty=3)
  
  myMapping  <<- mapping #ZZZ debug
  
  # give coordinates of intervals
  mapping_posRange <- range(mapping$pos, na.rm=TRUE)
  
  mapping_help <- mapping[which(mapping$svreg > 0), c("rname", "pos", "svreg")]
  retMat <- as.matrix(aggregate(mapping_help$pos, by=list(svreg=mapping_help$svreg), FUN=range, na.rm=TRUE))
  colnames(retMat) <- c("svreg", "lower", "upper")
  # how many BPs around the interval borders
  EXTRA_BP <- 150L
  retMat[,2] <- pmax(retMat[,2] - EXTRA_BP, mapping_posRange[1])
  retMat[,3] <- pmin(retMat[,3] + EXTRA_BP, mapping_posRange[2])
  
  nbrIntervalsInRegion0 <- NROW(retMat)
  
  # merging intervals
  retMat <- as.matrix(ranges(reduce(GenomicRanges::GRanges(seqnames=rep(chrom, NROW(retMat)),
                                                           #paste(retMat[,1], chrom, sep="_"),
                                                           ranges=IRanges(retMat[,2], retMat[,3])))))
  retMat[,2] <- retMat[,1] + retMat[,2] - 1
  retMat <- cbind(1:NROW(retMat), retMat)
  
  loginfo("Found %d intervals (%d intervals before merging) which are SV-like.", NROW(retMat), nbrIntervalsInRegion0)
  return(retMat)
}


#
# identifies intervals that are potential regions with SVs in CHROM
# first it looks for regions that are sufficiently covered
# Within  each covered region it delineates intervals that are SV-candidate (via clipped reads)
#+ ZZZ: we are actually hunting for differential SVs (so better look into reads from CTR sample as well)
####
findIntervals <- function(tum.bam){ #, ctr.bam){
  tum.bamFileInfo <- parseSimFile(tum.bam)
  #ctr.bamFileInfo <- parseSimFile(ctr.bam)
  
  # use bedtools to get coverage
  stopifnot( exists("BEDTOOLS") )
  
  tmpFile.tum <- tempfile(pattern="tum_", fileext=".cov")
  #tmpFile.ctr <- tempfile(pattern="ctr_", fileext=".cov")
  
  
  # -dz: only coverage >0, zero-based positions
  system(paste0(BEDTOOLS, " genomecov -dz -ibam ", tum.bam," -g  > ", tmpFile.tum),
         intern=FALSE, wait=TRUE)
  #system(paste0(BEDTOOLS, " genomecov -dz -ibam ", ctr.bam," -g  > ", tmpFile.ctr),
  #       intern=FALSE, wait=TRUE)
  
  # delete temporary files when function exits.
  on.exit(expr=file.remove(tmpFile.tum)) #, tmpFile.ctr))
  
  cov.tum <- read.delim(file=tmpFile.tum, header=FALSE)
  #cov.ctr <- read.delim(file=tmpFile.ctr, header=FALSE)
  
  names(cov.tum) <- c("chrom", "pos", "cov_tum")
  #names(cov.ctr) <- c("chrom", "pos", "cov_ctr")
  
  ### ZZZ single chromosome!!
  cov.tum <- cov.tum[cov.tum$chrom == CHROM,]
  #cov.ctr <- cov.ctr[cov.ctr$chrom == CHROM,]
  
  # look for gaps in coverage. Gaps should have minimum length [bp] (gap might also come from insertion-SV)
  MIN_GAP_LENGTH <- 500L
  
  # define intervals where there is read coverage and that are separated by a gap
  #cov.ctr$interv_ctr <- cut(cov.ctr$pos, breaks=cov.ctr$pos[c(1, which(diff(cov.ctr$pos)>MIN_GAP_LENGTH), NROW(cov.ctr))], include.lowest=TRUE, labels=FALSE)
  cov.tum$interv_tum <- cut(cov.tum$pos, breaks=cov.tum$pos[c(1, which(diff(cov.tum$pos)>MIN_GAP_LENGTH), NROW(cov.tum))], include.lowest=TRUE, labels=FALSE)
  
  
  
  # filter on those intervals that have enough coverage in TUMOR sample
  MIN_QUANT75_COVERAGE <- 10L # at least 30% of bases have coverage above 10
  minCov.tum <- aggregate(x=cov.tum[, "cov_tum"], by=list(interv_tum=cov.tum$interv_tum), FUN=function(x){quantile(x, probs=0.7)})
  names(minCov.tum)[NCOL(minCov.tum)] <- "cov_q75"
  minCov.tum <- minCov.tum[minCov.tum$cov_q75 >= MIN_QUANT75_COVERAGE,]
  loginfo("We have %d covered regions which fulfill minimum coverage requirements.", NROW(minCov.tum))
  
  #   cov.interv.ctr <- aggregate(x=cov.ctr[, "cov_ctr"],
  #                              by=list(interv_ctr=cov.ctr$interv_ctr),
  #                              FUN=function(x){quantile(x, probs=0.75)})
  #   names(cov.interv.ctr)[NCOL(cov.interv.ctr)] <- "cov_q75"
  #   cov.interv.ctr <- cov.interv.ctr[cov.interv.ctr$cov_q75 >= MIN_QUANT75_COVERAGE,]
  #   
  # filter on those intervals where minimum coverage is reached. 
  cov.tum2 <- merge(cov.tum, minCov.tum)
  cov.tum2 <- cov.tum2[order(cov.tum2$chrom, cov.tum2$pos),]
  #   cov.ctr2 <- merge(cov.ctr, cov.interv.ctr)
  
  #   cov.interv.ctr <- merge(cov.interv.ctr, aggregate(x=cov.ctr2[, "pos"], simplify=TRUE,
  #                                                     by=list(interv_ctr=cov.ctr2$interv_ctr),
  #                                                     FUN=min))
  #   names(cov.interv.ctr)[NCOL(cov.interv.ctr)] <- "interv_l"
  #   cov.interv.ctr <- merge(cov.interv.ctr, aggregate(x=cov.ctr2[, "pos"], simplify=TRUE, 
  #                                                     by=list(interv_ctr=cov.ctr2$interv_ctr),
  #                                                     FUN=max))
  #   names(cov.interv.ctr)[NCOL(cov.interv.ctr)] <- "interv_h"
  
  
  minCov.tum <- merge(minCov.tum, aggregate(x=cov.tum2[, "pos"], simplify=TRUE,
                                            by=list(interv_tum=cov.tum2$interv_tum),
                                            FUN=min))
  names(minCov.tum)[NCOL(minCov.tum)] <- "interv_l"
  minCov.tum <- merge(minCov.tum, aggregate(x=cov.tum2[, "pos"], simplify=TRUE, 
                                            by=list(interv_tum=cov.tum2$interv_tum),
                                            FUN=max))
  names(minCov.tum)[NCOL(minCov.tum)] <- "interv_h"
  
  
  # Assess covered intervals:
  # Look only into sufficiently covered intervals (at least in tumor sample)
  
  # further narrow the list of candidate intervals 
  # and also narrow the interval borders so that the borders have a meaning/are related to SV
  # check mapping for differential pattern (clipped reads and quality) with scanBam
  
  loginfo("Start looking into the mapping of those candidate covered regions.")
  
  
  # what fields to read out from BAM file
  what <- setdiff(scanBamWhat(), c("seq", "qual", "qwidth"))
  ##, isFirstMateRead=TRUE
  
  #ZZZ RangesList() needs to be properly used for multiple chromosomes! chrom as names
  ##here chr5: hard-coded because of NPM1
  #   myParam_ctr <- ScanBamParam(which=RangesList(chr5=IRanges(start=cov.interv.ctr$interv_l, end=cov.interv.ctr$interv_h)), flag=scanBamFlag(isPaired=TRUE, isProperPair=NA), what=what)
  myParam_tum <- ScanBamParam(which=RangesList(chr5=IRanges(start=minCov.tum$interv_l, end=minCov.tum$interv_h)), flag=scanBamFlag(isPaired=TRUE, isProperPair=NA), what=what)
  
  # have list of dataframes with the mapping per region which reaches minimum coverage in BAM file
  mappings.tum <<- lapply(scanBam(file=tum.bam, param=myParam_tum), as.data.frame)
  #mappings.ctr <- lapply(scanBam(file=ctr.bam, param=myParam_ctr), as.data.frame)
  
  #   # add name of covered region to data.frame ##seems not used/necessary
  #   for (i in 1:length(mappings.tum)){
  #     mappings.tum[[i]]$region <- names(mappings.tum)[i]
  #   }#rof
  
  
  #   #ZZZ Ideally, I would relate intervals betw tumor and control: find intervals in tumor and in control that cover (more or less) the same region. 
  #   #+Could be that I need to merge two intervals in -say- tumor probe (when a big deletion occured)
  #   #+For now, I look only at the tumor intervals and check clipped reads and quality
  #   
  #   # distribution of percentages of clipped reads with mapQ > 0 in covered regions
  #   clippedProp.interv.tum <- sapply(mappings.tum, function(mapping){sum(grepl("S", mapping$cigar) & mapping$mapq > 0) / NROW(mapping)})
  #   
  #   PERC_INTERV <- ceiling(length(mappings.tum)/5)
  #   # select those tumor intervals with highest  proportion of clipped reads
  #   highClippedRegions.tum <- sort(clippedProp.interv.tum, decreasing=TRUE)[1:PERC_INTERV]
  #   loginfo("The covered regions with highest clipped-reads proportion are %s.", names(highClippedRegions.tum))
  #   
  #   # robust variance
  #   qualVar.interv.tum <- sapply(mappings.tum, function(mapping){mad(mapping$mapq, na.rm=TRUE)})
  #   highMapQVarRegions.tum <- sort(qualVar.interv.tum, decreasing=TRUE)[1:PERC_INTERV]
  #   loginfo("The covered regions with highest mapQ-variability are %s.", names(highMapQVarRegions.tum))
  #   
  #   #covered regions that have relatively many clipped reads and relatively high qual variability
  #   candidateRegions1 <- intersect(names(highClippedRegions.tum), names(highMapQVarRegions.tum))
  #   
  #   mappings2.tum <- mappings.tum[candidateRegions1]
  #   loginfo("The following candidate covered regions remain %s.", candidateRegions1)
  
  # return list of intervals that look SV-like
  #return(lapply(mappings2.tum, findIntervalsInCandidateRegions))
  retList <- lapply(mappings.tum, findIntervalsInCandidateRegions)
  
  retMat <- NULL
  for (i in 1:length(retList)){
    retMat <- rbind(retMat, retList[[i]])
    #ZZZ rownames important?: bamFile + region?
  }
  
  #return(t(sapply(retList[! sapply(retList, is.null)], rbind)))
  return(retMat)
}



# wrapper to findIntervals() in different mapping BAM-files
# returns a matrix with coordinates per candidate interval
###
findIntervalsWrapper <- function(bamFiles){
  stopifnot( ! is.character(bamFiles) || length(bamFiles) > 0 )
  myResult <- NULL 
  samples <- character(0)
  
  for (i in bamFiles){
    loginfo("Find intervals for BAM file %s.", i)
    myCurrentInterval <- findIntervals(i)
    if ( NCOL(myCurrentInterval) > 0 ){
      samples <- c(samples, rep(i, NROW(myCurrentInterval)))
      myResult <- rbind(myResult, myCurrentInterval)
    }#fi
  }
  return(list(sample=samples, interval=myResult))
}


#' Replaces random positions of a given DNA sequence with a random nucleotide.
#' 
#' This function can be useful when simulating a virtual patient that has personal (germline) SNVs.
#' @author mkuhn, 201402
#' @param DNAseq the DNA sequence as \code{DNAString} object. Works also for rectangular \code{DNAStringSet} objects.
#' @param nbrSNVs number of SNVs (but could be silent, i.e. replace A with A)
#' @note ZZZ It was better when we set here a mutation rate, not the number of SNVs..
injectRandomSNVs <- function(DNAseq, nbrSNVs=1L){
  seqLen <- if ( inherits(DNAseq, "DNAString") ) length(DNAseq) else Biostrings::width(DNAseq)
  
  # nothing to do when empty sequence
  if ( all(seqLen == 0L) ) return(DNAseq)
  
  if ( is.numeric(nbrSNVs) && nbrSNVs > 0L ){
    nbrSNVs <- round(nbrSNVs)
    
    if ( inherits(DNAseq, "DNAString") && length(seqLen) == 1L ){
      myAt <- sample(seqLen, nbrSNVs) 
      myLetter <- Biostrings::DNAString(paste(sample(Biostrings::DNA_BASES, nbrSNVs, replace=TRUE), collapse=""))
    } else {
      stopifnot( inherits(DNAseq, "DNAStringSet") && length(table(Biostrings::width(DNAseq))) == 1L )
      
      seqWidth <- Biostrings::width(DNAseq)[1]
      h <- matrix(data = FALSE, nrow = length(DNAseq), ncol = seqWidth)
      for (i in 1:NROW(h)) h[i, sample(seqWidth, nbrSNVs, replace=FALSE)] <- rep(TRUE, nbrSNVs)
      myAt <- h
      
      
      myLetter <- replicate(n=length(DNAseq), paste(sample(Biostrings::DNA_BASES, nbrSNVs, replace=TRUE), collapse=""))
    }
    
    DNAseq <- Biostrings::replaceLetterAt(DNAseq, at=myAt, letter=myLetter )
  }
  
  return(DNAseq)
}




#' Function to abbreviate a long sequence as a string
#' @param seqObj a <code>XString</code> object
#' @param k is number of bases that can remain on either side
abbrevSeq <- function(seqObj, k=11L){
  if (is.null(seqObj)) return(NULL)
  
  seqObj.len <- length(seqObj)
  #logdebug("Got seq-object ::%s:: of length %d.", seqObj, seqObj.len)
  if ( seqObj.len <= 2*k + 3) return(seqObj)
  
  subseq(seqObj, start=k+1, end=seqObj.len-k) <- Biostrings::DNAString("...")
  return(seqObj)
}


#' Find local maxima on a function given on a discrete support (i.e. discrete set of function values).
#' 
#' @param x numeric vector of given function values.
#' @param threshold value the local maximum is required to exceed. Defaults to \code{NULL}.
#' @return a vector of indices on x of the positions of local maxima. NULL if no valid x
#' @references http://stackoverflow.com/questions/6836409/finding-local-maxima-and-minima
#' @export
localMax <- function(x, threshold=NULL){
  if (is.null(x) || length(x) == 0L) return(NULL)
  if (! is.numeric(x)) return(NULL)
  if (length(x) == 1L) return(1L)
  
  y <- c(TRUE, diff(x) > 0L)
  y <- cumsum(rle(y)$lengths)
  y <- y[seq.int(1L, length(y), 2L)] # 1, 3, 5, ... starting from first, every second y-value: end of TRUE blocks
  if (x[[1L]] == x[[2L]])  y <- y[-1L] # if first two values are equal we get initially for y: TRUE, FALSE, .. but x[1] is no max
  
  if (! is.null(threshold)) y <- intersect(y, which(x >= threshold))
  return(y)
}



#' Helper function to facilitate the sub-sampling of training and test data based on the status variable.
#' 
#' It collects the requested number of training data with a given proportion of outliers. 
#' The remaining observations are used for test data. The indices for a sample is given with the requested proportion.
#' \code{NA} can occur and will be preserved. The NA-treatment was motivated by test-data where initially the status was missing. In simulation it is known.
#' 
#' Bootstrap samples are not supported right now.
#' 
#' @param dat.status logical vector IN=\code{TRUE}, OUT=\code{FALSE}
#' @param train.n.out size of training data. If \code{NULL} or invalid value the maximal size for training data is assumed that guarantees the outlier proportion.
#' @param prop.outl requested proportion of outliers in the data
#' @return list of vector of indices, first for training set, then for test set with same proportion, often no cases are left for test set. NA-indices are separate category. In case of failure returns NULL.
#' @export
subsampleInd <- function(dat.status, train.n.out=NULL, prop.outl=.15){
  ##type=c("subsampling", "bootstrap"),  type <- match.arg(type)
  
  stopifnot( is.logical(dat.status) || (is.factor(dat.status) && identical(levels(dat.status), STATUS_LEVELS)) )
  stopifnot( prop.outl > 0L, prop.outl < 1L )
  
  # mkuhn, 20150425: transform status factor into a logical vector
  if (is.factor(dat.status)) dat.status <- (dat.status == STATUS_LEVELS[1])
  
  # mkuhn, 2015-12-28: added special NA-treatment
  na.ind <- which(is.na(dat.status)) 
  
  # calculation
  n.in <- sum(! is.na(dat.status)) # formely (w/o special NA treatment) length(dat.status)
  n.normal <- sum(dat.status, na.rm = TRUE) # nbr of TRUE cases
  n.outl <- n.in - n.normal
  
  
  # maximal training set
  if ( is.null(train.n.out) || train.n.out <= 0L || train.n.out >= n.in ){
    train.prop.outl <- n.outl / n.in
    #loginfo("Give maximum size that guarantees the outlier proportion of %f.", prop.outl)
    if (train.prop.outl >= prop.outl){
      # too many outliers in training data
      train.norm.ind <- sample(which(dat.status)) # take all normals available and permute
      n.train.ref <- ceiling(length(train.norm.ind) / (1L-prop.outl))
      train.outl.ind <- sample(which(!dat.status), size=min(n.outl, n.train.ref - length(train.norm.ind)))
    } else {
      # too many normals
      train.outl.ind <- sample(which(!dat.status)) # take all outliers available
      n.train.ref <- floor(length(train.outl.ind) / prop.outl)
      train.norm.ind <- sample(which( dat.status), size=min(n.normal, n.train.ref - length(train.outl.ind)))
    }
  } else {
    ## only a subset of the training set is requested
    
    if (train.n.out > n.in){
      logging::logwarn("Requested more output then there is input!")
      return(invisible(NULL))
    }
    
    # requested numbers
    n.train.normal <- ceiling(train.n.out*(1L-prop.outl))
    n.train.outl <- train.n.out - n.train.normal
    
    if ( n.train.normal > length(which(dat.status)) ){
      #logwarn("Would need %d normal cases but there are only %d normal cases  in data! Sorry.", n.train.normal, n.normal)
      return(invisible(NULL))
    }
    
    if (n.train.outl > length(which(!dat.status)) ){
      #logwarn("Would need %d oultier cases but not enough outlier cases in data! Sorry.", n.train.outl)
      return(invisible(NULL))
    }
    stopifnot( n.train.outl >= 1L )
    
    train.norm.ind <- sample(which(dat.status), size = n.train.normal) 
    train.outl.ind <- sample(which(! dat.status), size = n.train.outl)
  }#esle
  
  train.ind <- c(train.norm.ind, train.outl.ind)
  
  
  # remaining indices (with TRUE or FALSE at dat.status)
  remain.norm.ind <- setdiff(which( dat.status), train.norm.ind)
  remain.outl.ind <- setdiff(which(!dat.status), train.outl.ind)
  
  # enforce outlier proportion in remaining data (that can be used as test data)
  odds.outl <- prop.outl / (1L-prop.outl)
  if ( isTRUE(length(remain.norm.ind) * odds.outl  > length(remain.outl.ind)) ){
    ## nbr of outliers is limiting
    #loginfo("Test data: Outlier data limiting.")
    n.test <- floor(length(remain.outl.ind) / prop.outl)
    test.norm.ind <- sample(remain.norm.ind, n.test - length(remain.outl.ind))
    test.outl.ind <- remain.outl.ind
  } else {
    ## normal limiting
    #loginfo("Test data: Normal cases limiting.")
    n.test <- ceiling(length(remain.norm.ind) / (1L-prop.outl))
    test.norm.ind <- remain.norm.ind
    test.outl.ind <- sample(remain.outl.ind, n.test - length(remain.norm.ind))
  }
  test.ind <- c(test.norm.ind, test.outl.ind)
  
  return( list( train.ind=train.ind, test.ind=test.ind, na.ind=na.ind) )
}


#' Wrapper for read simulation utility.
#' 
#' Specify either \code{fCov} or \code{rCount} but not both.
#' @param fastaFile character. Fasta file with reference sequence to simulate reads from.
#' @param ngs.prop contains the sequencing settings to be used for the simulation
#' @param rCount total number of reads/read pairs to be generated
#' @param fCov the fold of read coverage to be simulated
#' @param paired logical flag if paired sequencing is requested.
#' @return logical flag for success status
readSimulatorWrapper <- function(fastaFile, prefixRead="SIM", prefixFasta, paired=TRUE, ngs.prop=NULL, fCov=NULL, rCount=NULL){
  
  stopifnot( ! (is.null(fCov) && is.null(rCount)), is.null(fCov) || is.null(rCount) ) # exactly one of both is not null!
  
  artOptions <- paste("--noALN", if (paired) "--paired", "--seqSys ", SEQ_SYSTEM, 
                      if (! is.null(fCov)) paste("--fcov", fCov), if (! is.null(rCount)) paste("--rcount", rCount))
  
  simCommand <- sprintf("%s  %s  --id %s --in %s --out %s --len %d  --mflen %d --sdev %d", 
          ART_EXE, artOptions, prefixRead, fastaFile, file.path(dirname(fastaFile), prefixFasta),
          ngs.prop$read.length, ngs.prop$pe.ins.mean, ngs.prop$pe.ins.sd)
  
  simCall <- system(simCommand, intern=FALSE)
  if ( simCall == 0L ) logging::loginfo("Reads simulated") else {
    logging::logerror("Read Simulation FAILED for with call %s", simCommand)
    return(invisible(FALSE))
  }
  
  invisible(TRUE)
}


#' BWA-mapper function wrapped in R.
#' @param bwa_opts character string that describes a BWA mapping profile. Currently, "clippy" with less penalties for clipping. Standard means all penalties have their default values.
#' @param bwa_optsExtra character string with extra options (except for std mapping profile, threads, smart pairing). For instance, -I for insert size distribution.
#' @param nbrThreads number of threads to use for mapping
#' @param interleaved logical flag if smart pairing (interleaved paired-end reads) is to be used. If \code{TRUE} only the first input file \code{inputFile1} is used that should contain both reads for pairs.
#' @return The BAM filename upon mapping, sorting and indexing success. Otherwise FALSE invisibly.
#' @note The BAM-flag "proper pair" is only set by BWA if enough reads are in the FASTQ files. Is pair rescue done?
#' @export
mapWrapper <- function(inputFile1, inputFile2=NULL, interleaved=FALSE, outputDir=dirname(inputFile1), bamFilename=NULL,
                       genomeRef=REF_GEN, bwa_opts=c("clippy", "standard"), bwa_optsExtra=NULL, nbrThreads=2L){
  bwa_opts <- match.arg(bwa_opts)
  
  # handle BWA options
  bwa_opts <- switch(bwa_opts,
                     clippy="-a -Y -L3 -U15", # -M", we do not use -M as we want to have supplementary alginments for clipped read ends!?
                     standard=,
                     "")
  
  #-C: adds FASTA comment to SAM as extra field [when it complies like BC:Z:AGCGA]
  #-a: all secondary alignments as separate entries in BAM file 
  #-Y: use soft clipping for supplementary alignments
  #-L: penalty for clipped reads  [default 5]
  #-U: penalty for unpaired read-pair [default 17]
  #-M: mark shorter split hits as secondary (and not as supplementary)
  #-t: number of threads
  BWA_OPT <- if (is.numeric(nbrThreads) && nbrThreads > 0L) paste0(bwa_opts, " -t", ceiling(nbrThreads)) else bwa_opts 
  
  if (! is.null(bwa_optsExtra) && nzchar(bwa_optsExtra)){
    BWA_OPT <- paste(BWA_OPT, bwa_optsExtra)
  }
  
  
  if ( is.null(inputFile1) || ! is.character(inputFile1) ){
    logerror("No valid input files given for 1st input FASTQ-file %s.", paste(inputFile1, collapse="::"))
    return(invisible(FALSE))
  }
  
  if ( length(inputFile1) > 1L && is.null(inputFile2) ){
    logging::loginfo("Map wrapper: taking two FASTQ files from argument inputFile1.")
    inputFile2 <- inputFile1[2L]
    inputFile1 <- inputFile1[1L]
  }
  
  
  
  if (isTRUE(interleaved)){
    logging::loginfo("Interleaved-mapping mode")
    if ( ! is.null(inputFile2) ) logging::logwarn("2nd input file is ignored!")
    inputFile2 <- NULL
    BWA_OPT <- paste(BWA_OPT, "-p")
    
  } else if (is.null(inputFile2) || ! is.character(inputFile2)){
    logging::loginfo("SE-mapping mode")
  } else {
    logging::loginfo("PE-mapping mode")
  }
  
  
  if (! checkFile(inputFile1) || (! is.null(inputFile2) && ! checkFile(inputFile2)) ){
    logerror("Input files are not valid. file1: %s and file2: %s.", inputFile1, inputFile2)
    return(invisible(FALSE))
  }
  
  
  
  
  # extract BAM filename
  if (is.null(bamFilename)){
    bamFilename <- inputFile1 %>%
      basename %>%
      stringr::str_replace(pattern = "[.](fq|subfq|fastq)$", replacement = "") %>% 
      stringr::str_replace(pattern = "_(R)*1$", replacement = "") %>% 
      paste0(".bam")
  }
  
  
  # BAM output file
  outputFile0 <- file.path(outputDir, bamFilename)
  
  bwaRet <- system(paste(BWA_EXE,"mem", BWA_OPT, genomeRef, inputFile1, inputFile2,
                         "|", SAMTOOLS_EXE,"view -@1 -Sb - > ", outputFile0), 
                   intern=FALSE)
  if ( bwaRet != 0 || ! checkFile(outputFile0) ){
    logerror("BWA mem failed for input %s and %s.. Sorry.", inputFile1, inputFile2)
    return(invisible(FALSE))
  } else logdebug("Mapping of %s and %s done.", inputFile1, inputFile2)
  
  
  ## SORT BAM file
  # sortRet <- system(paste(SAMTOOLS_EXE, "sort", outputFile, sub(".bam$", "", outputFile)), intern=FALSE)
  # if ( sortRet != 0 && checkFile(outputFile) ){
  #   logerror("BAM SORT failed on mapping %s.", outputFile)
  #   return(invisible(FALSE))
  # } else logdebug("BAM SORT done on mapping %s.", outputFile)
  on.exit(try(unlink(outputFile0)))
  outputFile1 <- Rsamtools::sortBam(outputFile0, destination = sub(".bam$", ".s", outputFile0))
  if (! checkFile(outputFile1)){
    logging::logerror("Failed to sort BAM file %s.", outputFile0)
    return(invisible(FALSE))
  }
  logging::loginfo("Sorted the BAM file and saved as %s.", outputFile1)
  
  
  
  ## INDEX BAM files
  Rsamtools::indexBam(outputFile1)
  logging::loginfo("Indexed the BAM file %s.", outputFile1)
  # # NORMAL probe mapping: indexing
  # indexRet <- system(paste(SAMTOOLS_EXE, "index", outputFile), intern=FALSE)
  # 
  # if (indexRet != 0 ){
  #   logerror("BAM INDEX failed on %s", outputFile)
  #   try(file.remove(outputFile), silent = TRUE)
  # } else logdebug("BAM indexing done on %s.", outputFile)
  
  return(outputFile1)
}



#' Hamming window
#' @param n length of window
#' @param normalized. boolean flag if weights should be normalized that they sum to n.
getHammingWindowWeights <- function(n, normalized=TRUE){
  weights <- 0.54 + 0.46 * cos( 2 * pi * seq(from=-n/2, length.out = n) / n)
  #weights <- 0.54 - 0.46 * cos( 2 * pi * seq(from=0, length.out = n) / n)
  
  if (normalized) weights <- weights * length(weights) / sum(weights)
  
  return(weights)
}



#' Graphical exploration of published SV-data from Database of Genome Variation.
#' @param dgvFile variants from the genome variant database, as an RDS object
exploreGenomeVariants <- function(dgvFile=file.path(BASEDIR, "refGen", "genomVar_hg19_20150723.rds")){
  if ( ! file.exists(dgvFile) ){
    logwarn("File %s not found. Sorry.", dgvFile)
    return(invisible(NULL))
  }
  
  #dgv <- read.delim(dgvFile)
  dgv <- readRDS(file=file.path(BASEDIR, "refGen", "genomVar_hg19_20150723.rds"))
  dgv %<>% dplyr::mutate_(variantlength=~end - start + 1L,
                          hasGene=~factor(genes != '', labels = c("no", "yes")))
  
  #ggplot(dgv, aes(x=variantsubtype, y=variantlength, fill=hasGene)) + geom_boxplot() + theme(axis.text.x=element_text(angle=45L, hjust=1L))
  ggplot(dgv, aes(x=variantsubtype, y=variantlength, fill=hasGene)) + geom_boxplot() + ylim(c(0,500)) + theme(axis.text.x=element_text(angle=45L, hjust=1L))
  ggplot(dgv, aes(x=variantsubtype, y=variantlength)) + geom_boxplot(varwidth = T) + ylim(c(0,500)) + theme(axis.text.x=element_text(angle=45L, hjust=1L))
  ggplot(dgv, aes(x=variantsubtype, y=variantlength, fill=hasGene)) + geom_boxplot() + ylim(c(0,200)) + theme(axis.text.x=element_text(angle=45L, hjust=1L))
}
