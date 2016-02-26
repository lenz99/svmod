# stuff for SVs and SV injection
#
# mkuhn, 2016-01-20
# There is another SV-injection function for virtual patients because it also uses SNV injection at the same time: there is a test routine to compare the two functions.
#+Ideally, there would be one function!


#' SV in a genomic context relative to a reference genome (S3 constructor).
#' 
#' The SV sequence for insertions is randomly drawn.
#' Strand is ignored. Translocations are not supported right now because we focus on local SVs or would show up as two independent SVs.
#' Field \code{posL} points to the left-most genomic position where the injected SV will produce a clipped-base peak. For duplication, this will be the start of the duplicated region.
#'
#' @param posL left-most position of SV in 1-based coordinates (like IRange). !!fixme!! For duplication it is the first (left-most) duplicated base (right after the original DNA-template to be duplicated)
#' @param seqLen numeric length of SV to be generated
#' @note Translocations are not directly supported. 
#' @export
SV <- function(type=SV_TYPES, chrom, posL, seqLen=SIM_SV_LENGTH_MEAN, refGen=REF_GEN){
  type <- match.arg(type)
  stopifnot( is.character(refGen), checkFile(refGen) )
  #stopifnot( ! is.null(seq) && ! isTRUE(generateSeq), is.null(seq) && isTRUE(generateSeq) )
  #stopifnot( ! isTRUE(generateSeq) || generateSeqLen >= 1L )
  
  stopifnot( seqLen >= 1L )
  seqLen <- round(seqLen)
  
  # sequence information from reference genome
  refGenHandle <- open(Rsamtools::FaFile(refGen))
  on.exit(try(close(refGenHandle), silent = TRUE))
  refGenInfo <- Rsamtools::seqinfo(refGenHandle)
  
  # checks on sequences
  stopifnot( chrom %in% seqnames(refGenInfo) )
  stopifnot( posL < GenomeInfoDb::seqlengths(refGenInfo)[chrom] )
  
  
  # # cast character to DNAString if necessary
  # if (is.character(seq)) seq <- Biostrings::DNAString(seq)
  
  # SV locus
  locus <- GenomicRanges::GRanges(seqnames = chrom, #if (is.null(seq)) generateSeqLen else nchar(seq)) ) 
                                  ranges = IRanges::IRanges(start = posL, width = seqLen)) 
  
  
  # reference sequence after the SV BP (useful for deletion and inversion)
  refSeq <- Biostrings::getSeq(refGenHandle, param=locus)[[1L]]

  # sequence that constitutes the SV
  svSeq <- switch(type,
                  insertion=Biostrings::DNAString(paste(sample(Biostrings::DNA_BASES, seqLen, replace=TRUE), collapse="")),
                  duplication=, # the duplicated sequence
                  deletion=refSeq, # the deleted sequence
                  inversion=Biostrings::reverseComplement(refSeq)
  )
  
   
  
  stopifnot( ! is.null(svSeq) , inherits(svSeq, "DNAString"), ! is.null(locus) )
  
  z <- list(
    type=type,
    # seq contains modifying sequence that makes up the SV, i.e. for inversion the revCompl sequence, for deletion the deleted part
    seq=svSeq,  
    locus=locus,
    refGen=refGen,
    refGenInfo = refGenInfo
  )
  class(z) <- "SV"
  return(z)
}

getChrom <- function(sv) GenomeInfoDb::seqlevels(sv[["locus"]])[1L]
getPosL <- function(sv) GenomicRanges::start(sv[["locus"]])
getLength <- function(sv) GenomicRanges::width(sv[["locus"]])
getLocus <- function(sv) sv[["locus"]]
getSVSeq <- function(sv) sv[["seq"]]
getShortName <- function(sv) paste("injSV", sv[["type"]], getLocus(sv), sep="_")
toString.SV <- function(x) {
  svGR <- x[["locus"]]
  svSeq <- getSVSeq(x)
  seqLen <- nchar(svSeq)
  seqStart <- as.character(Biostrings::subseq(svSeq, start=1L, end=min(seqLen, 10L)))
  if (seqLen > 10L){
    seqStart <- paste0(seqStart, "..")
  }
  sprintf("%s-SV [%s] of length %d at %s:%d relative to %s",
          x[["type"]], seqStart,
          GenomicRanges::width(svGR), GenomeInfoDb::seqlevels(svGR), GenomicRanges::start(svGR), basename(x[["refGen"]]) )
  }
print.SV <- function(x) print(toString.SV(x))

#' Get sequence of genome with injected SV at given region in reference coordinates.
#' 
#' We assume that the genomic region either completely covers the relevant SV-locus or not at all.
#' Through the SV the returned sequence can be shorter or longer than (for unbalanced SVs)  or equal (inversion) to the requested region.
#' For duplication the genomic region must start at least SV-length before breakpoint to encompass full SV context (i.e. the duplicated block twice)
#' 
#' Assumption here: relevant SV-locus is completely covered by requested context region or not at all (and then the reference seq is returned)
#' @param gr genomic region in reference sequence as \code{GRanges} object from where to extract the sequence
#' @return Reference sequence including the injected SV as \code{DNAString} or the unchanged reference sequence if SV not hit
#' @note The function could also give SV-context when requested region only partially covers the SV region, but is more cumbersome
getContextSeq <- function(sv, gr) {
  stopifnot( is(gr, "GRanges"), length(gr) == 1L )
  stopifnot( getChrom(sv) == GenomeInfoDb::seqlevels(gr)[1L] )
  
  #for duplication this becomes here both copies of the duplicated block
  isGainSV <- isTRUE(sv$type == 'duplication' || sv$type == 'insertion' )
  isLossSV <- isTRUE(sv$type == 'deletion')
  isNeutralSV <- isTRUE(sv$type == 'inversion')
  stopifnot( isGainSV + isLossSV + isNeutralSV == 1L )
  
  # for gains, by convention I consider the left (5') BP: once this left BP is covered, the SV-gain is included in the SV-context sequence. Does this matter?
  # relevant SV-locus:
  relSVLocus <- if ( isGainSV ){
    #GenomicRanges::resize(relSVLocus, width = 2L * GenomicRanges::width(svLocus), fix = "end")
    GenomicRanges::resize(getLocus(sv), width = 1L, fix = "start")
  } else {
    getLocus(sv)
  }
  
  
  # reference sequence of requested context region. In coordinates of refGen
  refSeq <- Biostrings::getSeq(Rsamtools::FaFile(sv[["refGen"]]), param=gr)[[1L]]
  
  # if there is no overlap of region with relevant SV-locus, return unchanged reference seq
  # if the relevant SV-locus does not lie completely within the requested context genomic region:
  overlapReg <- GenomicRanges::intersect(x = relSVLocus, y = gr)
  #if ( isGainSV && GenomicRanges::countOverlaps(query = relSVLocus, subject = gr, type = "within") < 1L ){
  if ( length(overlapReg) == 0L || GenomicRanges::width(overlapReg) == 0L ){
    #logging::logwarn("Requested genomic range does not cover SV-locus. Return refSeq at requested region!")
    warning("Requested genomic range does not cover SV-locus. Return refSeq at requested region!")
    return( refSeq )
  }
  
  
  # Assumption: relevant SV-locus is completely covered by gr!! this assumption makes function easier
  stopifnot( GenomicRanges::countOverlaps(query = relSVLocus, subject = gr, type = "within") == 1L )
  
  # refSeq bases before and after the relevant SV-locus
  grLen <- GenomicRanges::width(gr)
  nbrNtBeforeSV <- GenomicRanges::start(relSVLocus) - GenomicRanges::start(gr)
  nbrNtAfterSV <- GenomicRanges::end(gr) - GenomicRanges::end(relSVLocus) + isGainSV # gain-SVs have a length 1 relSVLocus: the first base of SV-locus. This needs to be include the first
  stopifnot( nbrNtBeforeSV <= grLen, nbrNtAfterSV <= grLen )
  
  seqBeforeSV <- if (nbrNtBeforeSV >= 1L) Biostrings::subseq(refSeq, end = nbrNtBeforeSV) else Biostrings::DNAString()
  seqAfterSV <- if (nbrNtAfterSV >= 1L) Biostrings::subseq(refSeq, start = -nbrNtAfterSV) else Biostrings::DNAString()
  # nbrNtBeforeSV <- GenomicRanges::start(svLocus) - GenomicRanges::start(gr)
  # seqBeforeSV <- Biostrings::subseq(refSeq, start = 1L, end = nbrNtBeforeSV)
  # # NO! for duplication, the SV-length must be taken twice ##(isDup + 1L) *
  # seqAfterSV <- Biostrings::subseq(refSeq, start = nbrNtBeforeSV + getLength(sv) + 1L)
  
  switch(sv[["type"]],
         insertion=,
         duplication=Biostrings::xscat(seqBeforeSV, getSVSeq(sv), seqAfterSV),
         deletion=Biostrings::xscat(seqBeforeSV, seqAfterSV), # del: SVseq is left out
         inversion=Biostrings::xscat(seqBeforeSV, getSVSeq(sv), seqAfterSV) # SVseq is already refComplemented
  )

}

#length.SV <- function(x) length(x[["seq"]]) # otherwise I get trouble with functions that use length for lists (like str)
#str.SV <- function(x) NextMethod()
##seqinfo.SV <- function(x) { on.exit(close(refGenFile)); refGenFile <- open(Rsamtools::FaFile(x[["refGen"]])); return(seqinfo(refGenFile)) }



#' Find real first position of mapped reads, also considering clipped bases.
#' 
#' Care should be taken to only give reads that are really mapped (where pos and cigar string has a real meaning)
#' @param cigarStr character vector
#' @param chrom character vector for chromosome. Use to build a GRanges object.
#' @param pos position of first mapped base
#' @return start and end position of the read bases on the reference as \code{GRanges} object 
cigarToPos <- function(cigarStr, chrom, pos){
  
  stopifnot( is.character(cigarStr), is.numeric(pos) )
  stopifnot( length(cigarStr) == length(pos) )
  
  if ( length(pos) == 0L ) return( NULL )
  
  # pos refers to first *mapped* base. field in presence of soft-clipped bases in beginning of reads.
  
  # # is NOT bullet proof! e.g. "8I92M" gives NA
  # startClippedBases <- stringr::str_replace(cigarStr, "S?[[:digit:]]+[M=X].*$", "")
  # startClippedBases[!nzchar(startClippedBases)] <- "0"
  # # I do expect a zero or a single number (because there should not be I or D before S). But to be really sure (e.g. could be hard-clipping H?)
  # startClippedBases <- stringr::str_extract(startClippedBases, "[[:digit:]]+$")
  # startClippedBases <- as.integer(startClippedBases)
  
  # extract number of soft clipped reads at the beginning of read
  startClippedBases <- cigarStr %>% 
    stringr::str_extract("^[[:digit:]]+S") %>% 
    stringr::str_replace("S$", "")
  startClippedBases[is.na(startClippedBases)] <- "0"
  startClippedBases <- as.integer(startClippedBases)
  
  endClippedBases <- cigarStr %>%
    stringr::str_extract("[[:digit:]]+S$") %>%
    stringr::str_replace("S$", "")
  endClippedBases[is.na(endClippedBases)] <- "0"
  endClippedBases <- as.integer(endClippedBases)
  # But to be really sure (e.g. could be hard-clipping H?)
  #endClippedBases <- stringr::str_extract(endClippedBases, "[[:digit:]]+$")
  
  if (any(is.na(startClippedBases)) | any(is.na(endClippedBases))){
    logging::logerror("NA in cigarToPos for cigar string %s.", paste(head(cigarStr[is.na(startClippedBases) | is.na(endClippedBases)]), collapse=" - "))
    startClippedBases[is.na(startClippedBases)] <- 0L
    endClippedBases[is.na(endClippedBases)] <- 0L
  }
  
  readStartPos <- pos - startClippedBases
  # Insertions do not count for width on reference
  # Deletions in the read do count, though
  readEndPos <- pos + GenomicAlignments::cigarWidthAlongReferenceSpace(cigarStr) - 1L + endClippedBases
  
  
  GenomicRanges::GRanges(seqnames = chrom,
                         ranges = IRanges::IRanges(start=readStartPos, end=readEndPos) )
}





#' Estimate a range where reads are affected through an SV-injection.
#' 
#' We can increase the affected region via the \code{margin=} parameter. Helpful as \code{scanBam} for a given region only shows reads that really cover the region with a base (not insert/gap).
#' Except for inversion, read pairs that span the SV-locus are affected by it.
#' With \code{margin=0} we get the relevant SV-area that is relevant to decide what read (pairs) are affected:
#' \itemize{
#'   \item insertion and duplication: margin around single left (5') breakpoint (BP)
#'   \item deletion: the SV locus
#'   \item inversion: the SV locus. For inversion (balanced SV), only reads are affected that really map above it.
#' }
#' @param sv SV-object that is to be inserted
#' @param margin numeric. Number of nucleotide base pairs to increase the range on either side of the SV-region as necessary depending on SV type
getAffectedRange <- function(sv, margin=0L) {
  svGR <- sv[["locus"]]
  stopifnot( margin >= 0L )
  
  affRange <- switch(sv[["type"]],
         insertion=,
         duplication= svGR %>%
           GenomicRanges::resize(width=1L, fix="start") %>%  ## only start-point matters for gain-SVs
           GenomicRanges::resize(width=2L*margin + 1L, fix = "center"),
         inversion= svGR, ## inversion is balanced: no margin effect, only bases directly above SV-bases are affected.
         deletion= svGR %>%  # deletion: start from borders of SV locus
           GenomicRanges::resize(width=GenomicRanges::width(svGR) + 2L*margin, fix = "center"),
         stop("This SV type is not supported yet!")
  )
  
  # fix start and end coordinates
  chromEndCoord <- GenomeInfoDb::seqlengths(sv[["refGenInfo"]])[getChrom(sv)]
  GenomicRanges::start(affRange) <- pmax(1L, GenomicRanges::start(affRange))
  GenomicRanges::end(affRange) <- pmin(chromEndCoord, GenomicRanges::end(affRange))
  
  return(affRange)
}




#' Simulate an SV-injection in a given read mapping.
#' 
#' Finds reads that directly or indirectly (via insert or mate read) are affected by the SV-injection.
#' We use a hack to have discordant read pairs in \code{GAlignmentPairs}.
#' Estimate "affected reads range" where a simulation will generate new read/readpairs over the affected region that incorporate the SV.
#' The remaining unaffected mapping in the vicinity to (i.e. around) the affected region is merged together.
#' 
#' 
#' The SV and the mapping should refer to the same reference. Simulated reads follow the NGS properties found in the mapping file.
#' @param bamFile mapping file in BAM format
#' @param sv SV-object
#' @param localMappingMargin numeric. number of basepairs to go beyond the SV to consider for local mapping of unaffected reads. Could depend on NGS-properties. 75\% of it are the margin for paired-end search.
#' @param ngsProp NGS property estimated from the whole BAM file. Used for SV-read simulation to match properties.
#' @param TL numeric. tumor load proportion of the SV to be injected: what percent of DNA molecules (and also reads) carries the SV to be injected?
#' @return list with affected reads and merged mapping that carries reads of the local mapping: remaining unaffected reads and SV-simulated reads
#' @note The SV and the mapping should refer to the same reference. For the time being, only read pairs (no SE) are simulated.
injectSVInMapping <- function(bamFile, sv, localMappingMargin, ngsProp, TL){
  
  if (FALSE){
    bamFile <- "~/projects/svmod/tests/testthat/resource/test_SV0.bam"
    sv <- SV(type = "insertion", chrom = "chr5", posL = 170837725L, seqLen = 51L, refGen = file.path(BASEDIR, "refGen", "hg19_normal.fa"))
    localMappingMargin <- 1e5L #2L * EXTRA_MARGIN_TARGET_REGIONS
    ngsProp <- ngs.prop(); TL <- 0.95
  }#fi
  
  stopifnot( checkFile(bamFile), TL > 0L, TL <= 1L )
  
  
  # ranges where SV has an effect: w/ and w/o extra margin
  svAffectedRange_0 <- getAffectedRange(sv, margin = 0L)  # no extra margin
  svAffectedRange_s <- getAffectedRange(sv, margin = READL+1L)  # single
  svAffectedRange_p <- getAffectedRange(sv, margin = localMappingMargin*.75+1L)  # pairs
  
  # range for local mapping: should in any case be increased around both sides of the SV
  svAffectedRange_l <- GenomicRanges::resize(svAffectedRange_0, # local
                                             width=GenomicRanges::width(svAffectedRange_0) + 2L*localMappingMargin + 1L, fix = "center")
  
  
  # what-field for BAM-queries
  myWhat <- setdiff(Rsamtools::scanBamWhat(), c("qwidth"))  #"qual", #"seq",  # SEQ  is needed later on
  
  
  # tmp dir
  bamBaseName0 <- paste(stringr::str_replace(basename(bamFile), ".bam$", ""), toString(parseSampleProp(bamFile)), sep="_")
  bamBaseName <-  paste0(bamBaseName0, "_INJ", stringr::str_sub(sv[["type"]], end=3L), getLength(sv), "bp_", getChrom(sv),":", getPosL(sv))#, TIMESTAMP)
  tmpMappingDir <- file.path(TMPDIR, bamBaseName0)
  try(dir.create(tmpMappingDir, showWarnings = FALSE, recursive = TRUE))
  
  stopifnot( file.exists(tmpMappingDir) )
  
  
  
  
  # find reads that are *mapped* in a broad region around SV
  # restrict to reads that are really affected by the SV:
  # + all reads that have mapped bases directly above the SV-locus are in svAffectedRange_0  _or_
  # + paired reads where SV lies within left-most and right-most end of read pair (including clipped bases)
  
  # 1/ find read pairing (single reads, read pairs -concordant and discordant- via hack into GAlignmentPairs) in broad candidate region
  # 2/ extend the range of reads (as clipped bases are not considered in BAM)
  # 3/ estimate which read pairs or single reads are affected: if left end to right end crosses SV-locus
  
  
  
  
  # mapped pairs ---------
  
  # concordantly mapped ---
  paramPaired <- Rsamtools::ScanBamParam(flag=Rsamtools::scanBamFlag(isPaired=TRUE, isUnmappedQuery = FALSE, hasUnmappedMate = FALSE, 
                                                                 isProperPair=NA, isSecondaryAlignment = FALSE, isDuplicate = FALSE),
                                     what=myWhat,
                                     which=svAffectedRange_l)
  
  # affectedReads_m <- Rsamtools::scanBam(file=bamFile, param=paramPaired)[[1L]]
  # # drop supplementary reads (Rsamtools does not handle that right now)
  # primReadsInd <- which(bitops::bitAnd(affectedReads_m$flag, FLAG_SECONDARY + FLAG_SUPPL) == 0L )
  # affectedReads_m <- lapply(affectedReads_m, function(x) x[primReadsInd])
  
  # GAlignementPairs would be nice, but alas:
  # only properly paired reads on opposite strands are allowed for now (says doc of GAlignments 1.6.1) discordant read pairs are dropped!
  # read in reads as pairs.
  # strandMode=0L means: no strand information for the pair. It is set finally below for both, concordant and discordant pairs together
  # strandMode: 2L for TruSeq recommended, does it have an impact for me?
  myStrandMode <- 0L
  affectedConcPairs_m <- suppressWarnings(GenomicAlignments::readGAlignmentPairs(file=bamFile, param = paramPaired, strandMode = myStrandMode))
  # drop read pairs that contain supplementary reads (Rsamtools does not handle that right now)
  primReadPairsInd <- which(
    bitops::bitAnd(GenomicRanges::mcols(GenomicAlignments::first(affectedConcPairs_m))$flag, FLAG_SECONDARY + FLAG_SUPPL) == 0L &
      bitops::bitAnd(GenomicRanges::mcols(GenomicAlignments::last(affectedConcPairs_m))$flag, FLAG_SECONDARY + FLAG_SUPPL) == 0L )

  affectedConcPairs_m <- affectedConcPairs_m[primReadPairsInd]
  
  
  
  # discordantly mapped pairs --- (need to come extra because of GAlignmentPairs limitation)
  paramPairedPlus <- Rsamtools::ScanBamParam(flag=Rsamtools::scanBamFlag(isPaired=TRUE, isUnmappedQuery = FALSE, hasUnmappedMate = FALSE,
                                                                        isMinusStrand = FALSE, isMateMinusStrand = FALSE,
                                                                        isProperPair=NA,
                                                                        isSecondaryAlignment = FALSE, isDuplicate = FALSE),
                                            what=myWhat,
                                            which=svAffectedRange_l)
  
  paramPairedMinus <- Rsamtools::ScanBamParam(flag=Rsamtools::scanBamFlag(isPaired=TRUE, isUnmappedQuery = FALSE, hasUnmappedMate = FALSE,
                                                                        isMinusStrand = TRUE, isMateMinusStrand = TRUE,
                                                                        isProperPair=NA,
                                                                        isSecondaryAlignment = FALSE, isDuplicate = FALSE),
                                            what=myWhat,
                                            which=svAffectedRange_l)
  
  
  affectedDiscordantPairs_m <- c(GenomicAlignments::readGAlignments(file=bamFile, param = paramPairedMinus),
                                 GenomicAlignments::readGAlignments(file=bamFile, param = paramPairedPlus) )
  
  
  
  primReadDiscPairsInd <- which(bitops::bitAnd(GenomicRanges::mcols(affectedDiscordantPairs_m)$flag, FLAG_SECONDARY + FLAG_SUPPL) == 0L)
  # discordant reads listed as individual reads (not paired)
  affectedDiscordantPairs_m <- affectedDiscordantPairs_m[primReadDiscPairsInd]
  
  logging::logdebug("There are %d discordantly mapped read ends w/o secondary and supplementary read ends.", length(affectedDiscordantPairs_m))
  
  
  # mkuhn, 2015-12-29: keep only reads that occur as a pair of two read ends
  discPairReadNamesToKeep <-  local({ h <- table(mcols(affectedDiscordantPairs_m)$qname); names(h)[which(h == 2L)] })
  affectedDiscordantPairs_m <- affectedDiscordantPairs_m[which(mcols(affectedDiscordantPairs_m)$qname %in% discPairReadNamesToKeep)]
  
  # expect even number of reads (as they are from pairs)
  nbrDiscReadEnds <- length(affectedDiscordantPairs_m)
  logging::logdebug("There are %d discordantly mapped read ends w/o secondary and supplementary read ends that occur as proper read pair", nbrDiscReadEnds)
  
  stopifnot( nbrDiscReadEnds %% 2L == 0L )
  
  # set results of paired reads
  resPairs0 <- if (nbrDiscReadEnds > 0L){
    # pair reads of discordant read pairs:
    # based on flag and qname ordering
    discPairData <- mcols(affectedDiscordantPairs_m)
    discPairData$isFirstInPair <- bitops::bitAnd(discPairData$flag, FLAG_FIRST_IN_PAIR) == FLAG_FIRST_IN_PAIR
    #discPairData[order(discPairData$isFirstInPair, discPairData$qname, decreasing = FALSE),]
    affectedDiscordantPairs_m <- affectedDiscordantPairs_m[order(discPairData$isFirstInPair, discPairData$qname, decreasing = FALSE)]
    try(rm(discPairData), silent = TRUE) #cleanup
    firstHalfInd <- 1L:(nbrDiscReadEnds %/% 2L)
    discProperFlag <- bitops::bitAnd(mcols(affectedDiscordantPairs_m[firstHalfInd])$flag, FLAG_PROPER_PAIR) == FLAG_PROPER_PAIR
    # mkuhn, 2015-12-17: hack to add same columns in mcols (groupid and mate_status) in order to be able to merge it with affectedConcPairs_m
    mcols(affectedDiscordantPairs_m)$groupid <- 1L + length(affectedConcPairs_m) + seq_along(affectedDiscordantPairs_m) # fake groupIDs (with an offset)
    mcols(affectedDiscordantPairs_m)$mate_status <- factor("mated", levels = c("mated", "ambiguous", "unmated"))
    # ordering was isFirstInPair=FALSE, then isFirstInPair=TRUE. Hence, first half has reads that are 2nd in pair.
    affectedDiscordantPairs_m <- GenomicAlignments::GAlignmentPairs(first = affectedDiscordantPairs_m[-firstHalfInd], last = affectedDiscordantPairs_m[firstHalfInd],
                                                                    strandMode = myStrandMode, isProperPair = discProperFlag)
    
    # combine candidate pairs (concordant and discordant) as GAlignmentPairs
    c(affectedConcPairs_m, affectedDiscordantPairs_m)
  } else affectedConcPairs_m
  GenomicAlignments::strandMode(resPairs0) <- myStrandMode # orientation for read pairs does not matter here!?
  
  
  # look into spanned range for each read pair
  extRangeFirst <- cigarToPos(chrom = getChrom(sv), cigarStr = GenomicAlignments::cigar(GenomicAlignments::first(resPairs0)),
                              pos = GenomicAlignments::start(GenomicAlignments::first(resPairs0)) )
  extRangeLast <- cigarToPos(chrom = getChrom(sv), cigarStr = GenomicAlignments::cigar(GenomicAlignments::last(resPairs0)),
                             pos = GenomicAlignments::start(GenomicAlignments::last(resPairs0)) )
  
  # build the union of first and last read range
  extRangePair <- GenomicRanges::punion(extRangeFirst, extRangeLast, fill.gap=TRUE)
  mcols(resPairs0)$extendedRange <- extRangePair
  
  # which read pairs span the SV-locus?
  readPairHitInd <- GenomicRanges::mcols(GenomicRanges::pintersect(x=extRangePair, y=svAffectedRange_0,
                                                                   drop.nohit.ranges=FALSE, ignore.strand = TRUE))$hit
  resPairs <- resPairs0[readPairHitInd]
  
  
  
  
  
  # single mapped reads (SE or PE with unmapped mate) ---------
  
  
  # PE-reads where only one mate of pair is mapped ---
  paramPairedSingle <- Rsamtools::ScanBamParam(flag=Rsamtools::scanBamFlag(isPaired=TRUE, isUnmappedQuery = FALSE, hasUnmappedMate = TRUE,
                                                                   isProperPair=NA,
                                                                   isSecondaryAlignment = FALSE, isDuplicate = FALSE),
                                       what=myWhat,
                                       which=svAffectedRange_s)  # no big margin, only read margin 
  affectedPairsSingle <- GenomicAlignments::readGAlignments(file=bamFile, param = paramPairedSingle)
  
  
  # single-end (SE) reads ---
  paramSE <- Rsamtools::ScanBamParam(flag=Rsamtools::scanBamFlag(isPaired=FALSE, isUnmappedQuery = FALSE, hasUnmappedMate = FALSE,
                                                                           isProperPair=FALSE,
                                                                           isSecondaryAlignment = FALSE, isDuplicate = FALSE),
                                               what=myWhat,
                                               which=svAffectedRange_s)  # no big margin, only read margin 
  affectedSE <- GenomicAlignments::readGAlignments(file=bamFile, param = paramSE)
  
  # all single mapped reads, from PE and from SE combined
  affectedSingle <- c(affectedPairsSingle, affectedSE)
  affectedSingle <- affectedSingle[which(bitops::bitAnd(GenomicRanges::mcols(affectedSingle)$flag, FLAG_SECONDARY + FLAG_SUPPL) == 0L)]
  
  # extend read ranges
  extRangeSingle <- cigarToPos(chrom = getChrom(sv), cigarStr = GenomicAlignments::cigar(affectedSingle),
                               pos = GenomicAlignments::start(affectedSingle) )
  mcols(affectedSingle)$extendedRange <- extRangeSingle
  
  resSingle <- if (! is.null(extRangeSingle) && length(extRangeSingle) > 0L){
    # find hits (after correcting for clipped bases at the end)
    readSingleHitInd <- GenomicRanges::mcols(GenomicRanges::pintersect(x=extRangeSingle, y=svAffectedRange_0,
                                                                       drop.nohit.ranges=FALSE, ignore.strand = TRUE))$hit
    
    affectedSingle[readSingleHitInd]
  } else affectedSingle
  
  
  
  
  
  
  # affected and unaffected local mapping ------------
  
  
  # names of reads that are affected by the SV
  # both read pairs names  (read pairs are only named once) and single end names and mixed
  affectedReadNames_p <- unique(mcols(GenomicAlignments::first(resPairs))$qname)
  affectedReadNames_s <- unique(mcols(resSingle)$qname)
  
  # mkuhn, 2016-01-10: tumor load:
  # only TL-proportion of the affected reads get simulated from SV-injected genome
  affectedReadNames_p <- sample(affectedReadNames_p, size = ceiling(length(affectedReadNames_p) * TL))
  affectedReadNames_s <- sample(affectedReadNames_s, size = ceiling(length(affectedReadNames_s) * TL))
  
  affectedReadNames <- unique(c(affectedReadNames_p, affectedReadNames_s))
  
  
  
  # filter remaining reads of local mapping into a separate BAM file (only mapped reads are filtered, because of ScanBamParam-which)
  localBamFile <- Rsamtools::filterBam(bamFile, destination = file.path(tmpMappingDir, paste0(bamBaseName, ".local.bam")),
                                       indexDestination = TRUE,
                                       param = Rsamtools::ScanBamParam(what=myWhat,
                                                                       which=svAffectedRange_l) # only local mapping!!
  )
  
  # remainingBamFile <- Rsamtools::filterBam(bamFile, destination = file.path(tmpMappingDir, paste0(bamBaseName, ".remaining.bam")),
  #                                          param = Rsamtools::ScanBamParam(what=myWhat,
  #                                                                          which=svAffectedRange_l),  # only local mapping!!
  #                                          filter=S4Vectors::FilterRules(list(qnameFilter=function(df) return(! df$qname %in% affectedReadNames))) )
  remainingBamFile <- Rsamtools::filterBam(localBamFile, destination = file.path(tmpMappingDir, paste0(bamBaseName, ".remaining.bam")),
                                           filter=S4Vectors::FilterRules(list(qnameFilter=function(df) return(! df$qname %in% affectedReadNames))) )
  
  
  
  
  # simulate new reads in affected region ----------
  
  if (length(affectedReadNames) > 0L){
    # combine GRanges region where to simulate new reads:
    allAffectedExtRanges <- c(mcols(resSingle)$extendedRange, mcols(resPairs)$extendedRange)
    # un-listing
    if (is.list(allAffectedExtRanges) && length(allAffectedExtRanges) == 1L) allAffectedExtRanges <- allAffectedExtRanges[[1L]]
    
    if ( ! is(allAffectedExtRanges, "GRanges") || length(allAffectedExtRanges) == 0L ){
      logging::logwarn("No affected reads have been found for %s.", sv)
      return(invisible(NULL))
    }
    
    # range of all affected reads (except crazy paired reads)
    affectedExtRange <- range(allAffectedExtRanges[width(allAffectedExtRanges) <= MAX_PAIR_WIDTH])
    stopifnot( length(affectedExtRange) == 1L )
    
    # go a little beyond SV-locus
    extendedSVLocus <- GenomicRanges::resize(getLocus(sv), width = getLength(sv) + READL + 1L, fix = "center")
    # mkuhn, 2016-01-10: <branch SV2> no additional extension for duplication (because we have rectified that posL is now left of the peaks)
    #+for duplication we need to cover both duplicated regions, hence need also to look SV-length ahead of SV-locus
    # if (sv$type == 'duplication') {
    # getLocus(sv) %>% 
    #   GenomicRanges::resize(width = 2L*getLength(sv), fix = "end") %>%  #grow to the left to include the first duplicate
    #   GenomicRanges::resize(width = 2L*getLength(sv) + READL + 1L, fix = "center") } else 
    
    
    # affectedExtRange should cover extended SV-locus:
    if ( GenomicRanges::countOverlaps(query=extendedSVLocus, affectedExtRange, type="within") != 1L ){
      logging::loginfo("extended SV locus %s not (completely) covered by extended range %s of affected reads and read pairs. Extending the range to SV locus.",
                       extendedSVLocus, affectedExtRange)
      affectedExtRange <- range(affectedExtRange, extendedSVLocus)
    }#fi
    
    
    # sequencing context of injected SV ---
    refSeqSV <- getContextSeq(sv = sv, gr = affectedExtRange) # new reference seq (including injected SV) for read simulation
    fastaFileSV <- file.path(tmpMappingDir, paste0(bamBaseName,".fasta"))
    Biostrings::writeXStringSet(x = Biostrings::DNAStringSet(refSeqSV) %>% setNames(getShortName(sv)), 
                                filepath = fastaFileSV)
    
    logging::loginfo("Seq Context for SV has length %d", length(refSeqSV))
    
    
    
    # simulate given number of new reads ---
    #ZZZ mkuhn, 2015-01-06: all reads are treated as if they are paired, with uniform coverage at affectedExtRange.
    #ZZZ with same characteristics as affected or local reads (SE/PE prop., insert size, orientation distribution, same coverage pattern?)
    # Better start with something like for separate handling of single and paired reads
    #nbrReadsSV_p <- ceiling(length(affectedReadNames_p) * (length(refSeqSV) + 1L) / (GenomicRanges::width(affectedExtRange)+1L))
    ##nbrReadsSV_s <- ceiling(length(affectedReadNames_s) * (length(refSeqSV) + 1L) / (GenomicRanges::width(affectedExtRange)+1L))
    
    # new number of reads: modify it by the ratio new length versus old length. E.g. insertion will get more reads, deletion less, inversion same number
    nbrReadsSV <- ceiling( length(affectedReadNames)  *  (length(refSeqSV) + 1L) / (GenomicRanges::width(affectedExtRange)+1L) )
    
    myFastaPrefix <- paste0("reads_", getShortName(sv), "_")
    readSimulatorWrapper(fastaFile = fastaFileSV, prefixRead = "simPE", 
                         prefixFasta = myFastaPrefix, paired = TRUE, ngs.prop = ngsProp, rCount = nbrReadsSV)
    
    
    myFastqFiles <- list.files(path = dirname(fastaFileSV), pattern = paste0("^", myFastaPrefix, ".*[12][.]fq$"), full.names = TRUE)
    stopifnot( length(myFastqFiles) == 2L )
    mappingSV <- mapWrapper(inputFile1 = myFastqFiles, interleaved = FALSE, bwa_optsExtra = getISDistForBWA(ngsProp))
    
    # ZZZ check here for secondary and supplementary reads in new mapping (>> could use that to get new patMappingInfo for injected BAM vs WT BAM quickly)
    
    
    # output ---------
    
    intermediateBamFiles <- c(remainingBamFile, mappingSV)
    
    # .. and merge new bam (in tmp/) with subtracted bam of unmodified reads from local mapping
    mergedBamFile <- file.path(tmpMappingDir, paste0(bamBaseName, ".new.bam"))
    logging::loginfo("Start with merging from local mapping from %s: unaffected remaining reads and new SV reads into %s", bamFile, mergedBamFile)
    # ZZZ this can crash in some cases (low TL) # maybe fixed now with if-test for affectedReadnames?!
    try({
      newBamFile <- Rsamtools::mergeBam(files = intermediateBamFiles, destination = mergedBamFile, overwrite=TRUE)
      logging::logdebug("Merging returned %s of class %s.", toString(newBamFile), class(newBamFile))
      # index new merged BAM file
      Rsamtools::indexBam(mergedBamFile) ##newBamFile
      # clean up intermediate FASTQ/BAM files..
      # ZZZ more files to delete?!
      try(unlink(x=c(myFastqFiles, intermediateBamFiles, paste0(intermediateBamFiles, ".bai"))), silent=TRUE)
    }, silent = FALSE)
    
    
    
    
  } else {
    
    logging::logwarn("No affected reads for BAM-file %s at %s", bamFile, toString(sv))
    mergedBamFile <- remainingBamFile
  }

  list(affectedReadNames=affectedReadNames, readPairs=resPairs, readSingle=resSingle, localBam=localBamFile, mergedBam=mergedBamFile)
}




