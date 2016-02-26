


#' Explore distribution of clipped bases across a mapping at a given region. Uses the BAM coordinates of read segments (left-most mapped base).
#' 
#' This function helps to explore how a region is actually covered.
#' It produces a plot showing the region of interest.
#' 
#' @param bamFile a BAM file. Typically, this will be a mapping from a (simulated) tumor sample.
#' @param windowSize integer value for the width of the rolling average window
#' @param startPos start coordinate for region
#' @param endPos end coordinate for region
#' @param SVstart offset to SV from startPos
#' @param focus.SV logical flag if we want to focus around the given breakpoint (BP) of the SV. Can also be an integer which already gives the width of the focus area. \code{TRUE} corresponds to 500 but minimally the length of the SV.
#' @param add.nbrQ1Reads flag if number of Q1-reads mapped at position should be added to the plot
#' @export
exploreClipping_BAM <- function(bamFile, chrom="chr5", startPos=1234567, endPos=startPos+500L, SVtype=NULL, SVlength=0, SVstart=0, 
                                windowSize=20L, all.Sbases=FALSE, add.nbrQ1Reads=TRUE, focus.SV=FALSE){
  
  
  if ( ! checkFile(bamFile) ){
    logerror("BAM file %s not found!", bamFile)
    return(invisible(FALSE))
  }
  
  
  
  GENOMIC_REGION <- paste0(chrom, ":", startPos,"-", endPos)
  loginfo("Look into mapping at region %s at BAM-file %s.", GENOMIC_REGION, basename(bamFile))
  
  #rga <- readGAlignmentsFromBam(bamFile)  #GenomicAlignments::
  
  # all paired reads that map into the genomic region of interest
  myParam <- ScanBamParam(flag=scanBamFlag(isPaired=TRUE, isUnmappedQuery = FALSE, isProperPair=NA, isNotPrimaryRead=FALSE),
                              what=setdiff(scanBamWhat(), c("qwidth")),  #"qual", #"seq",  #include SEQ for convenience with IGV
                              which=GRanges(seqnames=chrom, ranges=IRanges(start=startPos, end=endPos)) )
  
  mapping <- data.frame(scanBam(file=bamFile, param=myParam)[[1]])
  # filter on Q1 reads
  mapping <- mapping[mapping$mapq > 0L,]
  
  # sorted BAM file
  stopifnot( all( diff(mapping$pos) >= 0L ) )
  
  # no hard-clipping  (see BWA MEM -Y option)
  stopifnot( ! any(grepl("H", mapping$cigar)) )
  
  
  # sum of clipped bases within the Q1-reads, from either end
  mapping$Sbases_5p <- as.numeric(str_sub(str_extract(mapping$cigar, pattern = "^[[:digit:]]+S"), end=-2L))
  mapping$Sbases_3p <- as.numeric(str_sub(str_extract(mapping$cigar, pattern = "[[:digit:]]+S$"), end=-2L))
  
  mapping$Sbases_5p[is.na(mapping$Sbases_5p)] <- 0L
  mapping$Sbases_3p[is.na(mapping$Sbases_3p)] <- 0L
  
  loginfo("There are %d 5p (start) and %d 3p (end) S-bases in bamfile ::%s::", 
          sum(mapping$Sbases_5p, na.rm=TRUE), sum(mapping$Sbases_3p, na.rm=TRUE), basename(bamFile))
  
  SV.BP <- startPos + SVstart
  
  #mapping.sv <- filter(mapping, pos >= SV.BP - 120L, pos <= SV.BP + 120L)
  
  
  # aggregating per start position of mapped reads
  mappingCum <- group_by(mapping, pos) %>% summarise(nbrQ1Reads=n(), meanMapQ=mean(mapq, na.rm = TRUE), Sbases_5p=sum(Sbases_5p), Sbases_3p=sum(Sbases_3p))
  zeroPos <- setdiff(x=startPos:endPos, y=mappingCum$pos) # all positions that are no start position of mapped reads
  zeroMapping <- data.frame(pos=zeroPos, nbrQ1Reads=0, meanMapQ=NA, Sbases_5p=0, Sbases_3p=0)
  mappingCum <- rbind(mappingCum, zeroMapping)
  mappingCum <- arrange(mappingCum, pos) #sort by left-most position
  
  hammWeights <- getHammingWindowWeights(windowSize)
  #ZZZ I would like more smoothed curve for q1ReadPos, but longer window size leads to different length
  #windowSize.cov <- 2L * windowSize 
  #hammWeights.cov <- getHammingWindowWeights(windowSize.cov) 
  
  # filtered cumulated mapping data (with Hamming windows)
  mappingCumF <- data.frame(pos=roll_mean(mappingCum$pos, n=windowSize, weights = hammWeights),
                         nbrQ1Reads=roll_mean(mappingCum$nbrQ1Reads, n=windowSize, weights = hammWeights),
                         #meanMapQ=roll_mean(mappingCum$meanMapQ, n=windowSize, weights = hammWeights),
                         Sbases_5p=roll_mean(mappingCum$Sbases_5p, n=windowSize, weights = hammWeights),
                         Sbases_3p=roll_mean(mappingCum$Sbases_3p, n=windowSize, weights = hammWeights),
                         Sbases=roll_mean(mappingCum$Sbases_5p + mappingCum$Sbases_3p, n=windowSize, weights = hammWeights)
  )
  #mappingCumF$Sbases <- mappingCumF$Sbases_5p + mappingCumF$Sbases_3p
  
  if (is.logical(focus.SV)) focus.SV <- if (isTRUE(focus.SV)) max(500L, SVlength+3L) else 0L
  myXLim <- if ( focus.SV > 0L && SVlength > 0L ) SV.BP + c(-1, 1) * focus.SV else NULL 
    #as.numeric(c(SV.BP - SVlength - 50L, SV.BP + SVlength + 50L)) else NULL
  
  if ( isTRUE(all.Sbases)){
    plot( Sbases ~ pos, data=mappingCumF, type="l", ylab="accum. nbr clipped bases", xlim=myXLim )
  } else {
    plot( Sbases_5p ~ pos, data=mappingCumF, type="l", col="darkblue", ylab="accum. nbr clipped bases",
          xlim=myXLim, ylim=c(0.9, 1.05)*range(c(mappingCumF$Sbases_5p, mappingCumF$Sbases_3p)) )
    lines(Sbases_3p ~ pos, data=mappingCumF, type="l", lty=2, lwd=1.5, col="darkgreen")
    
    #     legend("topleft", c("5p (start)", "3p (end)"), lty=1:2, lwd=2, col = c("darkblue", "darkgreen"), title="Clipped", seg.len=3,
    #            cex=1, inset = c(0.01, 0.015))
  }
  
  if ( isTRUE(add.nbrQ1Reads)){
    par(new=TRUE)
    plot( nbrQ1Reads ~ pos, data=mappingCumF, type="l", axes=FALSE, ylab="", col="darkgrey", lwd=2, xlim=myXLim)
    axis(side = 4)
  }
  

  # add legend (if necessary)
  if ( ! isTRUE(all.Sbases) || isTRUE(add.nbrQ1Reads) ){
    if ( ! isTRUE(add.nbrQ1Reads)){
      stopifnot( ! isTRUE(all.Sbases) )
      legend("topleft", c("5p (start)", "3p (end)"), lty=1:2, lwd=2, col = c("darkblue", "darkgreen"), title="Clipped", seg.len=3,
           cex=1, inset = c(0.01, 0.015))
    } else if ( isTRUE(all.Sbases) ){
      legend("topleft", c("Clipped Bases", "Q1-Reads"), lty=1, lwd=2, col=c("black", "darkgrey"), title="Number Of", inset = c(0.01, 0.015) )
    } else {
      stopifnot( ! isTRUE(all.Sbases) )
      legend("topleft", c("Clipped 5p (start)", "Clipped 3p (end)", "Q1-Reads"), lty=c(1,2,1), lwd=2, col = c("darkblue", "darkgreen", "darkgrey"), title="Number Of", seg.len=3,
             cex=1, inset = c(0.01, 0.015))
    }
  }

  title(sub = paste("Window Size:", windowSize), cex.sub=0.8)
  
  if ( ! is.null(SVtype) && SVtype != 'no' ){
    title( paste0("Clipped Bases -- ", SVtype, "\n", basename(bamFile), if (focus.SV && SVlength > 0) paste0(" BP ",SV.BP, " (", SVlength, "bp)")) )
    abline(v=SV.BP, col="red", lty=4)
    if ( SVtype %in% c('inversion', 'deletion', 'duplication') ){
      abline(v=SV.BP - SVlength, col="grey", lty=4)
      abline(v=SV.BP + SVlength, col="grey", lty=4)
    }
  } else {
    title( paste("Clipped Bases -- no SV\n", basename(bamFile)))
  }
  
}





#' Explores distribution of clipped bases across a mapping at a given region. It uses the genomic coordinates of the involved bases. (_G for genomic).
#' 
#' This function helps to explore how a given region is covered by clipped bases.
#' It produces a plot showing the region of interest.
#' @param SVtype SV type that was simulated
#' @param SVlength SV length put into the simulation
#' @param SVstart SV start relative to target region start as put into the simulation
#' @param focusSV flag. Should the SV coordinates be given in the plot title?
#' @param showPeak flag. Should the predicted SV-region be highlighted?
#' @return Dataframe with 5p and 3p coordinates of candidate regions. As a side-effect it plots the clipped base distribution for the given genomic region.
#' @export
exploreClipping_G <- function(bamFile, chrom="chr5", startPos=1234567, endPos=startPos+500L, SVtype=c('no', 'inversion', 'deletion', 'duplication', 'insertion'), SVlength=0, SVstart=0,
                              windowSize=20L, focusSV=FALSE, showPeak=FALSE){
  SVtype <- match.arg(SVtype)
  
  # for quick run
  if (FALSE){
    chrom <- "chr5"
    bamFile <- file.path(getVirtualPatientPath(what="mapping", .sample.prop = sample.prop(cov=40L, tumorLoad=33)), "pat00048_T.bam")
    startPos <- 34655332; endPos <- 34656580; SVstart <- 219; SVlength <- 248
    SVtype <- "duplication"
    windowSize <- 20L; focusSV=TRUE; showPeak <- TRUE
    
    bamFile <- file.path(getVirtualPatientPath(what="mapping", .sample.prop = sample.prop(cov=40L, tumorLoad=90)), "pat00046_T.bam")
    startPos <- 80603764; endPos <- 80605857; SVstart <- 1017; SVlength <- 239
    SVtype <- "duplication"
    windowSize <- 20L; focusSV=TRUE; showPeak <- TRUE
    
    # insertion at the beginning
    bamFile <- file.path(getVirtualPatientPath(what="mapping", .sample.prop = sample.prop(cov=40L, tumorLoad=90)), "pat00050_T.bam")
    startPos <- 159833482; endPos <- 159835622
    SVtype <- "insertion"; SVstart <- 22L; SVlength <- 227L
    windowSize <- 20L; focusSV=TRUE; showPeak <- TRUE
    
    bamFile <- file.path(getVirtualPatientPath(what="mapping", .sample.prop = sample.prop(cov=40L, tumorLoad=90)), "pat00021_T.bam")
    startPos <- 175991672; endPos <- 175993744
    SVtype <- "insertion"; SVstart <- 620L; SVlength <- 240L
    windowSize <- 20L; focusSV=TRUE; showPeak <- TRUE
  }
  
  if ( ! checkFile(bamFile) ){
    logerror("BAM file %s not found!", bamFile)
    return(invisible(FALSE))
  }
  
  
  GENOMIC_REGION <- paste0(chrom, ":", startPos,"-", endPos)
  SV.BP <- startPos + SVstart
  loginfo("Look into mapping at region %s at BAM-file %s.", GENOMIC_REGION, basename(bamFile))
  
  clippedBasePos <- getClippedBasePositions(bamFile, chrom, targetStartPos = startPos, targetEndPos = endPos)
  
  # windowSize=31L, minThresholdAbs=2L, minThreshold=0.05, thresholdQuantile=0.95){ #
  clipAnalysis <- assessClippedBasesInRegion(bamFile=bamFile, chrom=chrom, .targetStartPos=startPos, .targetEndPos=endPos) #, .clippedBasePos = clippedBasePos)
  #clippedBasePos <- clipAnalysis[["clippedBasePos"]]
  clippPeaks <- clipAnalysis[["peakInfo"]] #[["peakLocation"]]
  nbrPeaks <- clippPeaks[["nbrPeaks"]]
  
  loginfo("Found %d peaks in clipped bases.", nbrPeaks)
  loginfo("Found %d entries for clippedBasePos.", NROW(clippedBasePos))
  
  # with(clippedBasePosition, matplot(x=pos, y=cbind(Sbases_5p, Sbases_3p), pch=16, cex=0.7))
  # with(clippedBasePosition, matplot(x=pos, y=sqrt(cbind(Sbases_5p, Sbases_3p)), pch=16, cex=0.7))
  # normalized
  with(clippedBasePos, matplot(x=pos, y=cbind(relSbases_5p, relSbases_3p), pch=16, cex=0.7, ylim=c(0, max(0.5, relSbases_5p, relSbases_3p)*1.05), ylab="Clipped Bases Proportion", bty='L'))
  legend(x=max(clippedBasePos$pos), y=0.4, legend = c("left", "right"), pch=16, cex=0.7, col=palette()[1:2], title = "Clipping", xpd = TRUE)
  
  if ( ! is.null(SVtype) && SVtype != 'no' ){
    title( paste0("Clipped Bases -- ", SVtype, "\n", basename(bamFile), if (focusSV && SVlength > 0) paste0("  BP ",SV.BP, " (", SVlength, "bp)")) )
    abline(v=SV.BP, col="darkblue", lty=4)
    if ( SVtype %in% c('inversion', 'deletion', 'duplication') ){
      #abline(v=SV.BP - SVlength, col="grey", lty=4)
      abline(v=SV.BP + SVlength, col="grey", lty=4)
    }
  } else {
    title( paste("Clipped Bases -- no SV\n", basename(bamFile)))
  }
  
  # show predicted clipped peak area (if showPeak)
  if (isTRUE(showPeak) && nbrPeaks > 0L){
    # for (i in seq_len(nbrPeaks)){
    #   polygon(x=rep(c(clippPeaks$pos.5p[i], clippPeaks$pos.3p[i]), each=2L), y=c(-0.02, 0.01, 0.01,-0.02), col = "lightblue")
    # }#rof
    polygon(x=rep(c(clippPeaks[["pos.5p"]], clippPeaks[["pos.3p"]]), each=2L), y=c(-0.02, 0.01, 0.01,-0.02), col = "lightblue")
  }
  return(clippPeaks)
}

