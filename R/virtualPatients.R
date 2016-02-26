
#' Applies a random SV of certain type and given length to an existing DNA sequence.
#' 
#' The SV start position is either given or randomly chosen so that it fully fits into the given DNA sequence.
#' This function can be useful when simulating a virtual patient with structural variations (SV).
#' 
#' In case of insertion the function inserts a random piece of DNA and the returned DNA string is \code{SVlength} bases longer.
#' A duplication will take the sequence starting at \code{SVstartPos} and duplicate it. \code{SVseq} will be the duplicated sequence.
#' For duplication, \code{SVstartPos} points to the first base of the DNA that comes twice.
#' A deletion will take out \code{SVlength} nucleotides starting at position \code{SVstartPos}. \code{SVseq} consists of the deleted sequences and
#' for deletion, \code{SVstartPos} points at the first base just before the deletion.
#' An inversion will not change the length of the DNA string. The inverted piece of DNA is given as \code{SVseq}.
#'
#' There is also the SV-S3-class that is used for SV-injection and it is based on reference genome. Here, we can provide the DNA-context
#' which can be useful to reflect also SNVs that are present. But the code would be simpler if we used only one version (the SV-class).
#'
#' @author mkuhn, Feb 2014
#' @param DNAseq the original DNA sequence as DNAString object
#' @param SVtype type of SV
#' @param SVstartPos starting position of SV in given DNA. The SV should fully fit in the given DNA. If NULL then a position will be sampled.
#' @param SVlength number of nucleotides that form the SV or \code{NULL} for auto-selection
#' @param SVlength_mean mean length of SVs aimed at (if target region allows for it). Default value is 175.
#' @param SVlength_sd standard deviation of SVs length aimed at (if target region allows for it). Default value is 15.
#' @param margin number of nucleotides that are guaranteed to be beyond the pure SV on either side. Important to have enough seq-context (mkuhn, 2016-02-05)
#' @return list describing the injected SV: the resulting modified DNA sequence \code{DNAseq}, the pure SV sequence \code{SVseq} et al. entries
injectMySVs <- function(DNAseq, SVtype=SV_TYPES, SVstartPos=NULL, SVlength=NULL, SVlength_mean=175L, SVlength_sd = 15L, margin=2L*READL){
  SVtype <- match.arg(SVtype)
  
  stopifnot( SVlength >= 1L )
  
  # input sequence must be long enough to incorporate a minimal deletion SV (that gets reduced)
  seqLen <- length(DNAseq)
  stopifnot( seqLen >= 2L * margin + SV_MIN_LENGTH + 1L )
  
  # mkuhn, 2016-02-06: choose an appropriate SVlength if not manually given
  if (is.null(SVlength) || SVlength < SV_MIN_LENGTH){
    # first try around aimed mean SV length
    SVlength <- min(SV_MAX_LENGTH, max(SV_MIN_LENGTH, round(rnorm(n = 1L, mean = SVlength_mean, sd = SVlength_sd))))
    # if SV too long than sample from what is available
    if ( SVlength >= seqLen - 2L * margin )
      SVlength <- min(SV_MAX_LENGTH, SV_MIN_LENGTH + sample(seqLen - 2L * margin - 1L - SV_MIN_LENGTH, size=1L))
    logging::loginfo("Chose SV-length of %d.", SVlength)
  }
  
  stopifnot( SVlength >= SV_MIN_LENGTH )
  stopifnot( SVlength < seqLen - 2L * margin )
  logging::loginfo("About to generate an SV of type %s of length %d into a sequence of length %d.", SVtype, SVlength, seqLen)
  
  # sample a start position. Not too close to the end (there was once a error message from replaceAt)
  SVstartPos <- if (is.null(SVstartPos)){
    logging::logdebug("Sampling start position of SV.")
    margin + sample(seqLen - SVlength - 2L * margin, size=1L) } else SVstartPos
  
  # relevant for deletion or inversion that we have enough sequence available
  stopifnot( SVstartPos + SVlength <= seqLen )  
  
  SVseq <- NULL
  
  switch(SVtype, 
         insertion={
           SVseq <- Biostrings::DNAString(paste(sample(Biostrings::DNA_BASES, SVlength, replace=TRUE), collapse=""))
           logging::loginfo("Insertion of %s at position %d out of %d.", SVseq, SVstartPos, seqLen)
           DNAseq <- Biostrings::replaceAt(DNAseq, at=SVstartPos, value=SVseq)
         },
         deletion={
           delRange <- IRanges::IRanges(SVstartPos, width=SVlength)
           SVseq <- Biostrings::extractAt(DNAseq, at=delRange)[[1]]
           logging::loginfo("Deletion of ::%s:: at position %d out of %d.", SVseq, SVstartPos, seqLen)
           DNAseq <- Biostrings::replaceAt(DNAseq, at=delRange)
           SVstartPos <- max(1L, SVstartPos - 1L) # mkuhn, 20160205: point just before deletion
           #logdebug("Deleted sequence is class %s has length %d.", class(SVseq), length(SVseq))
         },
         duplication={
           SVseq <- Biostrings::subseq(DNAseq, start=SVstartPos, width=SVlength)
           logging::loginfo("Duplication of %s at position %d out of %d.", SVseq, SVstartPos, seqLen)
           DNAseq <- Biostrings::replaceAt(DNAseq, at=SVstartPos, value=SVseq)
         },
         inversion={
           SVseq  <- Biostrings::reverseComplement(Biostrings::subseq(DNAseq, start=SVstartPos, width=SVlength))
           logging::loginfo("Inverted seq %s at position %d out of %d.", SVseq, SVstartPos, seqLen)
           DNAseq <- Biostrings::replaceAt(DNAseq, at=IRanges(SVstartPos, width=SVlength), value=SVseq )
         },
         stop("Not supported SV type: ", SVtype)
  )
  
  logging::loginfo("DNA after SV-injection has length %d. SV-start position is %d.", length(DNAseq), SVstartPos)
  
  return( list(DNAseq=DNAseq, SVtype=SVtype, SVlength=SVlength, SVstartPos=SVstartPos, SVseq=SVseq) )
}






#' Generates the genetic basis of a virtual patient. This is the first thing to do in order to have a simulation run.
#' 
#' It selects a collection of k genetic loci from a target file on a single fixed chromosome and generates corresponding genomic tumor and normal sequences.
#' These sequences are written out in two multiple FASTA file, one for normal sample and one for tumor sample.
#' 
#' The germline sequence (normal probe) of the patient contains private SNVs, but -for now- no SVs.
#' The tumor sequence starts off as a copy of the germline sequence but then a number of SNVs and SVs are added on top.
#' 
#' One patient differs from another patient only by the given parameters (rate of diffSVs, germline and somatic SNV rate).
#' If in a simulation run these parameters are not changed, the different patients are all the same and are actually like one patient.
#' 
#' @param patId the ID for the patient to be created. Defaults to 1. 
#' @param baseDir base directory of simulations
#' @param chrom chromosome, defaults to \code{CHROM}=chr5
#' @param nbrLoci number of requested loci where to record normal and mutated tumor sequence. Defaults to 10. If \code{NULL} or 0 or negative, use all available regions.
#' @param nbrDSV number of loci with differential SV. Defaults to 1.
#' @param minMargin integer minimal margin left and right of required regions
#' @param germSNVRate germline SNV mutation rate. The SNPs are superimposed on the (ref) sequence. It can happen that e.g. A => A. Hence, the net germline SNP rate is lower.
#' Default value of 0.0035 corresponds to a rate of 2.6 SNP per kilobase (based on dbSNP v129)
#' @param somSNVRate somatic SNV rate. Depends on cancer. For AML, Mardis (2009) reports 750 somatic point mutations that correponds here to 2e-5 (exome size 47.9Mb, factor 4/3)
#' @param SVlength_mean average SV-length aimed at for this patient. The actual length is drawn from a normal distribution or it is sampled from what is available
#' @param SVlength_sd standard deviation of SV length for this patient. Default value is 15.
#' @param targetBedFile path to BED-file that contains the coordinates of the target regions. Default is a TruSeq target file for hg19.
#' @return summary of simulated patient loci (where is a DSV). ZZZ Better as GRanges object?!
#' @export
createPatient <- function(patId=1L, nbrLoci=10L, nbrDSV=1L, germSNVRate=0.0035, somSNVRate=2e-5, SVlength_mean=175L, SVlength_sd=15L,
                          baseDir=BASEDIR, chrom = CHROM, minMargin=4L * READL,
                          genomeRef=file.path(baseDir, "refGen", "hg19_chr5.fa"),
                          targetBedFile=file.path(baseDir, "refGen", "TruSeq_exome_targeted_regions.hg19.bed")){
  
  
  # multi sequence FASTA file of enriched regions (+flanking bases) on chr5
  stopifnot( checkFile(genomeRef), checkFile(targetBedFile) )
  
  
  
  # Select random target regions ------
  
  # form enrichment GRanges object (from BED file)
  targetBed <- read.delim(file = targetBedFile, header = FALSE)
  stopifnot( NCOL(targetBed) == 6L )
  colnames(targetBed) <- c("chr", "start", "end", "id", "width", "strand")
  # BED-file is 0-based, GRanges are 1-based
  targetBed$start <- targetBed$start + 1L
  targetRegions <- GenomicRanges::makeGRangesFromDataFrame(targetBed, ignore.strand = TRUE, keep.extra.columns = FALSE)
  try(rm(targetBed), silent = TRUE)
  
  # use only a single fixed chromosome (defaults to CHROM=chr5) 
  targetRegions.chr <- targetRegions[GenomicRanges::seqnames(targetRegions) == chrom] %>% 
    GenomicRanges::reduce(drop.empty.ranges=TRUE, ignore.strand=TRUE, min.gapwidth=READL) #READL iso 1L
  
  # enforce some minimum length (therefore we do not widen the target regions and stay more close to "real")
  targetRegions.chr <- targetRegions.chr[GenomicRanges::width(targetRegions.chr) >= 3L * READL]
  try(rm(targetRegions), silent = TRUE)
  
  # mkuhn, 2016-02-06: extend target regions to have enough sequence context to have reads mapping left and right from target region
  # for a deletion the remaining DNA will be shortend by the SV-length, therefore we enforce a remaining target region of at least 2x (left and right) the margin
  targetRegions.chr.extended <- targetRegions.chr %>%
    GenomicRanges::resize(width = pmax(2L * minMargin + SV_MIN_LENGTH + 1L, GenomicRanges::width(.)), fix = "center")
  
  # targetRegions.chr <- targetRegions.chr[GenomicRanges::width(targetRegions.chr) > minTargetWidth] %>%
  #  GenomicRanges::reduce(drop.empty.ranges=TRUE, ignore.strand=TRUE, min.gapwidth=1L)
  #try(rm(targetRegions, targetRegions.chr), silent = TRUE)
  
  # # widen the intervals and reduce the regions (i.e. merge overlapping regions)
  # targetRegions.chr.reduced %<>% 
  #   # extend left and right to be sure to have enough sequence context for reads
  #   GenomicRanges::resize(width = GenomicRanges::width(.) + 2L * READL + 1L, fix = "center") %>% 
  #   GenomicRanges::reduce(drop.empty.ranges=TRUE, ignore.strand=TRUE, min.gapwidth=1L)
  
  
  
  nbrTargetRegions <- length(targetRegions.chr)
  
  # nbrLoci
  nbrLoci <- if ( is.null(nbrLoci) || ! is.numeric(nbrLoci) || nbrLoci <= 0L ){
    nbrTargetRegions
  } else {
    min(nbrLoci, nbrTargetRegions)
  }
  
  
  # nbrDSV (differential SV)
  stopifnot( is.numeric(nbrDSV), nbrDSV > 0L )
  if ( isTRUE(nbrDSV < 1L) ){
    nbrDSV <- ceiling(nbrDSV * nbrLoci)
  }
  
  stopifnot( nbrTargetRegions >= nbrLoci )
  stopifnot( nbrDSV <= nbrLoci )
  
  
  # random sample of nbrLoci target regions
  targetRegions.ind <- sample(nbrTargetRegions, size=nbrLoci)
  myTargets <- targetRegions.chr[targetRegions.ind]
  myTargets.extended <- targetRegions.chr.extended[targetRegions.ind]
  

    
  
  # Normal probe genome -----
  
  # start with baseline (normal) genome, sample random positions from target enrichment file
  # DNAStringSet object
  patLoci.N <- Rsamtools::getSeq(open(Rsamtools::FaFile(genomeRef)), param=myTargets.extended)
  names(patLoci.N) <- sapply(myTargets.extended, toString)
  
  patData <- data.frame(patId=patId, patIdStr=getPatStr(patId), chrom=chrom,
                        # 2016-02-11: detected a bug: I had used for targetStart/End-Pos myTargets
                        #+ (not extended, therefore the SV could lie outside of targetStartPos-targetEndPos)
                        targetStartPos=GenomicRanges::start(myTargets.extended),
                        targetEndPos=GenomicRanges::end(myTargets.extended), 
                        stringsAsFactors=FALSE )
  
  # default values for all loci to start with.
  #patData$patId <- rep(getPatStr(patId), nbrLoci)
  patData$status <- rep("germline", nbrLoci)
  patData$SVtype <- "no"
  patData$SVlength <- 0L
  patData$SVstart  <- 0L
  patData$SVseq <- ""
  
  # add germline point variations
  logging::logdebug("Patient %s [N]: Adding random SNVs.", patId)
  for (i in seq_along(patLoci.N)){
    patLoci.N[[i]] <- injectRandomSNVs(patLoci.N[[i]], nbrSNVs=rbinom(1, size=length(patLoci.N[[i]]), prob=germSNVRate))
  }#rof
  
  #possibly add germline SVs to normal reads here!? The method aims for *somatic* SVs, so germline SVs should not confirm!
  
  
  
  # Tumor probe genome -------
  
  # Tumor stems from normal sample
  #+start T as identical copy of N (including the germline SNVs)
  patLoci.T <- Biostrings::DNAStringSet(patLoci.N)
  
  
  # mkuhn, 20140624
  # add some extra tumor point mutations
  logging::logdebug("Patient %s [T]: Adding extra somatic SNVs in tumor probe genome.", patId)
  for ( i in seq_along(patLoci.T) ){
    patLoci.T[[i]] <- injectRandomSNVs(patLoci.T[[i]], nbrSNVs=rbinom(1, size=length(patLoci.T[[i]]), prob=somSNVRate))
  }#rof
  
  # put differential SV in the tumor genome
  # and record it (the simulated truth)
  logging::logdebug("Patient %s [T]: Adding %d differential SVs.", patId, nbrDSV)
  # SV with minimum length 
  svTypes <- sample(SV_TYPES, size = nbrDSV, replace = TRUE)
  #svLengths <- pmin(SV_MAX_LENGTH, pmax(SV_MIN_LENGTH, ceiling(rnorm(nbrDSV, mean=SVlength.mean, sd=SVlength.sd))))
  
  for (j in seq_len(nbrDSV) ){#along(svTypes) ){  #sample.int(length(patLoci.T), size = nbrDSV) ){ #
    #mySVtype <- "insertion" # fixed insertion here
    mySVtype <- svTypes[j]
    
    # j-th tumor SV
    tumorMut <- injectMySVs(DNAseq = patLoci.T[[j]], SVtype=mySVtype, SVstartPos = NULL, SVlength = NULL,
                            SVlength_mean = SVlength_mean, SVlength_sd = SVlength_sd, margin = minMargin) #, SVlength=svLengths[j]) 
    patLoci.T[[j]] <- tumorMut[["DNAseq"]] # overwrite tumor sequence
    logging::logdebug("%d th SVseq abbreviated is %s.", j, abbrevSeq(tumorMut$SVseq))
    names(patLoci.T)[[j]] <- paste0( names(patLoci.T)[[j]],"|",tumorMut$SVtype, "|len", tumorMut$SVlength, "|pos", tumorMut$SVstartPos,
                                    if (! is.null(tumorMut$SVseq)) paste0("|seq", abbrevSeq(tumorMut$SVseq)) )
    
    patData[j, c("status", "SVtype", "SVseq")] <- c("DSV", mySVtype, toString(tumorMut$SVseq)) #character stuff
    patData[j, c("SVlength", "SVstart")] <- c(tumorMut$SVlength, tumorMut$SVstartPos) #numeric stuff
  } #rof j
  
  
  
  # Write out FASTA data per patient -----
  
  # Tumor (T) and Normal (N)
  patFastaFile <- paste0(baseDir, "patient/", getPatStr(patId))
  loginfo("Patient %d: Write genomic information to files %s.", patId, patFastaFile)
  writeXStringSet(patLoci.N, filepath=paste0(patFastaFile, "Loci_N.fasta"))
  writeXStringSet(patLoci.T, filepath=paste0(patFastaFile, "Loci_T.fasta"))
  
  # return simulated truth
  return( patData )
}




#' Generates a pool of reads for a virtual patient, the so-called "read store".
#' 
#' The functions expects to find FASTA files that have the sequence information of patient's normal and tumor sample.
#' The sequencing characteristics are given by a \code{ngs.prop} object.
#' 
#' The read store should contain enough reads to find as many reads as needed to mix virtual samples from it.
#' Per patient four FASTQ files are generated (N_1.fq, N_2.fq, T_1.fq, T_2.fq)
#' The read names should not become too long. It seems this can cause trouble in the BAM file later.
#' Hence, the sequence information for mutated reads is not put fully in the read name any more.
#' 
#' @param patId the patient ID for which to put reads in the read store.
#' @param readCov the coverage of the read store. This should be high enough to take out as many reads as needed
#' @param myNGS.prop the sequencing properties that are valid for that store
#' @param patDir directory where the patient FASTA files are stored
#' @param outDir directory where the sequencing read files of the read store are to be saved
#' @return logical status invisibly.
#' @export
generatePatientReadStore <- function(patId, readCov=POOL_COV, myNGS.prop=ngs.prop(), baseDir=BASEDIR){
  
  # Patient FASTA file that contains genetic information of different loci per patient.
  patBaseFasta <- paste0(getVirtualPatientPath(basedir = baseDir, what="patient"), getPatStr(patId), "Loci_")
  patFasta_N <- paste0(patBaseFasta, "N.fasta")
  patFasta_T <- paste0(patBaseFasta, "T.fasta")
  
  logging::loginfo("Generate read store for patient %d from files %s and %s.", patId, patFasta_N, patFasta_T)
  
  if ( ! all(file.exists(patFasta_N, patFasta_T)) || ! all(checkFile(c(patFasta_N, patFasta_T))) ){
    logging::logwarn("FASTA files for patient %s not found.", patId)
    return(invisible(FALSE))
  }
  
  if (readCov < POOL_COV) logging::logwarn("Coverage for patient read store is set to %d. The default value is %d.", readCov, POOL_COV)
  
  ###
  ### Simulate reads
  ###
  
  #-p : paired
  #-na: noALN
  #-d : id
  #-l : -len
  #-f : fcov
  #-i : in
  #-o : out
  #-m : -mflen, mean insert length
  #-s : -sdev
  
  readstoreBase <- getVirtualPatientPath(basedir = BASEDIR, what = "readstore", .ngs.prop = myNGS.prop)
  if ( isTRUE(dir.create(readstoreBase, showWarnings=FALSE, recursive=TRUE)) ) logging::loginfo("Created readstore directory %s.", readstoreBase)
  
  patReads_N <- paste0(readstoreBase, getPatStr(patId), "_N_")
  patReads_T <- paste0(readstoreBase, getPatStr(patId), "_T_")
  
  
  logging::loginfo("Start with simulation of reads for readstore.")
  
  # ..for WT sequence
  simWTRet <- system(paste(ART_EXE, "--paired --noALN --id WT --in", patFasta_N, "--len", myNGS.prop$read.length,
                           "--fcov", readCov, "--mflen", myNGS.prop$pe.ins.mean, "--sdev", myNGS.prop$pe.ins.sd, "--out", patReads_N), intern=FALSE)
  if ( simWTRet == 0L ) logging::loginfo("WT Reads simulated") else { logging::logerror("WT Read Simulation FAILED"); return(invisible(FALSE)) }
  
  # ..for MUT sequence
  simMUTRet <- system(paste(ART_EXE, "-p -na --id MUT --in", patFasta_T, "--len", myNGS.prop$read.length,
                            "--fcov", readCov, "--mflen", myNGS.prop$pe.ins.mean, "--sdev", myNGS.prop$pe.ins.sd, "--out", patReads_T), intern=FALSE)
  if ( simMUTRet == 0L ) logging::loginfo("MUT Reads simulated") else { logging::logerror("MUT Read Simulation FAILED"); return(invisible(FALSE)) }
  
  return(invisible(TRUE))
}



#' Concatenate two text files.
#' 
#' This is useful when merging different FASTA files. 
#' Deletes the old two files upon success of concatenation.
#' It uses the UNIX-OS cat command.
#' 
#' @param f1 1st file to concatenate
#' @param f2 2nd file to concatenate
#' @param f3 Name of resulting new file
appendFiles <- function(f1, f2, f3){
  if ( ! checkFile(f1) || ! checkFile(f2)){
    logwarn("File %s or file %s not found.", f1, f2)
    return(1L)
  }
  
  
  catRet <- system(paste("cat ", f1, f2, "> ", f3), intern=FALSE)
  
  if (catRet == 0){
    logging::logdebug("Deleting two source files to concat.")
    file.remove(c(f1, f2))
  } else {
    logging::logwarn("Concatenating file %s and %s failed!", f1, f2)
  }
  
  return(invisible(catRet))
}





#' Draw a sample for a virtual patient by sampling from the reads in her readstore. 
#' 
#' To get mixed reads from a tumor and normal genome (e.g. for a tumor probe with not 100\% tumor load) this functions samples from the read store.
#' The sample sequences are written in a subdirectoy under 'sample/'.
#' 
#' The read IDs within a sample are generally \strong{not} unique.
#' 
#' @author MK, 2014-02
#' @param patId  A single patient ID. 
#' @param baseDir  Base directory where all the relevant simulation files are stored.
#' @param myNGS.prop  This ngs.prop object controls which read store to access, i.e. the sequencing setting (e.g. read length and paired end details) for the reads
#' @param mySample.prop The sample.prop object controls how the sample settings (e.g. tumor load and coverage) are set.
#' @return \code{TRUE} if drawing of tumor and normal sample succeeded, \code{FALSE} otherwise.
#' @export
drawSamples <- function(patId, baseDir=BASEDIR, myNGS.prop=ngs.prop(), mySample.prop=sample.prop()){
  
  #paste0(baseDir, "readstore/", toString(myNGS.prop))
  STORE_BASEDIR <- getVirtualPatientPath(basedir = baseDir, what = "readstore", .ngs.prop = myNGS.prop)
  if (! file.exists(STORE_BASEDIR)){
    logging::logwarn("Could not find the readstore directory %s.", STORE_BASEDIR)
    return(FALSE)
  }
  
  if (length(patId) == 0L){
    logging::logwarn("No valid patient Id given.")
    return(FALSE)
  }
  
  if (length(patId) > 1L){
    patId <- patId[1L]
    logging::logwarn("Using only first of given patient IDs, namely %s", patId)
  }
  
  patStr <- getPatStr(patId)
  readsN1 <- paste0(STORE_BASEDIR, patStr, "_N_1.fq")
  readsN2 <- paste0(STORE_BASEDIR, patStr, "_N_2.fq")
  readsT1 <- paste0(STORE_BASEDIR, patStr, "_T_1.fq")
  readsT2 <- paste0(STORE_BASEDIR, patStr, "_T_2.fq")
  
  
  if ( ! all(checkFile(c(readsN1, readsN2, readsT1, readsT2))) ){
    logging::logwarn("Could not find the FASTQ-files in the READSTORE for NGS-prop %s and patient %s.", toString(myNGS.prop), patId)
    return(invisible(FALSE))
  }
  
  # get total number of reads: expects a single filename
  nbrReadsN <- lineCount(readsN1) / 4L
  nbrReadsT <- lineCount(readsT2) / 4L
  
  if (nbrReadsN == 0L || nbrReadsT == 0L){
    logging::logerror("No reads found in read store for patient %s (file %s).", patId, readsN1)
    return(invisible(NULL))
  } else logging::logdebug("Patient %d: Found %d reads in [N] readstore.", patId, nbrReadsN)
  
  # mkuhn, 2016-04-07: check that we can serve the requested coverage from the readstore
  stopifnot( mySample.prop[["cov"]] <= POOL_COV )
  
  # mkuhn, 20140623: add some random fluctuation to requested coverage (with arbitrary ad-hoc sd-value)
  covInSample <- pmin(POOL_COV, pmax(3L, mySample.prop[["cov"]] + rnorm(2L, mean = 0L, sd = 2.567)))
  
  desiredReadsN <- ceiling(nbrReadsN * covInSample[1L] / POOL_COV)
  desiredReadsT <- ceiling(nbrReadsT * covInSample[2L] / POOL_COV)
  
  logging::logdebug("Patient %d: Request %d reads for [N] sample.", patId, desiredReadsN)
  
  outputDir <- getVirtualPatientPath(basedir = baseDir, what = "sample", .ngs.prop = myNGS.prop, .sample.prop = mySample.prop, create = TRUE)
  
  subSampleFq_pl <- paste0(system.file("exec", package="svmod"), "/SubSample_Fastq.pl")
  
  if (! checkFile(subSampleFq_pl) ){
    logging::logwarn("Failed to find PERL script %s for sub-sampling.", subSampleFq_pl)
    return(invisible(FALSE))
  }
  
  
  # NORMAL reads ----
  
  logging::loginfo("Start with sub-sampling from the read pool for normal (N) probe of patient %s.", patId)
  subsampleN <- system(paste(subSampleFq_pl, "--desired ", desiredReadsN," -1 ",readsN1, " -2 ", readsN2, " --output ", outputDir), intern=FALSE)
  
  if (subsampleN == 0L){
    logging::logdebug("Sub-Sampling for NORMAL probe for paitent %s succeeded.", patId)
  } else {
    logging::logwarn("Sub-Sampling failed for NORMAL probe for patient %s.", patId)
    return(invisible(FALSE))
  }
  
  
  
  # TUMOR reads ----
  
  logging::loginfo("Start with sub-sampling from the read pool for tumor (T) probe of patient %s.", patId)
  TL <- mySample.prop[["tumorLoad"]]
  if ( TL >= 100L ){
    logging::logdebug("100% tumor load at tumor probe of patient %d.", patId)
    subsampleT <- system(paste(subSampleFq_pl, "--desired ", desiredReadsT," -1 ",readsT1, " -2 ", readsT2, " --output ", outputDir), intern=FALSE)
    appRet <- 0L # appending status OK as no temporary file appending necessary
  } else {
    
    # TLprop as proportion
    TLprop <- round(TL / 100L, 4L)
    stopifnot( TLprop > 0L, TLprop <= 1L )
    
    desiredReadsT_T <- ceiling(TLprop*desiredReadsT)
    desiredReadsT_N <- desiredReadsT - desiredReadsT_T
    
    subsampleT_T <- system(paste(subSampleFq_pl, "--desired ", desiredReadsT_T," -1 ",readsT1, " -2 ", readsT2, " --output ", outputDir, " --name tmpT"), intern=FALSE)
    subsampleT_N <- system(paste(subSampleFq_pl, "--desired ", desiredReadsT_N," -1 ",readsN1, " -2 ", readsN2, " --output ", outputDir, " --name tmpN"), intern=FALSE)
    # combining both return codes
    subsampleT <- abs(subsampleT_T) + abs(subsampleT_N)
    
    # 1st read
    appRet1 <- appendFiles(f1=paste0(outputDir, sub(".fq", ".subfqtmpN", basename(readsN1))), f2=paste0(outputDir, sub(".fq", ".subfqtmpT", basename(readsT1))),
                           f3=paste0(outputDir, sub(".fq", ".subfq", basename(readsT1))) )
    # 2nd read
    appRet2 <- appendFiles(f1=paste0(outputDir, sub(".fq", ".subfqtmpN", basename(readsN2))), f2=paste0(outputDir, sub(".fq", ".subfqtmpT", basename(readsT2))),
                           f3=paste0(outputDir, sub(".fq", ".subfq", basename(readsT2))) )
    # combining both return codes from concat
    appRet <- abs(appRet1) + abs(appRet2)
  }
  
  if ( subsampleT == 0L && appRet == 0L ){
    logging::logdebug("Sub-Sampling for TUMOR probe succeeded.")
  } else {
    logging::logwarn("Sub-Sampling failed somewhere for TUMOR probe of patient %s.", patId)
    return(invisible(FALSE))
  }#fi
  
  return(invisible(TRUE))
}






#' This function does the mapping for samples of virtual patients.
#' 
#' As specification, it needs to know a certain sample-directory which is defined by its sample.prop and its ngs.prop.
#' It uses BWA MEM mapper (via a wrapper function) and creates sorted and indexed BAM files. 
#' For real sequencing data, the mapping is usually already done: The full BAM file serves as input. 
#' 
#' @author mkuhn, 201402
#' 
#' @param myNGS.prop the NGS properties (i.e. PE insert size)
#' @param mySample.prop the sample properties (i.e. tumor/normal ratio in tumor sample)
#' @param patId optional vector of numeric patientIds whoses reads should be mapped. NULL [default] means to take all patients from specified sample directory
#' @param baseDir base directory where all simulated (virtual) data of patients, samples and mapped reads are stored
#' @param nbrThreads number of threads for mapping
#' @param genomeRef path of reference file,
#' @export
doMapping <- function(myNGS.prop=ngs.prop(), mySample.prop=sample.prop(), patId=NULL, nbrThreads=1L,
                      baseDir=BASEDIR, genomeRef=REF_GEN){
  
  
  #GENOM_REF <- system.file("inst", "extdata", "refGen", "hg19_chr5.fa", package="svmod")
  stopifnot( checkFile(genomeRef) )
  stopifnot( is.numeric(nbrThreads) && nbrThreads > 0 )
  
  
  # sample input directory
  inputDir <- getVirtualPatientPath(basedir = baseDir, what = "sample", .sample.prop = mySample.prop, .ngs.prop = myNGS.prop) 
  
  if (! file.exists(inputDir)){
    if ( isTRUE(dir.create(inputDir, showWarnings = FALSE, recursive = TRUE)) )
      loginfo("Created directory %s to hold sequences of sample!", inputDir)
    else {
      logwarn("Directory containing sequenced samples %s already exists or does not exist and could not be created!", dQuote(inputDir))
      return(invisible(FALSE))
    }
  }#fi
  
  
  # all sampled fastq-files
  ff <- list.files(path=inputDir, pattern=".subfq$", full.names = TRUE)
  
  if (length(ff) == 0){
    logging::logwarn("No sequence files from samples found in directory %s.", dQuote(inputDir))
    return(FALSE)
  }
  
  # patients involved
  ffPats <- unique(substr(basename(ff), 1, nchar(getPatStr(1))))
  logging::loginfo("Found sequenced samples to patients [%s].", paste(ffPats, collapse=", "))
  
  if ( ! is.null(patId) && is.numeric(patId) && length(patId) >= 1 ){
    ffPats <- intersect(ffPats, getPatStr(patId))
    logging::loginfo("Restrict mapping to patients [%s].", paste(ffPats, collapse=", "))
  }
  
  # Create output directory
  outputDir <- getVirtualPatientPath(basedir = baseDir, "mapping", .sample.prop = mySample.prop, .ngs.prop = myNGS.prop, create=TRUE) 
  
  
  for (pat in ffPats){
    logdebug("Start mapping of sequenced samples of patient %s.", pat)
    patN.fqs <- grep(paste0(pat, "[_]N[_]"), ff, value=TRUE)
    patT.fqs <- grep(paste0(pat, "[_]T[_]"), ff, value=TRUE)
    
    if (length(patN.fqs) != 2 || length(patT.fqs) != 2){
      logging::logwarn("No two sequencing files (paired sequencing) found for sample (T or N) of patient %s.", pat)
      next
    }
    
    # mkuhn, 2015-03-24: using mapping wrapper function
    mapWrapper(inputFile1 = patN.fqs, outputDir = outputDir, genomeRef = genomeRef, nbrThreads = nbrThreads)
    mapWrapper(inputFile1 = patT.fqs, outputDir = outputDir, genomeRef = genomeRef, nbrThreads = nbrThreads)
    
  }#rof pat
  
  logging::loginfo("Finished mapping of files in directory %s.", dQuote(inputDir))
}


#' Returns path to a directory with relevant simulated patient data.
#' 
#' It returns an absolute path.
#' @param create flag to indicate if the directory should be created. Defaults to \code{FALSE}.
#' @return Path to requested directory.
#' @export
getVirtualPatientPath <- function(basedir=BASEDIR, what=c("performance", "feature", "mapping", "sample", "readstore", "patient"),
                                  .sample.prop=sample.prop(), .ngs.prop=ngs.prop(), create=FALSE){
  what   <- match.arg(what)
  
  pathStr <- file.path(basedir,
                    switch (EXPR=what,
                            performance = file.path("performance", toString(.sample.prop), toString(.ngs.prop), "/"),
                            feature = file.path("feature", toString(.sample.prop), toString(.ngs.prop), "/"),
                            mapping = paste0("mapping/", toString(.sample.prop), "/", toString(.ngs.prop), "/"),
                            sample = paste0("sample/", toString(.sample.prop), "/", toString(.ngs.prop), "/"),
                            readstore = paste0("readstore/", toString(.ngs.prop), "/"),
                            patient = "patient/"))
  
  if ( isTRUE(create) && ! isTRUE(file.exists(pathStr)) ){
    if (! isTRUE(dir.create(pathStr, showWarnings = FALSE, recursive = TRUE)))
      logging::logerror("Failed to create directory %s!", pathStr) else
        logging::loginfo("Created directory")
  }
  pathStr
}


#' Gets all available patient IDs from the 'patData.rds' file.
#' 
#' @param basedir the base directory where to start looking.
#' @export
getPatientIds <- function(basedir=BASEDIR){
  patDataRDSFile <- paste0(getVirtualPatientPath(basedir=basedir, what="patient"), "patData.rds")
  
  if (! file.exists(patDataRDSFile)){
    logwarn("No 'patData.rds' file found.")
    return(invisible(NULL))
  }
  
  patData <- readRDS(patDataRDSFile)
  return(unique(patData$patId))
}


#' Clean up data from virtual patients
#' 
#' This function is a convenience function to clean up stale data of simulated patients.
#' Deleting a patient data means also to delete everything downstream (NGS- and sample setting do *not* matter).
#' Deleting a readstore means also to delete every patient related data concerning his mapping and sampling with same NGS-setting. (sample-setting does *not* matter).
#' @param patId numeric patient ID
#' @param scope what data is to be deleted, it includes everything downstream. 'patient' will delete everything from the given patient.
#' @param myNGS.prop defines the NGS-setting for which data should be removed.
#' @param mySample.prop defines the sample settings for which data should be removed
#' @return Cumulative return code of the different calls to \code{\link{unlink}}. Zero means no error.
#' @export
cleanVirtualPatientData <- function(patId, scope=c("patient", "readstore" , "sample", "mapping", "feature", "performance"),
                                    basedir=BASEDIR, myNGS.prop=ngs.prop(), mySample.prop=sample.prop()){
  scope <- match.arg(scope)
  patId_str <- paste(patId, collapse = "-")
  loginfo("About to clean %s-data for patient %s.", scope, patId_str)
  
  
  # PERFORMANCE
  cleanPerformance <- unlink(paste0(getVirtualPatientPath(basedir=basedir, what = "performance", .sample.prop=mySample.prop, .ngs.prop=myNGS.prop), getPatStr(patId), "*rds"), force=TRUE)
  if (cleanPerformance == 0L) logdebug("Cleaned Performance/-data for %s.", patId_str) else logwarn("Failure when cleaning Performance/-files for patient %s.", patId_str)
  if (scope == 'performance') return(cleanPerformance)
  
  
  # FEATURE
  cleanFeature <- unlink(paste0(getVirtualPatientPath(basedir=basedir, what = "feature", .sample.prop=mySample.prop, .ngs.prop=myNGS.prop), getPatStr(patId), "*rds"), force=TRUE)
  if (cleanFeature == 0L) logdebug("Cleaned feature/-data for %s.", patId_str) else logwarn("Failure when cleaning feature/-files for patient %s.", patId_str)
  if (scope == 'feature') return(cleanPerformance + cleanFeature)
  
  # MAPPING
  cleanMapping <- unlink(paste0(getVirtualPatientPath(basedir=basedir, what = "mapping", .sample.prop=mySample.prop, .ngs.prop=myNGS.prop), getPatStr(patId), "*ba?"), force=TRUE)
  if (cleanMapping == 0L) logdebug("Cleaned mapping/-data for %s.", patId_str) else logwarn("Failure when cleaning mapping/-files for patient %s.", patId_str)
  if (scope == "mapping") return(cleanPerformance + cleanFeature + cleanMapping)
  
  # SAMPLE
  cleanSample <- unlink(paste0(getVirtualPatientPath(basedir=basedir, what="sample", .sample.prop=mySample.prop, .ngs.prop=myNGS.prop), getPatStr(patId), "*subfq"), force=TRUE)
  if (cleanSample == 0L) logdebug("Cleaned sample/-data for %s.", patId_str) else logwarn("Failure when cleaning sample/-files for patient %s.", patId_str)
  if (scope == "sample") return(cleanPerformance + cleanFeature + cleanMapping + cleanSample)
  
  # SAMPLE and MAPPING are specific to sample- and NGS-properties.
  # READSTORE does only depend on ngs-property, PATIENT stands above all.
  stopifnot( scope %in% c('patient', 'readstore'))
  
  # list of all files related to patient (for downstream deletions)
  allPatFiles <- list.files(path=basedir, pattern = paste0(getPatStr(patId), collapse = "|"), recursive = TRUE)
  logdebug("All remaining files related to patient %s in %s are:\n%s", patId_str, basedir, paste(allPatFiles, collapse= " :: "))
  
    
  # READSTORE
  if (scope == 'readstore'){
    
    #cleanReadstore <- unlink(paste0(getVirtualPatientPath(basedir=basedir, what="readstore", .ngs.prop=myNGS.prop), getPatStr(patId), "*fq"), force=TRUE)  #, .sample.prop=mySample.prop ## does not matter here
    #if (cleanReadstore == 0) logdebug("Cleaned readstore/-data for %s.", paste(patId, collapse = "-")) else logwarn("Failure when cleaning readstore/-files for patient %s.", paste(patId, collapse = "-"))
    
    # mkuhn, 2015-04-17: when readstore is deleted for a given NGS-setting then also delete all downstreams files in sample/ and mapping/ with this same NGS-setting independently of sample settings
    allPatFiles.readstore <- grep(pattern=paste0("^readstore|mapping|sample.+", toString(myNGS.prop)), allPatFiles, value = TRUE)
    loginfo("Delete the following files downstream of readstore: %s", paste(allPatFiles.readstore, collapse=" :: "))
    cleanReadstore_downstream <- unlink(file.path(basedir, allPatFiles.readstore), force=TRUE)
    if (scope == "readstore") return(cleanMapping + cleanSample + cleanReadstore_downstream)  ## + cleanReadstore
  }
  
  stopifnot( scope == 'patient' )
  
  # cleanPatient <- unlink(paste0(getVirtualPatientPath(basedir=basedir, what="patient"), getPatStr(patId), "*fasta"), force=TRUE)
  #loginfo("Delete the following files downstream of patient: %s", paste(allPatFiles.patient, collapse=" :: "))
  cleanPatient_downstream <- unlink(file.path(basedir, allPatFiles), force=TRUE)
  if (cleanPatient_downstream == 0L) logdebug("Cleaned patient/-data for %s.", patId_str) else logwarn("Some trouble when cleaning patient/-files for patient %s.", patId_str)
  
  # patient RDS data
  cleanPatient_rds <- 0L
  patDataRDSFile <- file.path(getVirtualPatientPath(basedir=basedir, what="patient"), "patData.rds")
  if (! checkFile(patDataRDSFile)){
    logwarn("Failed to open the 'patData.rds'-file.")
    cleanPatient_rds <- 17L # indicate error if patData.rds not there
  } else {
    loginfo("Drop patients %s from 'patData.rds'-file.", patId_str)
    patData <- readRDS(patDataRDSFile)
    patData.n <- NROW(patData)
    patData <- droplevels(patData[ ! patData$patId %in% patId, ])
    if ( NROW(patData) > 0L ){
      if (NROW(patData) < patData.n) saveRDS(patData, file = patDataRDSFile) 
    } else unlink(patDataRDSFile, force=TRUE)
    loginfo("Updated 'patData.rds' file.")
  }  
  
  return(invisible(cleanPerformance + cleanFeature + cleanMapping + cleanSample + cleanPatient_rds + cleanPatient_downstream))
}





#' Handle virtual patients identified by their IDs. 
#' 
#' This function combines some steps to generate sequencing data of virtual patients.
#' I.e. patient creation and its readstore and sampling. Mapping is delegated to another method (see \code{doMapping}).
#' 
#' @param patIds integer vector of virtual patient IDs
#' @param preserve logical flag: should existing virtual data be preserved and reused?
#' @param randomSeed integer. If greater 0 than set random seed
#' @param do.mapping logical flag: should all new sampled data be mapped?
#' @param ... additional arguments (\code{nbrLoci, nbrDSV, minMargin, SVlength_mean}) passed on to function \code{createPatient}
#' @return dataframe of all patients that were created. #ZZZ should that be stored as GRanges object or GRangesList object (one GRanges per patient) ?!
#' @export
processVirtualPatients <- function(patIds, basedir=BASEDIR, preserve=TRUE, do.mapping=TRUE, nbrCores=1L, randomSeed=0L,
                                   sample.prop=sample.prop(), ngs.prop=ngs.prop(), ...){
  
  if (length(patIds) < 1L || ! is.numeric(patIds)) return(NULL)
  
  if (is.numeric(randomSeed) && randomSeed > 0L){
    set.seed(ceiling(randomSeed))
    logging::loginfo("Set random seed in processVirtualPatients-call to %s.", ceiling(randomSeed))
  }
  
  # which patIds need mapping (because they got new sample fq-files)
  mapping.patIds <- numeric(length(patIds))
  
  
  # mkuhn, 20140624: problems with parallel execution in the executed Perl script
  # Perl script uses tmp-file 'temp.txt'. This is *not* parallel-/thread-safe.
  patData <- foreach (i=seq_along(patIds),  .combine = 'rbind') %do% {
    myPatId <- patIds[i]
    
    myPatDataEntries <- NULL
    if ( ! checkFile(f=file.path(getVirtualPatientPath(basedir = basedir, "patient", create = TRUE), paste0(getPatStr(myPatId), "Loci_T.fasta"))) || ! isTRUE(preserve) ){
      logging::loginfo("Creating patient %s with %s loci and %s DSVs", myPatId, if (exists("nbrLoci")) nbrLoci else "?", if (exists("nbrDSV")) nbrDSV else "?")
      myPatDataEntries <- createPatient(baseDir = myBaseDir, patId=myPatId, ...) #minMargin, nbrLoci=nbrLoci, nbrDSV=nbrDSV, SVlength_mean = SVLength_mean, SVlength_sd=SVLength_sd)
    }
    
    if ( ! checkFile(f=file.path(getVirtualPatientPath(basedir = basedir, "readstore", .ngs.prop = ngs.prop, create = TRUE), paste0(getPatStr(myPatId), "_T_1.fq"))) || ! isTRUE(preserve) ) {
      ## Generate read pool of reads from normal and tumor sample of patient  
      logging::loginfo("Creating read pool for patient %s.", myPatId)
      generatePatientReadStore(patId=myPatId, myNGS.prop = ngs.prop)
    }
    
    
    if ( ! checkFile(f=file.path(getVirtualPatientPath(basedir = basedir, "sample", .ngs.prop = ngs.prop, .sample.prop=sample.prop, create = TRUE), paste0(getPatStr(myPatId), "_T_1.subfq"))) || ! isTRUE(preserve) ) {
      ## generate the sequenced samples (Tumor and Normal)
      logging::loginfo("Draw samples for patient %s.", myPatId)
      drawSampleRes <- drawSamples(patId=myPatId, myNGS.prop=ngs.prop, mySample.prop=sample.prop)
      if ( isTRUE(drawSampleRes) ) mapping.patIds[i] <- myPatId else logging::logwarn("Failed to draw sample for patient %s", myPatId)
    }
    myPatDataEntries
  } #hcaerof
  
  
  mapping.patIds <- mapping.patIds[which(mapping.patIds > 0L)]
  
  if ( isTRUE(do.mapping) && length(mapping.patIds) > 0L ){
    # do mapping of the new patients in the current setting
    logging::loginfo("Start with mapping of %s new patients.", length(mapping.patIds))
    doMapping(patId=mapping.patIds, nbrThreads=max(1L, nbrCores), myNGS.prop = ngs.prop, mySample.prop = sample.prop) #seq_len(NBR_PATS))
  }#fi
  
  return(patData)
}


#' Extracts features from the simulated data of a virtual patient
#' 
#' The function uses the known true situation of the simulated SV, i.e. for simulated SVs the features are extracted around it.
#' For loci with no somatic SV, mapping characteristcs are randomly sampled in a similar fashion like for the somatic SV.
#' 
#' @param patId numeric patient ID
#' @param patData dataframe with rows per simulated locus (including status info) for all patients
#' @param sample.prop sample property of simulated data
#' @param ngs.prop sequencing property of simulated data
#' @param margin.bp integer value (in nucleotides) what window to the left and right of the SV is taken for feature extraction.
#' @return a dataframe of extracted features, including simulation context information (i.e. sample and sequencing information and status of locus)
#' @export
getVirtualPatientFeatures <- function(patId, patData, baseDir, sample.prop, ngs.prop, margin.bp=25L){
  
  mappingDir <- getVirtualPatientPath(baseDir, what = "mapping", .sample.prop=sample.prop, .ngs.prop=ngs.prop)
  pat.n.bamFile <- list.files(path=mappingDir, pattern=paste0("^", getPatStr(patId),"[_]N.bam$"), full.names = TRUE)
  pat.t.bamFile <- list.files(path=mappingDir, pattern=paste0("^", getPatStr(patId),"[_]T.bam$"), full.names = TRUE)
  
  stopifnot( length(pat.n.bamFile) == 1L && checkFile(pat.n.bamFile) )
  stopifnot( length(pat.t.bamFile) == 1L && checkFile(pat.t.bamFile) )
  
  
  # use simulation information
  thePatData <- patData[patData$patId == patId,]
  stopifnot( length(table(thePatData$chrom)) == 1L )
  patChrom <- as.character(thePatData$chrom[1L])

  myPat.nbrNoSVRegions <- sum(thePatData$SVlength == 0L)
  myPatSVlength.mean <- mean(thePatData$SVlength[which(thePatData$SVlength > 0L)])
  myPatSVlength.sd <- sd(thePatData$SVlength[which(thePatData$SVlength > 0L)])
  logging::loginfo("Data of sim-patient %s with %d target regions. %d SVs. %d non-SV regions.",
          patId, NROW(thePatData), sum(thePatData$SVlength>0),  myPat.nbrNoSVRegions)
  
  #loginfo("Mean length of SVs of this patient is %f. (%s)", myPatSVlength.mean, paste(thePatData$SVlength[which(thePatData$SVlength > 0)]))
  
  # Regions without a somatic SV:
  #+set also artificial SV characteristics for feature extraction later (aka sim-HACK?!)
  # mkuhn, 2016-02-04: I use SVlength as it fits into the target regions, at least SV_MIN_LENGTH
  ##thePatData[thePatData$SVlength == 0L, "SVlength"] <- ceiling(0.1 + abs(rnorm(myPat.nbrNoSVRegions, mean=myPatSVlength.mean, sd=myPatSVlength.sd)))
  
  for (i in seq_len(NROW(thePatData))){
    if ( thePatData[i, "SVstart"] > 0L ) next
    stopifnot( thePatData[i, "SVstart"] == 0L )
    seqLen <- thePatData[i, "targetEndPos"] - thePatData[i, "targetStartPos"] + 1L #mkuhn, 2016-02-04 +1 for included end-coordinate?!
    stopifnot( seqLen > 0L, seqLen - 2L * margin.bp > SV_MIN_LENGTH )
    
    thePatData[i, "SVlength"] <- min(SV_MAX_LENGTH, SV_MIN_LENGTH + sample(seqLen - 2L * margin.bp - SV_MIN_LENGTH, size = 1L))
    thePatData[i, "SVstart"]  <- margin.bp + sample(seqLen - 2L * margin.bp - thePatData[i, "SVlength"], size = 1L) #like for SV regions
  }# rof #sim-HACK
  
  # mkuhn, 2016-02-04: every SV is filled in
  stopifnot( all(thePatData$SVlength > 0L), all(thePatData$SVstart > 0L) )

  #   # insertions have only single breakpoint with respect to clipped base peaks
  #   clippingUpperBound <- ifelse(thePatData$SVtype == 'insertion', thePatData$targetStartPos + thePatData$SVstart + margin.bp,
  #                        thePatData$targetStartPos + thePatData$SVstart + thePatData$SVlength + margin.bp)
  thePatData$selectedStartPos <- pmax(thePatData$targetStartPos, thePatData$targetStartPos + thePatData$SVstart - margin.bp)
  thePatData$selectedEndPos <-  pmin(thePatData$targetEndPos, thePatData$targetStartPos + thePatData$SVstart + thePatData$SVlength + margin.bp)
  
  # filter those selected regions that have too few Q1-read
  thePatData$Q1filterOut <- 0L
  
  for (i in 1:NROW(thePatData)){
    myStartPos <- thePatData[i, "selectedStartPos"]
    myEndPos <- thePatData[i, "selectedEndPos"]
    
    SELECTED_REGION <- paste0(patChrom, ":", myStartPos, "-", myEndPos)
    rc.q1.N <- as.numeric(system(paste(SAMTOOLS_EXE, "view -c -q1 ", pat.n.bamFile, SELECTED_REGION), intern=TRUE))
    rc.q1.T <- as.numeric(system(paste(SAMTOOLS_EXE, "view -c -q1 ", pat.t.bamFile, SELECTED_REGION), intern=TRUE))
    # Should I use these functions to get already clipped read information?
    #     baseN <- getClippedBasePositions(bamFile = pat.n.bamFile,
    #                             targetStartPos = myStartPos, targetEndPos = myEndPos)
    #     baseT <- getClippedBasePositions(bamFile = pat.t.bamFile,
    #                                      targetStartPos = myStartPos, targetEndPos = myEndPos)
    # Used here before: sum(baseN$cov)
    if (rc.q1.N < 5L * (myEndPos - myStartPos) / READL) thePatData$Q1filterOut[i] <- 1L
    if (rc.q1.T < 3L * (myEndPos - myStartPos) / READL) thePatData$Q1filterOut[i] <- 2L
  }#rof
  logdebug("For patient %d there are %d entries of which %d entries pass Q1-filter.", patId, NROW(thePatData), length(which(thePatData$Q1filterOut==0L)))
  thePatData <- dplyr::filter_(thePatData, ~Q1filterOut == 0L)
  
  
  ## collect clipping info per target region
  clipInfos <- vector(mode = "list", length = NROW(thePatData))
  for (i in 1:NROW(thePatData)){
    
    # windowSize=31L, minThresholdAbs=2L, minThreshold=0.05, thresholdQuantile=0.95){ 
    clipInfos[[i]] <- assessClippedBasesInRegion(bamFile = pat.t.bamFile, chrom = patChrom, #thePatData[i, "chrom"],
                                         .targetStartPos = thePatData[i, "targetStartPos"], .targetEndPos = thePatData[i, "targetEndPos"])[["peakInfo"]]
  }#rof
  
  # target regions come directly from the simulation where we used a target enrichment file and sampled some target regions.
  pat.targetRegions <- IRanges::IRanges(start=thePatData$targetStartPos, end = thePatData$targetEndPos)
  pat.regions <- IRanges::IRanges(start=thePatData$selectedStartPos, end = thePatData$selectedEndPos)
  
  
  # get mapping info for patient
  patMappingInfo <- getPatientMappingInfo(bamFile_WT = pat.n.bamFile, bamFile_MUT = pat.t.bamFile)
  
  loginfo("Found %d candidate regions for patient %d.", NROW(pat.regions), patId)
  # parallelism happened here in the regions, but gave corruption problems (maybe because of external system-calls, like samtools??)
  ffm <- getFeaturesFromMappings(patId = patId, bamFile_WT = pat.n.bamFile, bamFile_MUT = pat.t.bamFile, patMappingInfo = patMappingInfo, 
                                 chrom = patChrom, regions = pat.regions, targetRegions = pat.targetRegions, clipInfos=clipInfos) #, statusInfo=thePatData$status)
  
  # mkuhn, 20150425: add context-information of patient simulation to feature data. We should not loose any feature entry
  nbrFeatures <- NROW(ffm)
  ffm <- dplyr::inner_join(x=dplyr::select_(thePatData,  ~patId, ~chrom, ~targetStartPos, ~targetEndPos, ~status, ~SVtype),
                           y=ffm, by=c(patId = "patId", chrom = "chrom", targetStartPos = "targetStartPos", targetEndPos = "targetEndPos"))
  stopifnot(nbrFeatures == NROW(ffm))
  # fixing the status labels: I could do that once and for all at the basis at patData!? ZZZ 
  stopifnot( is.character(ffm$status), identical(unique(ffm$status), c("DSV", "germline")) )
  stopifnot( is.character(ffm$SVtype), "no" %in% unique(ffm$SVtype) )
             
  ffm$status <- factor(ffm$status, levels = c("germline", "DSV"), labels = STATUS_LEVELS)
  ffm$SVtype <- relevel(factor(ffm$SVtype), ref = "no")
  
  # mkuhn, 20150422: patient and locus info first, then sample and NGS-scenario values from the simulation
  ffm1 <- ffm %>%
    dplyr::select_(~patId, ~chrom, ~startPos, ~endPos, ~targetStartPos, ~targetEndPos, ~status, ~SVtype)
  ffm2 <- cbind(ffm1,
                margin=rep(margin.bp, NROW(ffm)),
                cov=rep(sample.prop[["cov"]], NROW(ffm)),
                tumorLoad=rep(sample.prop[["tumorLoad"]], NROW(ffm)),
                read.length=rep(ngs.prop[["read.length"]], NROW(ffm)),
                pe.ins.mean=rep(ngs.prop[["pe.ins.mean"]], NROW(ffm)),
                pe.ins.sd=rep(ngs.prop[["pe.ins.sd"]], NROW(ffm)),
                ffm[, !names(ffm) %in% names(ffm1)])
  
  return(ffm2)
}


