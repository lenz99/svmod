# mkuhn, 20130828
# init script for svmods-project
# it sets some default values and logging stuff
#####

### set up the logging
logging::basicConfig()
logging::setLevel(level='INFO')  #set level of root logger
logging::setLevel('INFO', getHandler('basic.stdout')) # set level of handler


#' @export
PERF_PREFIX <- "perf"

# fixed chromosome to test the package
CHROM <- "chr5"

#' read length used in the simulations
READL <- 100L

#' minimum legnth of SVs that we are looking at
SV_MIN_LENGTH <- 50L
SV_MAX_LENGTH <- 300L

#' SIM_SV_MEAN_LENGTH
SIM_SV_LENGTH_MEAN <- 100L
SIM_SV_LENGTH_SD <- 9L
stopifnot( SIM_SV_LENGTH_MEAN > SV_MIN_LENGTH + 3L * SIM_SV_LENGTH_SD )

#' margin in bp to widen the target regions
EXTRA_MARGIN_TARGET_REGIONS <- 3L * READL
#' minimal width of covered regions to be found
MIN_COV_REGION_WIDTH <- 3L * READL

#' minimal requested coverage for WT to be considered a "covered region"
MINCOV_WT <- 15L

#' maximum extended width of a read pair to count for the affected region to simulate reads in SV-injection
MAX_PAIR_WIDTH <- 2L * READL + 1000L #ad-hoc upper bound for max insert size betw the two read ends

#' sequencing system (for read simulation ART)
SEQ_SYSTEM <- "HS20"

#' Supported types of structural variations.
#' @export
SV_TYPES <- c("insertion", "deletion", "duplication", "inversion")



# coverage from the read store
POOL_COV <- 150L


# BAM flags
FLAG_PAIRED <- 1L
FLAG_PROPER_PAIR <- 2L
FLAG_READ_UNMAPPED <- 4L
FLAG_MATE_UNMAPPED <- 8L
FLAG_FIRST_IN_PAIR <- 64L
FLAG_SECOND_IN_PAIR <- 128L
FLAG_SECONDARY <- 256L  #old name in Rsamtools: isNotPrimaryRead 
FLAG_SUPPL <- 2048L
FLAG_READ_REV <- 16L
FLAG_MATE_REV <- 32L


#' Labels used for 2-class classification for outlier detection
#' @export
STATUS_LEVELS <- c("IN", "OUT")


FEAT_COLS0 <- c("patId", "chrom", "startPos", "endPos", "targetStartPos", "targetEndPos")
FEAT_COLS1 <- c("status", "SVtype", "SVlength", "SVmargin")
FEAT_COLS <- c(FEAT_COLS0, FEAT_COLS1)

LIN_BASE <- "~/seq1" #"/data/sequencing"

EXEC_BASEDIR_MAC <- "/Users/kuhnmat/.bin/"
EXEC_BASEDIR_LIN <- file.path(LIN_BASE, "tools/")

ON_MACKIE <- file.exists(EXEC_BASEDIR_MAC) && isTRUE(file.info(EXEC_BASEDIR_MAC)$isdir)

EXEC_BASEDIR <- if (ON_MACKIE) EXEC_BASEDIR_MAC else EXEC_BASEDIR_LIN

BASEDIR_EXT <- "/Volumes/External_Toshiba/SVdata/"
BASEDIR_MAC <- "/Users/kuhnmat/SVdata/"
BASEDIR_LIN <- file.path(LIN_BASE,"SVdata/")

#' Base directory where the data (of simulation) is stored.
#' 
#' Later on this should be more flexible so that svmod can be called easily with any file structure
#' @export
BASEDIR <- if ( file.exists(BASEDIR_EXT) && isTRUE(file.info(BASEDIR_EXT)$isdir) ) BASEDIR_EXT else
  if ( file.exists(BASEDIR_MAC) && isTRUE(file.info(BASEDIR_MAC)$isdir) ) BASEDIR_MAC else BASEDIR_LIN

stopifnot( file.exists(EXEC_BASEDIR) )
stopifnot( file.exists(BASEDIR) )

logdebug("Base directory is ::%s:: with executable base directory ::%s::", BASEDIR, EXEC_BASEDIR)

#' path to reference genome
REF_GEN <- file.path(BASEDIR, "refGen", "hg19_normal.fa")
stopifnot( checkFile(REF_GEN) )

#' @export
TMPDIR <- file.path(BASEDIR, "tmp")
if ( ! file.exists(TMPDIR) ){
  try(dir.create(TMPDIR))
}

stopifnot( file.exists(TMPDIR))

ART_EXE <- paste0(EXEC_BASEDIR, "art/art_illumina")
BWA_EXE <- paste0(EXEC_BASEDIR, "bwa/bwa")
#' @export
SAMTOOLS_EXE <- paste0(EXEC_BASEDIR, "samtools/samtools")
BEDTOOLS <- paste0(EXEC_BASEDIR, "bedtools/bedtools") # "bedtools2/bin/bedtools")

stopifnot( checkFile(SAMTOOLS_EXE), checkFile(BEDTOOLS), checkFile(BWA_EXE), checkFile(ART_EXE) )


REF_GEN <- file.path(BASEDIR, "refGen", "hg19_normal.fa")
stopifnot( checkFile(REF_GEN) )

#' ngs.prop constructor.
#' 
#' Properties of the sequencing preparation and sequencing run. It is about the read length and paired-end insert size distribution.
#' 
#' @param read.length the read length of the Illumina sequencing machine
#' @export
ngs.prop <- function(read.length=100L, pe.ins.mean=200L, pe.ins.sd=20L){
  z <- list(read.length=read.length, pe.ins.mean=pe.ins.mean, pe.ins.sd=pe.ins.sd)
  class(z) <- "ngs.prop"
  return(z)
}

#' @export
toString.ngs.prop <- function(x){ stri <- paste0("PEm",x$pe.ins.mean, "sd", x$pe.ins.sd, "_R",x$read.length); stri }
#' @export
print.ngs.prop <- function(x){ print(toString(x)) }
#' Parse string that encodes a NGS property.
#' @param s string representation of a NGS property or a directory path
#' @export
parseNGSProp <- function(s){
  r <- NULL
  
  ngsMatch <- regexec("PEm[[:digit:]]+sd[[:digit:]]+_R[[:digit:]]+", s)[[1]]
  if (ngsMatch > -1L){ #grepl(ngsPattern, s)){
    # NGS match
    s <- substr(s, start = ngsMatch, stop = ngsMatch + attr(ngsMatch, "match.length")-1L)
    stopifnot(grepl("^PEm", s))
    s <- substr(s, 4, nchar(s)) #PEm
    match.m <- regexec("^[[:digit:]]+", s)[[1]]
    if ( match.m == 1){
      pe.m <- strtoi(substr(s, 1, attr(match.m, "match.length")), base=10)
      s <- substr(s, attr(match.m, "match.length")+1+2, nchar(s)) #sd
      
      match.sd <- regexec("^[[:digit:]]+", s)[[1]]
      if ( match.sd == 1){
        pe.sd <- strtoi(substr(s, 1, attr(match.sd, "match.length")), base=10)
        s  <- substr(s, attr(match.sd, "match.length")+1+2, nchar(s)) #_R
        
        match.r <- regexec("^[[:digit:]]+$", s)[[1]]
        if (match.r == 1){
          r <- strtoi(s, base=10)
          
          r <- ngs.prop(read.length=r, pe.ins.mean=pe.m, pe.ins.sd=pe.sd)
        }
      }
    }
  }#fi ngsMatch
  return(r)
}

#' Use insert size distribution from NGS-property for BWA calling (at least for FR-oriented reads)
getISDistForBWA <- function(x){
  paste0("-I", x[["pe.ins.mean"]], ",",x[["pe.ins.sd"]])
}


#' Sample properties: tumor load in sample and coverage.
#' 
#' 
#' Coverage is listed here although it has a link to NGS setting. Still, it is \emph{not} treated as a NGS property
#' because the read-store in simulations is organized by NGS properties and coverage does not belong there.
#' Another case are the SV characteristics of the tumor. How frequent are the SVs of the tumor and what mean length do these SV have? But this can be seen in the patient data.
#' 
#' It is a S3-constructor. 
#' @export
sample.prop <- function(cov=60, tumorLoad=90){
  z <- list(cov=cov, tumorLoad=tumorLoad)
  class(z) <- "sample.prop"
  return(z)
}

#' @export
toString.sample.prop <- function(x){ stri <- paste0("cov", x$cov, "_TL",x$tumorLoad); stri }
#' @export
print.sample.prop <- function(x){ print(toString(x)) }

#' Parse string that encodes a sample property.
#' @param  s string representation of a sample property
#' @export
parseSampleProp <- function(s){
  r <- NULL
  
  sampleMatch <- regexec("cov[[:digit:]]+_TL[[:digit:]]+", s)[[1]]
  if (sampleMatch > -1L) {
    # sample match
    s <- substr(s, start = sampleMatch, stop = sampleMatch + attr(sampleMatch, "match.length")-1L)
    stopifnot(grepl("^cov", s))
    
    s <- substr(s, 4, nchar(s))
    match.cov <- regexec("^[[:digit:]]+", s)[[1]]
    if (match.cov == 1L){
      myCov <- strtoi(substr(s, 1, attr(match.cov, "match.length")), base=10L)
      s <- substr(s, attr(match.cov, "match.length")+1+3, nchar(s))
      match.tl <- regexec("^[[:digit:]]+$", s)[[1]]
      if (match.tl == 1L){
        myTL <- strtoi(s, base=10L)
      
        r <- sample.prop(cov=myCov, tumorLoad=myTL)
      }# fi match TL
    }# fi match cov
  }# fi sampleMatch
  return(r)
}





