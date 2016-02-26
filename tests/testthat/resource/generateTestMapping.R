# mkuhn, 2015-12-14
# script to generate NGS-reads and mapping for testing.
# 
# Input: FASTQ file "ill_reads_pe_is300_hg19_art.fastq" 
# ++ with simulated PE-reads with concordant and discordant read pairs, cf function startWithSimulatedReadPairs().
# Output: sorted indexed BAM file

library(GenomicRanges)
library(Rsamtools)
library(ShortRead)

library(svmod)

# for reproducible data
set.seed(1605)


TEST_DIR <- "~/projects/svmod/tests/testthat/resource/"
stopifnot( file.exists(TEST_DIR) )
#BASEDIR <- "~/SVdata"
REF_GEN <- file.path(BASEDIR, "refGen", "hg19_normal.fa")
TARGET_BED <- file.path(BASEDIR, "refGen", "TruSeq_exome_targeted_regions.hg19.bed")


# get handcrafted read-pairs
CHROM <- "chr5"
READL <- 100L
READL_LONG <- 170L

# # NPM1 two last targets
# 170834703       170834778       chr5:170834704-170834778:NPM1
# 170837530       170837887       chr5:170837531-170837887:NPM1
POS <- GRanges(seqnames = CHROM, ranges = IRanges(start=170837750L, end=170837750L))
SV_POS <- resize(POS, width = 51L, fix = "center")




QUAL_CHARS <- strsplit("?@ABCDEFGH", split="")[[1L]]
readQual0  <- ShortRead::FastqQuality(paste0(sample(QUAL_CHARS, size = READL, replace = TRUE), collapse="")) #, rep(NBR_FASTQ_ENTRIES))
readQual0L <- ShortRead::FastqQuality(paste0(sample(QUAL_CHARS, size = READL_LONG, replace = TRUE), collapse="")) #, rep(NBR_FASTQ_ENTRIES))


# insert size 300
grPE1 <- GRanges(seqnames = CHROM, ranges = IRanges(start = c(start(SV_POS)-125L, end(SV_POS)+ 25L, #SVgap
                                                              start(SV_POS) -25L, end(SV_POS)+125L, #SVread
                                                              start(SV_POS) -85L), #SVread 
                                                    width=READL))
grSE1 <- GRanges(seqnames = CHROM, ranges = IRanges(start = start(SV_POS)-floor(.85*READL_LONG), width = READL_LONG))
                 #configName = c(rep("SVgap", 2L), rep("SVread", 3L)))
## longer SE read to map?!

# pseudo gene at chr8?!
#grALT <- GRanges(seqnames = c(CHROM, "chr8"), ranges = IRanges(start = c(170837700, 62115916), width=READL))
#Biostrings::getSeq(open(Rsamtools::FaFile(REF_GEN)), param=grALT)

readSeqRaw <- c(Biostrings::getSeq(open(Rsamtools::FaFile(REF_GEN)), param=grPE1),
                Biostrings::getSeq(open(Rsamtools::FaFile(REF_GEN)), param=grSE1) )
#mcols(readSeqRaw) <- mcols(grPE1)


randReads <- DNAStringSet(Biostrings::DNAString(paste(sample(DNA_BASES, size = READL, replace = TRUE), collapse="")))
randReads <- Biostrings::append(randReads, DNAStringSet(Biostrings::DNAString(paste(sample(DNA_BASES, size = READL, replace = TRUE), collapse=""))))


# reads for interleaved FASTQ
readSeq <- c(
  readSeqRaw[1], Biostrings::reverseComplement(readSeqRaw[2]), #SV gap, PE normal
  Biostrings::reverseComplement(readSeqRaw[2]),readSeqRaw[1],  #SV gap, PE normal
  readSeqRaw[3], Biostrings::reverseComplement(readSeqRaw[4]), #SVread, PE normal
  readSeqRaw[1], randReads[1], #SV gap, PE randR2
  readSeqRaw[3], randReads[2], #SVread, PE randR2
  readSeqRaw[1],                    #SV gap, SE
  reverseComplement(readSeqRaw[2]), #SV gap, SE
  readSeqRaw[3],                    #SVread, SE
  readSeqRaw[5],                    #SVread, SE
  readSeqRaw[6],                    #SVread, SE (longer so it properly maps)
  readSeqRaw[1], readSeqRaw[2], #SV gap, PE discordant: +/+
  readSeqRaw[3], readSeqRaw[4]  #SVread, PE discordant: +/+
)


readNames <- c("PE1_SVgap_norm/1", "PE1_SVgap_norm/2", "PE2_SVgap_norm/1", "PE2_SVgap_norm/2", "PE3_SVread_norm/1", "PE3_SVread_norm/2",
               "PE4_SVgap_randR2/1", "PE4_SVgap_randR2/2", "PE5_SVread_randR2/1", "PE5_SVread_randR2/2",
               "SE1_SVgap/1", "SE2_SVgap/1", "SE3_SVread/1", "SE4_SVread/1", "SE5_SVread/1",
               "PE5_SVgap_disc/1", "PE5_SVgap_disc/2", "PE6_SVread_disc/1", "PE6_SVread_disc/2")

stopifnot( length(readNames) == length(readSeq) )

# building up the qualities
readQual <- readQual0
for (i in seq_along(readNames)){
  if (i == 1L) next
  readQual <- append(readQual, if (width(readSeq)[i] == READL) readQual0 else readQual0L)
}

srObj <- ShortRead::ShortReadQ(sread = readSeq,
                                quality = readQual,
                                id = Biostrings::BStringSet(readNames) )

targetFastqFile <- file.path(TEST_DIR, "test_SV0.fastq")
# starting point: group of paired end reads: to let BWA estimate insert size distribution
stopifnot( file.copy(from = file.path(TEST_DIR, "ill_reads_pe_is300_hg19_art.fastq"), to = targetFastqFile, overwrite = TRUE) )
# append my few reads on chr5
ShortRead::writeFastq(srObj, file = targetFastqFile, compress=FALSE, mode = "a") 

# map all reads
mapWrapper(inputFile1 = targetFastqFile, smartPairing = TRUE) #, bwa_optsExtra = "-I 300")





#' Starting point are some paired-end reads: FR but also FF orientation.
#' This is to have the mate-rescue working for BWA
#' This is for documentation how I arrived at the starting point of 600 PE-reads (half reversed by simply reversing the second read which gives a smaller insert size estimation)
startWithSimulatedPairedReads <- function(){
  
  # call to art to simulate paired reads: this was terminated after a while
  #art_illumina -na -ss HS20  -p -l 100 -m 300 -s 1 --rcount 10 -i ~/SVdata/refGen/hg19_normal.fa -o reads_dat
  # 
  # merge PE reads into a single file
  #seqtk mergepe reads_dat* > out.fastq
  #
  # renamed  out.fastq to ill_reads_pe_is300_hg19_art.fastq
  
  # adding FF-oriented read-pairs
  PE_ART_SIM_FILE <- "ill_reads_pe_is300_hg19_art.fastq"
  illR <- ShortRead::readFastq(dirPath = ".", pattern=PE_ART_SIM_FILE)
  isFirstRead <- 1:length(illR) %% 2 == 1L
  
  
  newReads <- c(sread(illR)[isFirstRead],
                Biostrings::reverseComplement(sread(illR)[!isFirstRead]))
  
  newNames <- Biostrings::BStringSet(paste0("REV",
                                            c(as.character(id(illR)[isFirstRead]), 
                                              as.character(id(illR)[!isFirstRead]))
  ))
  
  
  srObj <- ShortRead::ShortReadQ(sread = newReads[order(newNames)],
                                 quality = quality(illR),
                                 id = newNames[order(newNames)] )
  # add REV-PE reads to the simulated FASTQ-file  
  ShortRead::writeFastq(srObj, file = PE_ART_SIM_FILE, compress=FALSE, mode = "a") 
  
}
