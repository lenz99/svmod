library(svmod)

context("Inject SV")

baseDir <- BASEDIR

#ZZZ check more end range ?!
test_that("Test that 'cigar to pos'-conversion works.", {
  myTestChrom <- "chr5"
  expect_equal(cigarToPos("100M",         chrom=myTestChrom, pos=1000L), GenomicRanges::GRanges(seqnames = myTestChrom, IRanges::IRanges(1000L, 1099L)) )
  expect_equal(cigarToPos("4S96M",        chrom=myTestChrom, pos=1000L), GenomicRanges::GRanges(seqnames = myTestChrom, IRanges::IRanges(996L, 1095L)) )
  expect_equal(cigarToPos("4S95M1S",      chrom=myTestChrom, pos=1000L), GenomicRanges::GRanges(seqnames = myTestChrom, IRanges::IRanges(996L, 1095L)) )
  expect_equal(cigarToPos("4S29M4I13M4S", chrom=myTestChrom, pos=1000L), GenomicRanges::GRanges(seqnames = myTestChrom, IRanges::IRanges(996L, 1045L)) )
})

# test_that("Test that insertion-SV work", {
#   expect_equal(sqrt(49), 7L)
# })


# the reason to keep 2 SV-functions is that for tumor probe creation we also add SNV (and that is not part of SV-object)
test_that("Test that the SV-injection at tumor creation is the same as for via injection of SV-object.", {
  genomeRef <- file.path(baseDir, "refGen", "hg19_chr5.fa")
  
  myTargetReg <- GenomicRanges::GRanges(seqnames = "chr5", ranges = IRanges::IRanges(start=76707825, end = 76708160))
  myTargetSeq <- Rsamtools::getSeq(open(Rsamtools::FaFile(genomeRef)), param=myTargetReg)[[1L]]
  
  mySVLen <- 51L
  
  set.seed(123L)
  insMut <- injectMySVs(myTargetSeq, SVtype="insertion", SVstartPos = 11L, SVlength=mySVLen)
  set.seed(123L)
  insSV <- SV(type = "insertion", chrom = "chr5", posL = GenomicRanges::start(myTargetReg) + 11L - 1L, seqLen = mySVLen, refGen = genomeRef)
  
  expect_identical(toString(insMut[["SVseq"]]), toString(getSVSeq(insSV)))
  expect_identical(toString(insMut[["DNAseq"]]), toString(getContextSeq(insSV, gr = myTargetReg)))
  
  
  dupMut <- injectMySVs(myTargetSeq, SVtype="duplication", SVstartPos = 11L, SVlength=mySVLen)
  dupSV <- SV(type = "duplication", chrom = "chr5", posL = GenomicRanges::start(myTargetReg) + 11L - 1L, seqLen = mySVLen, refGen = genomeRef)
  
  expect_identical(toString(dupMut[["SVseq"]]), toString(getSVSeq(dupSV)))
  expect_identical(toString(dupMut[["DNAseq"]]), toString(getContextSeq(dupSV, gr = myTargetReg)))
  
  
  delMut <- injectMySVs(myTargetSeq, SVtype="deletion", SVstartPos = 11L, SVlength=mySVLen)
  delSV <- SV(type = "deletion", chrom = "chr5", posL = GenomicRanges::start(myTargetReg) + 11L - 1L, seqLen = mySVLen, refGen = genomeRef)
  
  expect_identical(toString(delMut[["SVseq"]]), toString(getSVSeq(delSV)))
  expect_identical(toString(delMut[["DNAseq"]]), toString(getContextSeq(delSV, gr = myTargetReg)))
  
  
  invMut <- injectMySVs(myTargetSeq, SVtype="inv", SVstartPos = 11L, SVlength=mySVLen)
  invSV <- SV(type = "inv", chrom = "chr5", posL = GenomicRanges::start(myTargetReg) + 11L - 1L, seqLen = mySVLen, refGen = genomeRef)
  
  expect_identical(toString(invMut[["SVseq"]]), toString(getSVSeq(invSV)))
  expect_identical(toString(invMut[["DNAseq"]]), toString(getContextSeq(invSV, gr = myTargetReg)))
  
})
