


context("SVs and SV-injection")

test_that("Cigar to pos stuff", {
  myTestCigarStr <- c("100M", "8S9M10S", "8I92M", "8I92M4S", "4S9M", "10M4I3M2S")
  myTestResult <- cigarToPos(cigarStr = myTestCigarStr, chrom = "chromTest", pos = rep(100L, length(myTestCigarStr)))
  
  expect_equal( GenomicRanges::start(myTestResult), 100L - c(0, 8, 0, 0, 4, 0) )
  expect_equal( GenomicRanges::end(myTestResult), 100L + c(100L, 9L, 92L, 92L, 9L, 13L) -1L + c(0, 10, 0, 4, 0, 2) )
})

test_that("Test that SV and extracting an SV-context works as expected", {
  # a random SV location
  myPosL <- 3600000L
  mySVLocus <- GenomicRanges::GRanges(seqnames = "chr5",
                                      ranges = IRanges::IRanges(start = myPosL, width = 3L))
  # myContext1 covers the complete SV-locus
  myContext1 <- GenomicRanges::resize(mySVLocus, width = GenomicRanges::width(mySVLocus) + 4L, fix = "center")
  # myContext2 covers the SV-breakpoint (posL) but not the complete SV-locus
  myContext2 <- GenomicRanges::shift(mySVLocus, shift = -1L)
  # myContext3: just the left breakpoint at posL
  myContext3 <- GenomicRanges::resize(mySVLocus, width = 1L, fix = "start")
  
  expect_identical(GenomicRanges::width(myContext1), 3L + 4L)
  expect_identical(GenomicRanges::start(myContext1), myPosL - 2L)
  expect_identical(GenomicRanges::width(myContext2), 3L)
  expect_identical(GenomicRanges::start(myContext2), myPosL - 1L)
  expect_identical(GenomicRanges::width(myContext3), 1L)
  expect_identical(GenomicRanges::start(myContext3), myPosL)
  
  
  SVins <- SV(type="ins", chrom="chr5", posL = myPosL, seqLen = 3L)
  SVdup <- SV(type="dup", chrom="chr5", posL = myPosL, seqLen = 3L)
  SVdel <- SV(type="del", chrom="chr5", posL = myPosL, seqLen = 3L)
  SVinv <- SV(type="inv", chrom="chr5", posL = myPosL, seqLen = 3L)
  
  expect_identical(length(SVins$seq), 3L) # anything length 3
  expect_identical(toString(SVdup$seq), "AGC") # the duplicated sequence
  expect_identical(toString(SVdel$seq), "AGC")
  expect_identical(toString(SVinv$seq), "GCT") #revComp
  
  expect_identical(getLocus(SVins), mySVLocus )
  expect_identical(getLocus(SVdup), mySVLocus )
  expect_identical(getLocus(SVdel), mySVLocus )
  expect_identical(getLocus(SVinv), mySVLocus )
  
  
  # get SV context for a region that covers the relevant SV region
  #
  #    G C A ** A G C ** C C
  
  # insertion is by convention added at left BP, shifting the original reference seq at SV-locus (i.e. behind the left BP) to the right
  expect_identical( toString(getContextSeq(sv = SVins, gr = myContext1)), paste0("CA", SVins$seq, "AGC", "CC"))
  # for duplication, posL points to the left-most base of the duplicated DNA
  expect_identical( toString(getContextSeq(sv = SVdup, gr = myContext1)), paste0("CA", SVdup$seq, "AGC", "CC"))
  expect_identical( toString(getContextSeq(sv = SVdel, gr = myContext1)), paste0("CA", "CC"))
  expect_identical( toString(getContextSeq(sv = SVinv, gr = myContext1)), paste0("CA", SVinv$seq, "CC"))
  
  # get SV context for a region that does not cover the relevant SV region but lies completely off-target
  expect_warning( getContextSeq(sv = SVins, gr = GenomicRanges::shift(myContext1, shift = 50L)) )
  expect_identical( suppressWarnings(toString(getContextSeq(sv = SVins, gr = GenomicRanges::shift(myContext1, shift = 50L)))), "CGCCCAA")
  
  # Context 2: partially covering the SV-locus
  expect_identical(toString(getContextSeq(sv = SVins, gr = myContext2)), paste0("A", SVins$seq, "AG"))
  expect_identical(toString(getContextSeq(sv = SVdup, gr = myContext2)), paste0("A", SVdup$seq, "AG"))
  expect_error(getContextSeq(sv = SVdel, gr = myContext2))
  expect_error(getContextSeq(sv = SVinv, gr = myContext2))
  
  # Context 3: single base at posL
  expect_identical(toString(getContextSeq(sv = SVins, gr = myContext3)), paste0(SVins$seq, "A"))
  expect_identical(toString(getContextSeq(sv = SVdup, gr = myContext3)), paste0(SVdup$seq, "A"))
})

