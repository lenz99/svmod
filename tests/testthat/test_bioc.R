
context("Bioconductor functionality")


test_that("Test that reading BAMs works as expected", {
  
  BAM_FILENAME <- "resource/test_SV0.bam"
  expect_equivalent( file.exists(BAM_FILENAME), TRUE)
  
  bamFileS <- Rsamtools::BamFile(BAM_FILENAME)  #open() seems to disturb readGAlignment.. calls
  bamFileP <- Rsamtools::BamFile(BAM_FILENAME, asMates = TRUE)
  
  
  # middle position (potential SV) for hand-crafted read-pairs
  posGR <- GenomicRanges::GRanges(seqnames="chr5", ranges=IRanges::IRanges(start=170837725L,
                                                                           end=170837775L))
  
  myParam <- Rsamtools::ScanBamParam(flag=Rsamtools::scanBamFlag(isPaired=NA, isUnmappedQuery = FALSE,
                                                                 isProperPair=NA, isSecondaryAlignment = FALSE),
                                     what=setdiff(Rsamtools::scanBamWhat(), c("qwidth", "qual")),  # #"seq",  # SEQ  is needed later on
                                     which=posGR)
  
  
  expect_equal(getReadCountsFromMapping(bamFileS, minMapQ = 0L), 619L)
  expect_equal(getReadCountsFromMapping(bamFileS, minMapQ = 1L), 594L)
  expect_equal(getReadCountsFromMapping(bamFileS, minMapQ = 0L, regionStr = toString(posGR)), 4L) # four read ends map over posGR
  expect_equal(getReadCountsFromMapping(bamFileS, minMapQ = 1L, regionStr = toString(posGR)), 3L) # SE4_SVread has mapping qual:0
  
  
  # read with asMates=FALSE
  mappingScanBamS <- Rsamtools::scanBam(file=bamFileS, param=myParam)[[1L]]
  mappingGAlignmentS <- GenomicAlignments::readGAlignments(file=bamFileS, param = myParam)
  expect_warning(GenomicAlignments::readGAlignmentPairs(file=bamFileS, param = myParam), "not concordant")
  mappingGAlignmentPairsS <- suppressWarnings(GenomicAlignments::readGAlignmentPairs(file=bamFileS, param = myParam))
  
  # read with asMates=TRUE
  mappingScanBamP <- Rsamtools::scanBam(file=bamFileP, param=myParam)[[1L]]
  mappingGAlignmentP <- GenomicAlignments::readGAlignments(file=bamFileP, param = myParam)
  expect_warning(GenomicAlignments::readGAlignmentPairs(file=bamFileP, param = myParam), "not concordant")
  mappingGAlignmentPairsP <- suppressWarnings(GenomicAlignments::readGAlignmentPairs(file=bamFileP, param = myParam))
  
  mappingScanBamS.len <- lengths(mappingScanBamS)["qname"]
  mappingScanBamP.len <- lengths(mappingScanBamP)["qname"]
  
  
  expect_equivalent(mappingScanBamS.len, 4L)
  expect_equivalent(length(mappingGAlignmentS), mappingScanBamS.len)
  expect_equivalent(length(mappingGAlignmentPairsS), 1L) # 1 concordant read pair
  
  # asMates=TRUE will return both read pairs (even if only one has been matched)
  expect_equivalent(mappingScanBamP.len, 6L) # 2 read-pairs, 2 single reads
  expect_equivalent(length(mappingGAlignmentP), mappingScanBamP.len)
  
  expect_identical(mappingGAlignmentPairsS, mappingGAlignmentPairsP)
  
  
  
  # # bam file with a read-pair
  # bamFileS <- Rsamtools::BamFile("resource/test.s.bam")  #open() seems to disturb readGAlignment.. calls
  # bamFileP <- Rsamtools::BamFile("resource/test.s.bam", asMates = TRUE)
  # 
  # # end=78602583: pos within the insert of a paired-end read
  # # end=78602683: pos overlaps one read end of a paired-end read
  # posGR <- GenomicRanges::GRanges(seqnames="chr5", ranges=IRanges::IRanges(start=78602580,
  #                                                                          end=78602683))
  # 
  # myParam <- Rsamtools::ScanBamParam(flag=Rsamtools::scanBamFlag(isPaired=TRUE, isUnmappedQuery = FALSE,
  #                                                                isProperPair=NA, isSecondaryAlignment = FALSE),
  #                                  what=setdiff(Rsamtools::scanBamWhat(), c("qwidth", "qual")),  # #"seq",  # SEQ  is needed later on
  #                                  which=posGR)
  # 
  # 
  # # read with asMates=FALSE
  # mappingScanBamS <- Rsamtools::scanBam(file=bamFileS, param=myParam)[[1L]]
  # mappingGAlignmentS <- GenomicAlignments::readGAlignments(file=bamFileS, param = myParam)
  # mappingGAlignmentPairsS <- GenomicAlignments::readGAlignmentPairs(file=bamFileS, param = myParam)
  # 
  # # read with asMates=TRUE
  # mappingScanBamP <- Rsamtools::scanBam(file=bamFileP, param=myParam)[[1L]]
  # mappingGAlignmentP <- GenomicAlignments::readGAlignments(file=bamFileP, param = myParam)
  # mappingGAlignmentPairsP <- GenomicAlignments::readGAlignmentPairs(file=bamFileP, param = myParam)
  # 
  # expect_equivalent(lengths(mappingScanBamS)["qname"], 1L)
  # expect_equivalent(length(mappingGAlignmentS), 1L)
  # expect_equivalent(length(mappingGAlignmentPairsS), 1L)
  # 
  # # asMates=TRUE will return both read pairs (even if only one has been matched)
  # expect_equivalent(lengths(mappingScanBamP)["qname"], 2L)
  # expect_equivalent(length(mappingGAlignmentP), 2L)
  # expect_equivalent(length(mappingGAlignmentPairsP), 1L)
  # 
  # expect_identical(mappingGAlignmentPairsS, mappingGAlignmentPairsP)
  # 
  
  
  
  # cleanup
  #try(expr = {close(bamFileS); close(bamFileP)}, silent = TRUE)
  
  
})
