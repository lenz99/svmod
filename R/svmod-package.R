#' SVmod. A package to find differential SV in tumor/normal samples through statistical machine-learning 'model'.
#' 
#' Uses classifiers methods (like one-class SVM)  to predict somatic (=differential) SVs. It is tuned to work with exome NGS-data.
#' 
#' @name svmod
#' @docType package
#' @useDynLib svmod
#' @import Biostrings bitops doParallel foreach ggplot2 GenomicRanges GenomicAlignments IRanges logging mlr RcppRoll Rsamtools rtracklayer XVector
#' @importFrom dplyr group_by select filter arrange n summarise full_join
#' @importFrom magrittr %<>% %>% 
NULL

