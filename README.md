# svmod
Search small somatic structural variations (SV) in tumor/normal exome data

`svmod` implements a machine-learning approach to SV-finding, more specifically, it is about small (50-300bp) somatic (tumor-only) SVs that are detectable in exome sequencing data.


## Requirements

`svmod` has been tested on Mac and on Linux.
Please check the file `init.R` where important directory paths are set.

It is an R-package that requires some Bioconductor packages and third-party tools. 

### Third-party tools
  * BWA
  * ART
  * samtools
  * bedtools
  * a reference genome fasta

### Bioconductor packages
  * GenomicAlignments
  * GenomicFiles
  * GenomicRanges
  * Rsamtools
  * rtracklayer

