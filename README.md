# Find small somatic structural variants in tumor/normal exome data

Structural variants (SV) are genomic variation at the DNA-level that affects more than 50 nt. Typically, SVs range in size from one kilobase up to few megabases. `svmod` implements a machine-learning approach to SV-finding, more specifically, it is about small (50-300bp) somatic (tumor-only) SVs that are detectable in exome sequencing data. You can find more on this approach at [our publication](http://rdcu.be/niI9).

Currently, `svmod` is still in a prototype stage with many rough edges (see below).. So, you have been warned!


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


## First Steps

The `svmod` package ships with some R-scripts in its `exec` directory that make the functionality of the package more accessible.


### Simulation

1. `sim_patGenerator.R` initiates a simulation by creating tumor and normal genomic background of virtual patients. Somatic SVs will be introduced for each patient according to user settings. NGS reads are simulated in a given coverage and subsequently mapped to the reference.
2. `sim_extractFeatures.R` extracts features of the local mapping patterns from a simulated patient's mapping. Feature extraction can directly use the information of the simulation where somatic SV reside (`--useSimInfo`) or rely on clipped base screening. Use the former to evaluate only the confirmation step while the latter mimicks the situation of real data when no somatic SV information is available. Feature data is saved to disk as R-dataframe.
3. `benchmarkLearners.R` trains and tunes different classification methods, namely k-NN, random forest, elastic net, SVM and logistic regression. These classifiers are then evaluated, either by repeatedly subsampling from training data or by using the screen-positive regions of the patient as natural hold-out data.


### Real data
For real patient tumor/normal samples the entry point are mapping files.

1. `real_extractFeatures.R` script runs the screening step and extracts features of local mapping patterns. Consult this script to understand the naming conventions that are in place to identify the tumor and normal mapping BAM-files. Currently, test data (i.e. screen-positive regions) is annotated according to the FLT3-ITD that is discussed in the publication.
2. `benchmarkLearners.R` to train, tune and evaluate different learners. See above at the simulation section.


## Known bugs/limitations
1. No clean interface yet to run `svmod`. In particular for real patient data you have to fiddle with script internas.
2. NGS-read simulation is hard-coded to use the Illumina HiSeq 2000. Hence, the simulated reads are limited to 100bp.
