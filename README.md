# methylpy
Welcome to the home page of methylpy, a pyhton-based analysis pipeline for
* (single-cell) (whole-genome) bisulfite sequencing data
* (single-cell)  NOMe-seq data

## What's in methylpy
* Processing bisulfite sequencing data and NOMe-seq data - from raw reads (fastq) to methylation/open chromatin read out
* Calling differentially methylated regions (DMRs) across 2 or more samples


## Compiling rms.cpp
* Most cases
g++ -O3 -l gsl -l gslcblas -o run_rms_tests.out rms.cpp

* Ubuntu 16.04
g++ -o run_rms_tests.out rms.cpp `gsl-config --cflags —libs`
