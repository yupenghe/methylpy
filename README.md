# methylpy
Welcome to the home page of methylpy, a pyhton-based analysis pipeline for
* (single-cell) (whole-genome) bisulfite sequencing data
* (single-cell)  NOMe-seq data

## What can methylpy do
Processing bisulfite sequencing data and NOMe-seq data
* fast and flexible pipeline for both single-end and paired-end data
* all the way from raw reads (fastq) to methylation state and/or open chromatin readouts
* including options for read trimming, quality filter and PCR duplicate removal
* accept compressed input and generate compressed output
* support post-bisulfite adaptor tagging (PBAT) data

Calling differentially methylated regions (DMRs)
* DMR calling at single cytosine level
* support comparison across 2 or more samples/groups
* conservative and accurate
* useful feature for dealing with low-coverage data by combining data of adjacent cytosines

## Installation

## Compiling rms.cpp
* Most cases
g++ -O3 -l gsl -l gslcblas -o run_rms_tests.out rms.cpp

* Ubuntu 16.04
g++ -o run_rms_tests.out rms.cpp `gsl-config --cflags —libs`
