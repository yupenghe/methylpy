# methylpy
Welcome to the home page of methylpy, a pyhton-based analysis pipeline for
* (single-cell) (whole-genome) bisulfite sequencing data
* (single-cell)  NOMe-seq data

## Note
* Currently, methylpy only supports python2.7 but support for python3 is coming.
* methylpy has major changes compared to previous version. Please double check whether your code is compatible.

## What can methylpy do?
Processing bisulfite sequencing data and NOMe-seq data
* fast and flexible pipeline for both single-end and paired-end data
* all the way from raw reads (fastq) to methylation state and/or open chromatin readouts
* also support getting readouts from alignment (BAM file)
* including options for read trimming, quality filter and PCR duplicate removal
* accept compressed input and generate compressed output
* support post-bisulfite adaptor tagging (PBAT) data

Calling differentially methylated regions (DMRs)
* DMR calling at single cytosine level
* support comparison across 2 or more samples/groups
* conservative and accurate
* useful feature for dealing with low-coverage data by combining data of adjacent cytosines

## Installation
#### step 1 - download methylpy and set up environment variable
Enter the directory where you would like to install methylpy and run
```
git clone https://github.com/yupenghe/methylpy.git
```
Next, the methylpy folder needs to be included in the python search path. Also, methylpy executable will easy to use if its path is included in `PATH`. To do these, if the path of methylpy (not the methylpy/methylpy folder) is /YOUR/PATH/methylpy/, please include the below code in the $HOME/.bashrc file:
```
export PYTHONPATH=/YOUR/PATH/methylpy/:$PYTHONPATH
export PATH=/YOUR/PATH/methylpy/bin/:$PATH
```
and then do
```
source $HOME/.bashrc
```
If the below code gives no error, the setup is successful.
```
methylpy
```

#### step 2 - install dependencies
methylpy is written in python2 so obviously python2 needs to be installed.
methylpy also depends on two python modules, [numpy](http://www.numpy.org/) 
and [scipy](https://www.scipy.org/).
The easiest way to resolve these dependencies is to install [anaconda](https://www.anaconda.com/download/)

In addition, some features of methylpy depend on several publicly available tools
* [cutadapt](http://cutadapt.readthedocs.io/en/stable/installation.html) (>=1.12) for raw read trimming
* [bowtie](http://bowtie-bio.sourceforge.net/index.shtml) and/or [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) for alignment
* [samtools](https://github.com/samtools/samtools) (>=1.3) for alignment results manipulation
* [Picard](https://broadinstitute.github.io/picard/index.html) (>=2.10.8) for removal of PCR duplicates

Lastly, if paths to cutadapt, bowtie/bowtie2 and samtools are included in `PATH` variable,
methylpy can run these tools directly. Otherwise, the paths have to be passed to methylpy as augments. 
Path to Picard needs to be passed to methylpy as a parameter to run PCR duplicate removal.

#### Compiling rms.cpp
* Most cases
g++ -O3 -l gsl -l gslcblas -o run_rms_tests.out rms.cpp

* Ubuntu 16.04
g++ -o run_rms_tests.out rms.cpp `gsl-config --cflags —libs`

## Using methylpy

