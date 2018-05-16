methylpy
========

Welcome to the home page of methylpy, a pyhton-based analysis pipeline
for

-  (single-cell) (whole-genome) bisulfite sequencing data
-  (single-cell) NOMe-seq data
-  differential methylation analysis

methylpy is available at
`github <https://github.com/yupenghe/methylpy>`__ and
`PyPI <https://pypi.python.org/pypi/methylpy/>`__.

Note
====

-  Version 1.3 has major changes on options related to mapping. A new
   aligner, minimap2, is supported starting in this version. To
   accommodate this new features, ``--bowtie2`` option is replaced with
   ``--aligner``, which specifies the aligner to use. The parameters of
   ``--build-reference`` function are modified as well.
-  methylpy only considers cytosines that are in uppercase in the genome
   fasta file (i.e. not masked)
-  methylpy was initiated by and built on the work of `Mattew D.
   Schultz <https://github.com/schultzmattd>`__
-  beta version of
   `tutorial <https://github.com/yupenghe/methylpy/blob/methylpy/tutorial/tutorial.md>`__
   is released!

What can methylpy do?
=====================

Processing bisulfite sequencing data and NOMe-seq data
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

-  fast and flexible pipeline for both single-end and paired-end data
-  all the way from raw reads (fastq) to methylation state and/or open
   chromatin readouts
-  also support getting readouts from alignment (BAM file)
-  including options for read trimming, quality filter and PCR duplicate
   removal
-  accept compressed input and generate compressed output
-  support post-bisulfite adaptor tagging (PBAT) data

Calling differentially methylated regions (DMRs)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

-  DMR calling at single cytosine level
-  support comparison across 2 or more samples/groups
-  conservative and accurate
-  useful feature for dealing with low-coverage data by combining data
   of adjacent cytosines

What you want to do
===================

-  `Install methylpy <#install-methylpy>`__
-  `Test methylpy <#test-methylpy>`__
-  `Process data <#process-data>`__
-  `Call DMRs <#call-dmrs>`__
-  `Additional functions for data
   processing <#additional-functions-for-data-processing>`__
-  `Cite methylpy <#cite-methylpy>`__

run ``methylpy -h`` to get a list of functions.

Install methylpy
================

Step 1 - Download methylpy
^^^^^^^^^^^^^^^^^^^^^^^^^^

Easiest way of installing methylpy will be through PyPI by running
``pip install methylpy``. The command ``pip install --upgrade methylpy``
updates methylpy to latest version. Alternatively, methylpy can be
installed through github: enter the directory where you would like to
install methylpy and run

::

    git clone https://github.com/yupenghe/methylpy.git
    cd methylpy/
    python setup.py install

If you would like to install methylpy in path of your choice, run
``python setup.py install --prefix=/USER/PATH/``. Then, try ``methylpy``
and if no error pops out, the setup is likely successful. See `Test
methylpy <#test-methylpy>`__ for more rigorious test. Last, processing
large dataset will require large spare space for temporary files.
Usually, the default directory for temporary files will not meet the
need. You may want to set the ``TMPDIR`` environmental variable to the
(absolute) path of a directory on hard drive with sufficient space (e.g.
``/YOUR/TMP/DIR/``). This can be done by adding the below command to
``~/.bashrc file``: ``export TMPDIR=/YOUR/TMP/DIR/`` and run
``source ~/.bashrc``.

Step 2 - Install dependencies
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

python is required for running methylpy. Both python2 (>=2.7.9) and
python3 (>=3.6.2) will work. methylpy also depends on two python
modules, `numpy <http://www.numpy.org/>`__ and
`scipy <https://www.scipy.org/>`__. The easiest way to get these
dependencies is to install
`anaconda <https://www.anaconda.com/download/>`__.

In addition, some features of methylpy depend on several publicly
available tools (not all of them are required if you only use a subset
of methylpy functions). \*
`cutadapt <http://cutadapt.readthedocs.io/en/stable/installation.html>`__
(>=1.9) for raw read trimming \*
`bowtie <http://bowtie-bio.sourceforge.net/index.shtml>`__ and/or
`bowtie2 <http://bowtie-bio.sourceforge.net/bowtie2/index.shtml>`__ for
alignment \* `samtools <https://github.com/samtools/samtools>`__ (>=1.3)
for alignment result manipulation. Samtools can also be installed using
conda ``conda install -c bioconda samtools`` \*
`Picard <https://broadinstitute.github.io/picard/index.html>`__
(>=2.10.8) for PCR duplicate removal \* java for running Picard (its
path needs to be included in ``PATH`` environment variable) . \*
`wigToBigWig <http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/wigToBigWig>`__
for converting methylpy output to bigwig format

Lastly, if paths to cutadapt, bowtie/bowtie2, samtools and wigToBigWig
are included in ``PATH`` variable, methylpy can run these tools
directly. Otherwise, the paths have to be passed to methylpy as
augments. Path to Picard needs to be passed to methylpy as a parameter
to run PCR duplicate removal.

Optional step - Compile rms.cpp
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

DMR finding requires an executable
``methylpy/methylpy/run_rms_tests.out``, which was compiled from C++
code ``methylpy/methylpy/rms.cpp``. In most cases, the precompiled file
can be used directly. To test this, simply run execute
``methylpy/methylpy/run_rms_tests.out``. If help page shows, recompiling
is not required. If error turns up, the executable needs to be
regenerated by compiling ``rms.cpp`` and this step requires
`GSL <https://www.gnu.org/software/gsl/>`__ installed correctly. In most
linux operating system, the below commands will do the job

::

    cd methylpy/methylpy/
    g++ -O3 -l gsl -l gslcblas -o run_rms_tests.out rms.cpp

In Ubuntu (>=16.04), please try the below commands first.

::

    cd methylpy/methylpy/
    g++ -o run_rms_tests.out rms.cpp `gsl-config --cflags —libs`

Lastly, the compiled file ``run_rms_tests.out`` needs to be copied to
the directory where methylpy is installed. You can get the directory by
running the blow commands in python console (``python`` to open a python
console):

::

    import methylpy
    print(methylpy.__file__[:methylpy.__file__.rfind("/")]+"/")

Test methylpy
=============

To test whether methylpy and the dependencies are installed and set up
correctly, run

::

    wget http://neomorph.salk.edu/yupeng/share/methylpy_test.tar.gz
    tar -xf methylpy_test.tar.gz
    cd methylpy_test/
    python run_test.py

The test should take around 3 minutes, and progress will be printed on
screen. After the test is started, two files ``test_output_msg.txt`` and
``test_error_msg.txt`` will be generated. The former contains more
details about each test and the later stores error message (if any) as
well as additional information.

If test fails, please check ``test_error_msg.txt`` for the error
message. If you decide to submit an issue regarding test failure to
methylpy github page, please include the error message in this file.

Process data
============

Please see
`tutorial <https://github.com/yupenghe/methylpy/blob/methylpy/tutorial/tutorial.md>`__.
for more details.

Step 1 - Build converted genome reference
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Build bowtie/bowtie2 index for converted genome. Run
``methylpy build-reference -h`` to get more information. An example of
building mm10 mouse reference index:

::

    methylpy build-reference \
        --input-files mm10_bt2/mm10.fa \
        --output-prefix mm10_bt2/mm10 \
        --bowtie2 True

Step 2 - Process bisulfite sequencing and NOMe-seq data
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Function ``single-end-pipeline`` is For processing single-end data. Run
``methylpy single-end-pipeline -h`` to get help information. Below code
is an example of using methylpy to process single-end bisulfite
sequencing data. For processing NOMe-seq data, please use
``num_upstr_bases=1`` to include one base upstream cytosine as part of
cytosine sequence context, which can be used to tease out GC sites.

::

    methylpy single-end-pipeline \
        --read-files raw/mESC_R1.fastq.gz \
        --sample mESC \
        --forward-ref mm10_bt2/mm10_f \
        --reverse-ref mm10_bt2/mm10_r \
        --ref-fasta mm10_bt2/mm10.fa \
        --num-procs 8 \
        --remove-clonal True \
        --path-to-picard="picard/"

An command example for processing paired-end data. Run
``methylpy paired-end-pipeline -h`` to get more information.

::

    methylpy paired-end-pipeline \
        --read1-files raw/mESC_R1.fastq.gz \
        --read2-files raw/mESC_R2.fastq.gz \
        --sample mESC \
        --forward-ref mm10_bt2/mm10_f \
        --reverse-ref mm10_bt2/mm10_r \
        --ref-fasta mm10_bt2/mm10.fa \
        --num-procs 8 \
        --remove-clonal True \
        --path-to-picard="picard/"

If you would like methylpy to perform binomial test for teasing out
sites that show methylation above noise level (which is mainly due to
sodium bisulfite non-conversion), please check options ``--binom-test``
and ``--unmethylated-control``.

Output format
^^^^^^^^^^^^^

Output file(s) are (compressed) tab-separated text file(s) in allc
format. "allc" stands for all cytosine (C). Each row in an allc file
corresponds to one cytosine in the genome. An allc file contain 7
mandatory columns and no header. Two additional columns may be added
with ``--add-snp-info`` option when using ``single-end-pipeline``,
``paired-end-pipeline`` or ``call-methylation-state`` methods.

+---------+----------+----------+--------+
| index   | column   | example  | note   |
|         | name     |          |        |
+=========+==========+==========+========+
| 1       | chromoso | 12       | with   |
|         | me       |          | no     |
|         |          |          | "chr"  |
+---------+----------+----------+--------+
| 2       | position | 18283342 | 1-base |
|         |          |          | d      |
+---------+----------+----------+--------+
| 3       | strand   | +        | either |
|         |          |          | + or - |
+---------+----------+----------+--------+
| 4       | sequence | CGT      | can be |
|         | context  |          | more   |
|         |          |          | than 3 |
|         |          |          | bases  |
+---------+----------+----------+--------+
| 5       | mc       | 18       | count  |
|         |          |          | of     |
|         |          |          | reads  |
|         |          |          | suppor |
|         |          |          | ting   |
|         |          |          | methyl |
|         |          |          | ation  |
+---------+----------+----------+--------+
| 6       | cov      | 21       | read   |
|         |          |          | covera |
|         |          |          | ge     |
+---------+----------+----------+--------+
| 7       | methylat | 1        | indica |
|         | ed       |          | tor    |
|         |          |          | of     |
|         |          |          | signif |
|         |          |          | icant  |
|         |          |          | methyl |
|         |          |          | ation  |
|         |          |          | (1 if  |
|         |          |          | no     |
|         |          |          | test   |
|         |          |          | is     |
|         |          |          | perfor |
|         |          |          | med)   |
+---------+----------+----------+--------+
| 8       | (optiona | 3,2,3    | number |
|         | l)       |          | of     |
|         | num\_mat |          | match  |
|         | ches     |          | baseca |
|         |          |          | lls    |
|         |          |          | at     |
|         |          |          | contex |
|         |          |          | t      |
|         |          |          | nucleo |
|         |          |          | tides  |
+---------+----------+----------+--------+
| 9       | (optiona | 0,1,0    | number |
|         | l)       |          | of     |
|         | num\_mis |          | mismat |
|         | matches  |          | ches   |
|         |          |          | at     |
|         |          |          | contex |
|         |          |          | t      |
|         |          |          | nucleo |
|         |          |          | tides  |
+---------+----------+----------+--------+

Call DMRs
=========

This function will take a list of compressed/uncompressed allc files
(output files from methylpy pipeline) as input and look for DMRs. Help
information of this function is available via running
``methylpy DMRfind -h``.

Below is the code of an example of calling DMRs for CG methylation
between two samples, ``AD_HT`` and ``AD_IT`` on chromosome 1 through 5
using 8 processors.

::

    methylpy DMRfind \
        --allc-files allc/allc_AD_HT.tsv.gz allc/allc_AD_IT.tsv.gz \
        --samples AD_HT AD_IT \
        --mc-type "CGN" \
        --chroms 1 2 3 4 5 \
        --num-procs 8 \
        --output-prefix DMR_HT_IT

Please see
`tutorial <https://github.com/yupenghe/methylpy/blob/methylpy/tutorial/tutorial.md>`__
for details.

Additional functions for data processing
========================================

Extract cytosine methylation state from BAM file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The ``call-methylation-state`` function allows users to get cytosine
methylation state (allc file) from alignment file (BAM file). It is part
of the data processing pipeline which is especially useful for getting
the allc file from alignment file from other methylation data pipelines
like bismark. Run ``methylpy call-methylation-state -h`` to get help
information. Below is an example of running this function. Please make
sure to remove ``--paired-end True`` or use ``--paired-end False`` for
BAM file from single-end data.

::

    methylpy call-methylation-state \
        --input-file mESC_processed_reads_no_clonal.bam \
        --paired-end True \
        --sample mESC \
        --ref-fasta mm10_bt2/mm10.fa \
        --num-procs 8

Get methylation level for genomic regions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Calculating methylation level of certain genomic regions can give an
estimate of the methylation abundance of these loci. This can be
achieved using the ``add-methylation-level`` function. See
``methylpy add-methylation-level -h`` for more details about the input
format and available options.

::

    methylpy add-methylation-level \
        --input-tsv-file DMR_AD_IT.tsv \
        --output-file DMR_AD_IT_with_level.tsv \
        --allc-files allc/allc_AD_HT_1.tsv.gz allc/allc_AD_HT_2.tsv.gz \
            allc/allc_AD_IT_1.tsv.gz allc/allc_AD_IT_2.tsv.gz \
        --samples AD_HT_1 AD_HT_2 AD_IT_1 AD_IT_2 \
        --mc-type CGN \
        --num-procs 4

Merge allc files
^^^^^^^^^^^^^^^^

The ``merge-allc`` function can merge multiple allc files into a single
allc file. It is useful when separate allc files are generated for
replicates of a tissue or cell type, and one wants to get a single allc
file for that tissue/cell type. See ``methylpy merge-allc -h`` for more
information.

::

    methylpy merge-allc \
        --allc-files allc/allc_AD_HT_1.tsv.gz allc/allc_AD_HT_2.tsv.gz \
        --output-file allc/allc_AD_HT.tsv.gz \
        --num-procs 1 \
        --compress-output True

Filter allc files
^^^^^^^^^^^^^^^^^

The ``filter-allc`` function is for filtering sites by cytosine context,
coverage etc. See ``methylpy filter-allc -h`` for more information.

::

    methylpy filter-allc \
        --allc-file allc/allc_AD_HT_1.tsv.gz \
        --output-file allc/allCG_AD_HT_1.tsv.gz \
        --mc-type CGN \
        --min-cov 2 \
        --compress-output True

Index allc files
^^^^^^^^^^^^^^^^

The ``index-allc`` function allows creating index file for each allc
file. The index file can be used for speeding up allc file reading
similar to the .fai file for .fasta file. See ``methylpy index-allc -h``
for more information.

::

    methylpy index-allc \
        --allc-files allc/allc_AD_HT_1.tsv.gz allc/allc_AD_HT_2.tsv.gz \
        --num-procs 2 \
        --no-reindex False

Convert allc file to bigwig format
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The ``allc-to-bigwig`` function generates bigwig file from allc file.
Methylation level will be calculated in equally divided non-overlapping
genomic bins and the output will be stored in a bigwig file. See
``methylpy allc-to-bigwig -h`` for more information.

::

    methylpy allc-to-bigwig \
        --allc-file results/allc_mESC.tsv.gz \
        --output-file results/allc_mESC.bw \
        --ref-fasta mm10_bt2/mm10.fa \
        --mc-type CGN \
        --bin-size 100  

Quality filter for bisulfite sequencing reads
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Sometimes, we want to filter out reads that cannot be mapped confidently
or are likely from under-converted DNA fragments. This can be done using
the ``bam-quality-filter`` function. See
``methylpy bam-quality-filter -h`` for parameter inforamtion.

For example, below command can be used to filter out reads with less
than 30 MAPQ score (poor alignment) and with mCH level greater than 0.7
(under-conversion) if the reads contain enough (at least 3) CH sites.

::

    methylpy bam-quality-filter \
        --input-file mESC_processed_reads_no_clonal.bam \
        --output-file mESC_processed_reads_no_clonal.filtered.bam \
        --ref-fasta mm10_bt2/mm10.fa \
        --min-mapq 30 \
        --min-num-ch 3 \
        --max-mch-level 0.7 \
        --buffer-line-number 100

Reidentify DMRs from existing result
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

methylpy is able to reidentify-DMR based on the result of previous
DMRfind run. This function is especially useful in picking out DMRs
across a subset of categories and/or with different filters. See
``methylpy reidentify-DMR -h`` for details about the options.

::

    methylpy reidentify-DMR \
        --input-rms-file results/DMR_P0_FBvsHT_rms_results.tsv.gz \
        --output-file results/DMR_P0_FBvsHT_rms_results_recollapsed.tsv \
        --collapse-samples P0_FB_1 P0_FB_2 P0_HT_1 P0_HT_2 \
        --sample-category P0_FB P0_FB P0_HT P0_HT \
        --min-cluster 2

Cite methylpy
=============

If you use methylpy, please cite >Matthew D. Schultz, Yupeng He, John
W.Whitaker, Manoj Hariharan, Eran A. Mukamel, Danny Leung, Nisha
Rajagopal, Joseph R. Nery, Mark A. Urich, Huaming Chen, Shin Lin, Yiing
Lin, Bing Ren, Terrence J. Sejnowski, Wei Wang, Joseph R. Ecker. Human
Body Epigenome Maps Reveal Noncanonical DNA Methylation Variation.
Nature. 523(7559):212-216, 2015 Jul.
