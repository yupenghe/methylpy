# Tutorial for methylpy
A pyhton-based analysis pipeline for
* (single-cell) (whole-genome) bisulfite sequencing data
* (single-cell)  NOMe-seq data

## Installation
Enter the directory where you would like to install methylpy and run
`git clone https://github.com/yupenghe/methylpy.git`

## methylpy dependencies
methylpy depends on two python modules, numpy and scipy
methylpy also depends on 
cutadapt (>=1.12) for trimming
samtools (>=1.3)
bowtie and/or bowtie2
picard (>=2.10.8)

To install numpy and scipy, the easiest way is to install anaconda. If you have the methylpy cloned and have these dependencies installed, the last step is to include methylpy into python library search path list. Supposing that your methylpy clone is at /YOUR/PATH/methylpy/, then you can include the command below in your ~/.bashrc file:
“export PYTHONPATH=/YOUR/PATH/methylpy/:$PYTHONPATH”
Also, make sure the paths to cutadapt, samtools, bowtie and/or bowtie2 are included in PATH. Otherwise, you will have to pass these path to methylpy functions when running the pipeline.

## Checking
To check whether your installation is successful, try:
```
python2.7
>>>import methylpy.call_mc
>>>import methylpy.DMRfind
```
If you get no error, there is a good chance that methylpy is correctly installed. However, these commands do not check the dependencies. 

## Running methylpy for paired-end MethylC-seq read mapping
Step -1 Set directory for temporary files
Before running methylpy, please set the directory for storing temporary files to be a directory with large space by adding below command to ~/.bashrc file: 
export TMPDIR=/YOUR/TMP/DIR/

### Step 0 Prepare (converted) reference genome for mapping
Include the sequence of lambda phage genome (as control) into the reference genome fasta file. For example: 
cat mm10.fa chrL.fa > mm10_with_lambda.fa
chrL.fa can be downloaded from [here](http://neomorph.salk.edu/yupeng/share/chrL.fa).

### Step 1 Build bowtie/bowtie2 index for converted genome
Build the bowtie/bowtie2 index for bisulfite-converted genome. Note that the genome index used for WGBS data processing is different from the one used for ChIP-seq or genome sequence data processing.  
```
from methylpy.call_mc import build_ref
#fasta file(s) of genome                                                        
#Specifying multiple files like input_files=['chr1.fa','chr2.fa',...,'chrY.fa','chrL.fa'] also #work
input_files=['mm10_with_lambda.fa']
#Prefix of output files                                                      
output='mm10'

#Build bowtie2 index
build_ref(input_files,output,bowtie2=True)
	
#Build bowtie index
build_ref(input_files,output,bowtie2=False)
```	
	
### Step 2.1 Create script to run methylpy (paired-end data)
Create a run_call_mc.py file like for each sample: 
```	
from methylpy.call_mc import run_methylation_pipeline_pe
	
#Sample name
sample="Test_WGBS"
	
#Path to fastq files
fastq_R1 = ["lib1/*_1.fastq"," lib2/*_1.fastq"] #read 1
fastq_R2 = ["lib1/*_2.fastq"," lib2/*_2.fastq"] #read 2
	
#Assoicate fastq files with corresponding library names
libraries=["lib1","lib2"]
	
#Genome (take mm10 as example)
f_ref = "/home/users/mm10/mm10_f"
r_ref = "/home/users/mm10/mm10_r"
ref_fasta = "/home/users/mm10/mm10.fa"
	
#Number of processors
num_procs = 8
	
#Adapter to trim
adapter_seq_R1='AGATCGGAAGAGCACACGTCTGAAC' adapter_seq_R2='AGATCGGAAGAGCGTCGTGTAGGGA'
	
#Control (unmethylated lambda genome) to estimate bisulfite conversion efficiency #which will be measured by bisulfite non-conversion rate
m_control="chrL:"
	

read1_files = fastq_R1, 
read2_files = fastq_R2, 
libraries=libraries, 
sample=sample,
    forward_reference = f_ref,
		reverse_reference = r_ref,
			reference_fasta = ref_fasta,
				unmethylated_control=m_control, 
					path_to_samtools="",
						path_to_bowtie="",                    
							num_procs=num_procs,
	trim_reads=True,
	path_to_cutadapt="", 
	adapter_seq_R1= adapter_seq_R1,
	adapter_seq_R2= adapter_seq_R2,
	sig_cutoff=0.01, 
		binom_test=True, 
			bh=True,
				sort_mem='1G',
					path_to_MarkDuplicates=”/YOURPATH/picard/”
)
```

### Step 2.2 Create script to run methylpy (for single-end data)
```
#Start of the script: run_call_mc.py
from methylpy.call_mc import run_methylation_pipeline

#Sample name
sample="Test_WGBS"

#Path to fastq files
files=["lib1/*.fastq","lib2/*.fastq"]

#Assoicate fastq files with corresponding library names
libraries=["lib1","lib2"]

#Genome (take mm10 as example)
f_ref = "/home/users/mm10/mm10_f"
r_ref = "/home/users/mm10/mm10_r"
ref_fasta = "/home/users/mm10/mm10.fa"

#Number of processors
num_procs = 8

#Adapter to trim
adapter = 'AGATCGGAAGAGCTCGTATGCC'

#Control (unmethylated lambda genome) to estimate bisulfite conversion efficiency which will be measured by bisulfite non-conversion rate
m_control="chrL:"

run_methylation_pipeline(
files=files, 
libraries=libraries, 
sample=sample,
forward_reference = f_ref,
reverse_reference = r_ref,
reference_fasta = ref_fasta,
unmethylated_control=m_control, 
path_to_samtools="", 
path_to_bowtie="",
num_procs=num_procs, 
path_to_cutadapt ="", 
adapter_seq=adapter_seq,
sig_cutoff=0.01, 
binom_test=True, 
bh=True,
sort_mem='1G',
path_to_MarkDuplicates=”/YOURPATH/picard/”
)
```
### Step 3 Running methylpy and read output
Copy the script to output directory and run it with command “python2.7 run_call_mc.py > log 2> err”. Printed output will be store in the log file and error messages will be store in the err file. If everything is correct, the allc files (final output from methylpy) will be created in the directory. allc files are tab-delimited files containing seven columns, which correspond to chromosome, position, strand, methylation class, methylation count, total coverage, and methylation call.
						 
### Step 4 Remove some temporary files to free up space
Remove *mpileup* files in the output directory. These intermediate files in total take usually huge space (hundred gigabytes) and they can be regenerated using bam file
						 
### Step 5 Calling differentially methylated regions
You can use the below scripts to get regions with different methylation patterns in different samples.
```
#Start of the script: run_call_mc.py	
from methylpy.DMRfind import DMRfind

#Samples included in the comparison
samples = ["Test_WGBS_1”,” Test_WGBS_2”]

#Regions included in the DMR analysis
chrom = map(str,range(1,20))
chrom.extend('X')
chrom.extend('Y')
region_dict={}
for chr in chrom:
	region_dict[chr]=[0,5000000000000000]

#methylation type
mc_type=['CGN']

#Number of processors
num_procs=4

path_to_allc = "/path/where/you/store/all/allc files/”
output_prefix = "DMR_CG_Test_WGBS"

DMRfind(mc_type=mc_type,region_dict=region_dict,samples=samples,
	path_to_allc= path_to_allc,
	num_sims=3000, 
	num_sig_tests=100,
	use_mc_status = False,
	num_procs = num_procs,
	save_result = output_prefix,
	dmr_max_dist=250)                                                                              
```
							 
							 
### Hard drive space and memory requirement
Hard drive space usually space that is three times of the size of all uncompressed fastq files is enough for running the pipeline. 
Memory
methylpy is memory efficient but the sort function (from linux system) could potentially occupy memory. This can be control by the “sort_mem” option (see below). For example, if a node has 24G memory and 16 CPUs, the “sort_mem” option should be set to be “1G” if all 16 CPUs are used. In that case, methylpy will take a most 16 * 1G = 16G memory. Please set the value of sort_mem according to the capacity of your computer/cluster.
