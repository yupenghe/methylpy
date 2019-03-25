import sys
import os
import multiprocessing
import subprocess
import scipy.stats as sci
from scipy.stats.mstats import mquantiles
from methylpy.utilities import print_checkpoint, print_error, print_warning
from methylpy.utilities import split_fastq_file
from methylpy.utilities import split_fastq_file_pbat
from methylpy.utilities import open_allc_file,index_allc_file
from methylpy.utilities import read_allc_index,bgzip_allc_file
from methylpy.utilities import check_call_mc_dependencies
import pdb
import shlex
import itertools
import re
import glob
import io as cStr
import bisect
import gzip
import math

def run_methylation_pipeline(read_files, sample,
                             forward_reference, reverse_reference, reference_fasta,
                             libraries = None,
                             unmethylated_control=None,
                             path_to_output="", sig_cutoff=0.01,
                             num_procs=1, sort_mem="500M",
                             num_upstr_bases=0, num_downstr_bases=2,
                             generate_allc_file=True,
                             generate_mpileup_file=True,
                             compress_output=True,
                             bgzip=False,
                             path_to_bgzip="",
                             path_to_tabix="",
                             binom_test=False, min_cov=2,
                             trim_reads=True, path_to_cutadapt="",
                             pbat=False,check_dependency=True,
                             path_to_aligner="",
                             aligner="bowtie2",aligner_options=None,
                             merge_by_max_mapq=False,min_mapq=30,
                             remove_clonal=False,keep_clonal_stats=True,
                             path_to_picard="",java_options="-Xmx20g",
                             path_to_samtools="",
                             remove_chr_prefix=True,
                             add_snp_info=False,
                             adapter_seq="AGATCGGAAGAGCACACGTCTG",
                             max_adapter_removal=None,
                             overlap_length=None, zero_cap=None,
                             error_rate=None, min_qual_score=10,
                             min_read_len=30,
                             keep_temp_files=False,
                             min_base_quality=1):

    """
    read_files is a list of all the fastq files you'd like to run through the pipeline.
        Note that globbing is supported here (i.e., you can use * in your paths)

    libraries is a list of library IDs (in the same order as the files list) indiciating which
        libraries each set of fastq files belong to. If you use a glob, you only need to indicate
        the library ID for those fastqs once (i.e., the length of files and libraries should be
        the same)

    sample is a string indicating the name of the sample you're processing. It will be included
        in the output files.

    forward_reference is a string indicating the path to the forward strand reference created by
        build_ref

    reverse_reference is a string indicating the path to the reverse strand reference created by
        build_ref

    reference_fasta is a string indicating the path to a fasta file containing the sequences
        you used for mapping
        input is the path to a bam file that contains mapped bisulfite sequencing reads

    unmethylated_control is the name of the chromosome/region that you want to use to estimate
        the non-conversion rate of your sample, or the non-conversion rate you'd like to use.
        Consequently, control is either a string, or a decimal.
        If control is a string then it should be in the following format: "chrom:start-end".
        If you'd like to specify an entire chromosome simply use "chrom:"

    remove_clonal is a boolean indicating that you want to remove clonal reads (PCR duplicates).
        If true, picard.jar should be available in folder specified in path_to_picard.

    path_to_picard is a string of the path to "picard.jar". "picard.jar" is assumed to be
        in your path if this option isn't used.

    path_to_samtools is a string indicating the path to the directory containing your
        installation of samtools. Samtools is assumed to be in your path if this is not
        provided.

    path_to_aligner is a string indicating the path to the folder in which bowtie resides. Bowtie
        is assumed to be in your path if this option isn't used.

    aligner_options is a list of strings indicating options you'd like passed to the aligner.

    num_procs is an integer indicating how many num_procs you'd like to run this function over.

    trim_reads is a boolean indicating that you want to have reads trimmed by cutadapt.

    path_to_cutadapt is the path to the cutadapt execuatable. Otherwise this is assumed to be in 
        your path.

    adapter_seq is the sequence of an adapter that was ligated to the 3' end. The adapter itself and
        anything that follows is trimmed.

    max_adapter_removal indicates the maximum number of times to try to remove adapters. Useful when
        an adapter gets appended multiple times.

    overlap_length is the minimum overlap length. If the overlap between the read and the adapter is
        shorter than LENGTH, the read is not modified. This reduces the no. of bases trimmed purely
        due to short random adapter matches.

    zero_cap causes negative quality values to be set to zero (workaround to avoid
        segmentation faults in BWA).

    error_rate is the maximum allowed error rate (no. of errors divided by the length of the
        matching region). Default: 0.1

    min_qual_score allows you to trim low-quality ends from reads before adapter removal. The
        algorithm is the same as the one used by BWA (Subtract CUTOFF from all qualities; compute
        partial sums from all indices to the end of the sequence; cut sequence at the index at
        which the sum is minimal).

    min_read_len indicates the minimum length a read must be to be kept. Reads that are too short even
        before adapter removal are also discarded. In colorspace, an initial primer is not counted.

    sig_cutoff is a float indicating the adjusted p-value cutoff you wish to use for determining 
        whether or not a site is methylated

    min_cov is an integer indicating the minimum number of reads for a site to be tested.

    binom_tests indicates that you'd like to use a binomial test, rather than the alternative method
        outlined here: https://bitbucket.org/schultzmattd/methylpy/wiki/Methylation%20Calling

    keep_temp_files is a boolean indicating that you'd like to keep the intermediate files generated
        by this function. This can be useful for debugging, but in general should be left False.

    bowtie2 specifies whether to use the bowtie2 aligner instead of bowtie

    sort_mem is the parameter to pass to unix sort with -S/--buffer-size command

    path_to_output is the path to a directory where you would like the output to be stored.
        The default is the same directory as the input fastqs.

    min_base_quality is an integer indicating the minimum PHRED quality score for a base to be
        included in the mpileup file (and subsequently to be considered for methylation calling).
    """

    if check_dependency:
        check_call_mc_dependencies(path_to_samtools,
                                   trim_reads,
                                   path_to_cutadapt,
                                   path_to_aligner,
                                   remove_clonal,
                                   path_to_picard)

    if libraries is None:
        libraries = ["libA"]

    if not isinstance(libraries, list):
        if isinstance(libraries, str):
            libraries = [libraries]
        else:
            exit("libraries must be a list of string(s)")

    if len(libraries) == 1 and len(read_files) > 1:
        uniq_library = libraries[0]
        libraries = [uniq_library for ind in range(len(read_files))]

    #Default bowtie option
    if aligner_options is None:
        if aligner.lower() == "minimap2":
            aligner_options = ["-ax","sr","--secondary=no"]
        elif aligner.lower() == "bowtie":
            aligner_options = ["-S", "-k 1", "-m 1", "--chunkmbs 3072",
                               "--best", "--strata", "-o 4", "-e 80",
                               "-l 20", "-n 0"]
            aligner_options.append("--phred33-quals")
        else: # bowtie2
            aligner_options = []
            aligner_options.append("--phred33-quals")

    # CASAVA >= 1.8
    quality_base = 33

    if len(path_to_samtools) != 0:
        path_to_samtools += "/"

    if len(path_to_aligner) != 0:
        path_to_aligner += "/"

    # output path
    if len(path_to_output) != 0:
        path_to_output += "/"
        if not os.path.exists(path_to_output):
            try:
                os.makedirs(path_to_output)
            except:
                print_error("  Failed to create output folder!")

    expanded_file_list = []
    expanded_library_list = []

    total_input = 0
    total_unique = 0
    total_clonal = 0
    for path, library in zip(read_files, libraries):
        glob_list = glob.glob(path)
        for filen in glob_list:
            expanded_file_list.append(filen)
            expanded_library_list.append(library)

    for current_library in set(libraries):
        library_files = [filen for filen, library
                         in zip(expanded_file_list, expanded_library_list)
                         if library == current_library]

        #deal with actual filename rather than path to file
        lib_input, lib_unique = run_mapping(current_library, library_files, sample,
                                            forward_reference, reverse_reference, reference_fasta,
                                            path_to_output=path_to_output,
                                            path_to_samtools=path_to_samtools,
                                            path_to_aligner=path_to_aligner,
                                            aligner=aligner,
                                            aligner_options=aligner_options,
                                            merge_by_max_mapq=merge_by_max_mapq,
                                            pbat=pbat,
                                            min_mapq=min_mapq,
                                            num_procs=num_procs,
                                            trim_reads=trim_reads,
                                            path_to_cutadapt=path_to_cutadapt,
                                            adapter_seq=adapter_seq,
                                            max_adapter_removal=max_adapter_removal,
                                            overlap_length=overlap_length, zero_cap=zero_cap,
                                            quality_base=quality_base,
                                            error_rate=error_rate,
                                            min_qual_score=min_qual_score,
                                            min_read_len=min_read_len,
                                            keep_temp_files=keep_temp_files,
                                            sort_mem=sort_mem)
        total_input += lib_input
        total_unique += lib_unique

        ## Remove clonal reads
        if remove_clonal == True:
            lib_clonal = remove_clonal_bam(input_bam=path_to_output+sample+"_"+
                                           str(current_library)+"_processed_reads.bam",
                                           
                                           output_bam=path_to_output+sample+"_"+
                                           str(current_library)+"_processed_reads_no_clonal.bam",
                                           
                                           metric=path_to_output+sample+"_"+
                                           str(current_library)+".metric",
                                           
                                           is_pe=False,
                                           path_to_picard=path_to_picard,
                                           java_options=java_options)

            subprocess.check_call(shlex.split("rm "+path_to_output+sample+"_"+
                                              str(current_library)+"_processed_reads.bam"))
            if not keep_clonal_stats:
                subprocess.check_call(shlex.split("rm "+" "+path_to_output+sample+"_"+
                                                  str(current_library)+".metric"))

            total_clonal += lib_clonal
    print_checkpoint("There are "+str(total_input)+" total input reads")
    print_checkpoint("There are "+str(total_unique)+" uniquely mapping reads, " +
                     str(float(total_unique) / total_input*100)+" percent remaining")

    if remove_clonal == True:
        total_non_clonal = total_unique - total_clonal
        print_checkpoint("There are " + str(total_non_clonal) + " non-clonal reads, " +
                         str(float(total_non_clonal) / total_input*100) + " percent remaining")
        ## Merge bam files to get final bam file
        library_files = [path_to_output+sample+"_"+str(library)+"_processed_reads_no_clonal.bam"
                         for library in set(libraries)]
        if len(library_files) > 1:
            merge_bam_files(library_files, path_to_output+sample+"_processed_reads_no_clonal.bam", path_to_samtools)
            subprocess.check_call(shlex.split("rm "+" ".join(library_files)))
        else:
            subprocess.check_call(shlex.split("mv "+library_files[0]+" "+
                                              path_to_output+sample+"_processed_reads_no_clonal.bam"))
    ## If not removing clonal reads
    else:
        library_files = [path_to_output+sample+"_"+str(library)+
                         "_processed_reads.bam" for library in set(libraries)]
        if len(library_files) > 1:
            merge_bam_files(library_files, path_to_output+sample+"_processed_reads.bam", path_to_samtools)
            subprocess.check_call(shlex.split("rm "+" ".join(library_files)))
        else:
            subprocess.check_call(shlex.split("mv "+library_files[0]+" "+path_to_output+sample+"_processed_reads.bam"))

    if generate_allc_file:
        print_checkpoint("Begin calling mCs")
        if remove_clonal == True:
            output_bam_file = path_to_output+sample+"_processed_reads_no_clonal.bam"
        else:
            output_bam_file = path_to_output+sample+"_processed_reads.bam"
        call_methylated_sites(output_bam_file,
                              sample,
                              reference_fasta,
                              unmethylated_control=unmethylated_control,
                              sig_cutoff=sig_cutoff,
                              num_procs=num_procs,
                              num_upstr_bases=num_upstr_bases,
                              num_downstr_bases=num_downstr_bases,
                              generate_mpileup_file=generate_mpileup_file,
                              compress_output=compress_output,
                              bgzip=bgzip,
                              path_to_bgzip=path_to_bgzip,
                              path_to_tabix=path_to_tabix,
                              min_mapq=min_mapq,
                              min_cov=min_cov,
                              binom_test=binom_test,
                              remove_chr_prefix=remove_chr_prefix,
                              sort_mem=sort_mem,
                              path_to_files=path_to_output,
                              path_to_samtools=path_to_samtools,
                              add_snp_info=add_snp_info,
                              min_base_quality=min_base_quality,
                              keep_temp_files=keep_temp_files)
    print_checkpoint("Done")

def run_mapping(current_library, library_files, sample,
                forward_reference, reverse_reference, reference_fasta,
                path_to_output="",
                path_to_samtools="", path_to_aligner="",
                aligner="bowtie2",
                aligner_options=None,merge_by_max_mapq=False,
                min_mapq=30,pbat=False,
                num_procs=1, trim_reads=True, path_to_cutadapt="",
                adapter_seq="AGATCGGAAGAGCACACGTCTG",
                max_adapter_removal=None, overlap_length=None, zero_cap=None,
                quality_base=None, error_rate=None,
                min_qual_score=10, min_read_len=30,
                keep_temp_files=False,
                sort_mem="500M"):
    """
    This function runs the mapping portion of the methylation calling pipeline.

    current_library is the ID that you'd like to run mapping on.

    libraries is a list of library IDs (in the same order as the files list) indiciating which
        libraries each set of fastq files belong to. If you use a glob, you only need to indicate
        the library ID for those fastqs once (i.e., the length of files and libraries should be
        the same)

    sample is a string indicating the name of the sample you're processing. It will be included
        in the output files.

    forward_reference is a string indicating the path to the forward strand reference created by
        build_ref

    reverse_reference is a string indicating the path to the reverse strand reference created by
        build_ref

    reference_fasta is a string indicating the path to a fasta file containing the sequences
        you used for mapping

    path_to_samtools is a string indicating the path to the directory containing your
        installation of samtools. Samtools is assumed to be in your path if this is not
        provided.

    path_to_aligner is a string indicating the path to the folder in which bowtie resides. Bowtie
        is assumed to be in your path if this option isn't used.

    aligner_options is a list of strings indicating options you'd like passed to bowtie
    
    num_procs is an integer indicating how many num_procs you'd like to run this function over

    trim_reads is a boolean indicating that you want to have reads trimmed by cutadapt.

    path_to_cutadapt is the path to the cutadapt execuatable. Otherwise this is assumed to be in your
        path.

    adapter_seq is the sequence of an adapter that was ligated to the 3' end. The adapter itself and
        anything that follows is trimmed.

    max_adapter_removal indicates the maximum number of times to try to remove adapters. Useful when
        an adapter gets appended multiple times.

    overlap_length is the minimum overlap length. If the overlap between the read and the adapter is
        shorter than LENGTH, the read is not modified. This reduces the no. of bases trimmed purely
        due to short random adapter matches.

    zero_cap causes negative quality values to be set to zero (workaround to avoid segmentation faults
        in BWA).

    quality_base is the offset for quality scores. In other words, assume that quality values are
        encoded as ascii(quality + QUALITY_BASE). The default (33) is usually correct, except for
        reads produced by some versions of the Illumina pipeline, where this should be set to 64.

    error_rate is the maximum allowed error rate (no. of errors divided by the length of the matching
        region). Default: 0.1

    min_qual_score allows you to trim low-quality ends from reads before adapter removal. The algorithm
        is the same as the one used by BWA (Subtract CUTOFF from all qualities; compute partial sums
        from all indices to the end of the sequence; cut sequence at the index at which the sum is minimal).

    min_read_len indicates the minimum length a read must be to be kept. Reads that are too short even
        before adapter removal are also discarded. In colorspace, an initial primer is not counted.

    keep_temp_files is a boolean indicating that you'd like to keep the intermediate files generated
        by this function. This can be useful for debugging, but in general should be left False.

    bowtie2 specifies whether to use the bowtie2 aligner instead of bowtie

    sort_mem is the parameter to pass to unix sort with -S/--buffer-size command
    """

    #Default bowtie option
    if aligner_options is None:
        if aligner.lower() == "minimap2":
            aligner_options = ["-ax", "sr","--secondary=no"]
        elif aligner.lower() == "bowtie":
            aligner_options = ["-S", "-k 1", "-m 1", "--chunkmbs 3072",
                               "--best", "--strata", "-o 4", "-e 80",
                               "-l 20", "-n 0"]
            aligner_options.append("--phred33-quals")
        else: # bowtie2
            aligner_options = []          
            aligner_options.append("--phred33-quals")

    # CASAVA >= 1.8
    quality_base = 33

    if len(path_to_output) != 0:
        path_to_output += "/"

    total_unique = 0

    file_name = sample+"_"+str(current_library)
    file_path = path_to_output+file_name

    print_checkpoint("Begin splitting reads for "+file_name)
    if pbat:
        total_input = split_fastq_file_pbat(num_procs, library_files, file_path+"_split_")
    else:
        total_input = split_fastq_file(num_procs, library_files, file_path+"_split_")

    if trim_reads:
        print_checkpoint("Begin trimming reads for "+file_name)
        quality_trim([file_path+"_split_"+str(i) for i in range(0, num_procs)],
                     output=[file_path+"_split_trimmed_"+str(i)
                             for i in range(0, num_procs)],
                     adapter_seq=adapter_seq,
                     error_rate=error_rate,
                     quality_base = quality_base,
                     min_qual_score=min_qual_score,
                     min_read_len=min_read_len,
                     input_format="fastq",
                     num_procs=num_procs,
                     max_adapter_removal=max_adapter_removal,
                     overlap_length=overlap_length,
                     zero_cap=zero_cap,
                     path_to_cutadapt=path_to_cutadapt)

        subprocess.check_call(shlex.split("rm "+" ".join([file_path+"_split_"+str(i)
                                                          for i in range(0,num_procs)])))

        print_checkpoint("Begin converting reads for "+file_name)
        if num_procs > 1:
            pool = multiprocessing.Pool(num_procs)
            for inputf, output in zip([file_path+"_split_trimmed_"+str(i) for i in range(0, num_procs)],
                                     [file_path+"_split_trimmed_converted_"+str(i)
                                      for i in range(0, num_procs)]):
                pool.apply_async(convert_reads,(inputf,output))
            pool.close()
            pool.join()
        else:
            for inputf, output in zip([file_path+"_split_trimmed_"+str(i) for i in range(0, num_procs)],
                                     [file_path+"_split_trimmed_converted_"+str(i)
                                      for i in range(0, num_procs)]):
                                          convert_reads(inputf,output)
        subprocess.check_call(shlex.split("rm "+
                                          " ".join([file_path+"_split_trimmed_"+str(i)
                                                    for i in range(0,num_procs)])))
        input_fastq = [file_path+"_split_trimmed_converted_"+str(i) for i in range(0, num_procs)]
    else:
        print_checkpoint("No trimming on reads")
        print_checkpoint("Begin converting reads for "+file_name)
        if num_procs > 1:
            pool = multiprocessing.Pool(num_procs)
            for inputf, output in zip([file_path+"_split_"+str(i) for i in range(0, num_procs)],
                                     [file_path+"_split_converted_"+str(i) for i in range(0, num_procs)]):
                pool.apply_async(convert_reads, (inputf, output))
            pool.close()
            pool.join()
        else:
            for inputf, output in zip([file_path+"_split_"+str(i) for i in range(0, num_procs)],
                                     [file_path+"_split_converted_"+str(i) for i in range(0, num_procs)]):
                convert_reads(inputf, output)
        subprocess.check_call(shlex.split("rm "+" ".join([file_path+"_split_"+str(i)
                                                          for i in range(0, num_procs)])))
        input_fastq = [file_path+"_split_converted_"+str(i) for i in range(0, num_procs)]

    #Run bowtie
    if aligner.lower() == "minimap2":
        print_checkpoint("Begin Running minimap2 for "+current_library)
    elif aligner.lower() == "Bowtie":
        print_checkpoint("Begin Running Bowtie for "+current_library)
    else:
        print_checkpoint("Begin Running Bowtie2 for "+current_library)
    total_unique = run_alignment(current_library,
                                 input_fastq,
                                 sample,
                                 forward_reference, reverse_reference, reference_fasta,
                                 path_to_output=path_to_output,
                                 aligner=aligner,
                                 aligner_options=aligner_options,
                                 merge_by_max_mapq=merge_by_max_mapq,
                                 min_mapq=min_mapq,
                                 path_to_aligner=path_to_aligner, num_procs=num_procs,
                                 keep_temp_files=keep_temp_files,
                                 sort_mem=sort_mem)

    #subprocess.check_call(shlex.split("rm " + " ".join(input_fastq)))

    return total_input, total_unique

def merge_bam_files(input_files,output,path_to_samtools=""):
    """
    This function will merge several bam files and create the correct header.
    
    input_files is a list of files produced by collapse_clonal reads. In other words, they're assumed
        to be named like <sample>_<processed_reads>_<lib_id>_no_clonal.bam
    
    output is the name of the merged bam file
    
    path_to_samtools is a string indicating the path to the directory containing your 
        installation of samtools. Samtools is assumed to be in your path if this is not
        provided.
    """ 
    f=open("header.sam",'w')
    subprocess.check_call(shlex.split(path_to_samtools+"samtools view -H "+input_files[0]),stdout=f)
    for filen in input_files:
        f.write("@RG\tID:" + filen[:filen.rindex(".bam")] + "\tLB:" +
                filen +
                "\tSM:NA" + "\n")
    f.close()
    subprocess.check_call(shlex.split(path_to_samtools+"samtools merge -r -h header.sam "+ output +" "+" ".join(input_files)))
    subprocess.check_call(["rm", "header.sam"])

def build_ref(input_files,
              output,
              buffsize=100,
              aligner="bowtie2",
              path_to_aligner="",
              num_procs=1):
    """
    Creates 2 reference files: one with all C's converted to T's, and one with all G's converted to A's
    
    input_files is a list of files to build a reference from
    
    output is the prefix of the two output reference files that will be created
    
    buffsize is the number of bytes that will be read in from the reference at once    
    """
    if len(path_to_aligner) !=0:
        path_to_aligner+="/"

    if not isinstance(input_files, list):
        if isinstance(input_files, str):
            input_files = [input_files]
        else:
            sys.exit("input_files must be a list of strings")
    #outf = Convert all C to T. outr = convert all G to A  
    with open(output+"_f.fasta", 'w') as outf, open(output+"_r.fasta", 'w') as outr:
        for filen in input_files:
            f = open(filen, 'r')
            line = f.read(buffsize)
            header = False #indicates when you are currently reading in a header
            while line:
                for base in line:
                    if header==True:
                        if base=="\n": #when you encounter a newline, you are no longer in a header
                            header=False
                        outf.write(base)
                        outr.write(base)
                    else:
                        if base==">": #all headers begin with >
                            header=True
                            outf.write(base)
                            outr.write(base)
                        elif base=="C":
                            outf.write("T")
                            outr.write(base)
                        elif base=="c":
                            outf.write("t")
                            outr.write(base)
                        elif base=="G":
                            outf.write(base)
                            outr.write("A")
                        elif base=="g":
                            outf.write(base)
                            outr.write("a")
                        else:
                            outf.write(base)
                            outr.write(base)
                line = f.read(buffsize)
            f.close()

    # minimap2
    if aligner.lower() == "minimap2":
        subprocess.check_call([path_to_aligner+"minimap2",
                               "-t",str(num_procs),
                               "-d",output+"_f.mmi",
                               output + "_f.fasta"])
        subprocess.check_call([path_to_aligner+"minimap2",
                               "-t",str(num_procs),
                               "-d",output+"_r.mmi",
                               output + "_r.fasta"])
        return 0
        
    # bowtie2
    base_cmd = path_to_aligner+"bowtie2-build -f "
    # bowtie
    if aligner.lower() == "bowtie":
        base_cmd = path_to_aligner+"bowtie-build -f "

    if num_procs > 1:
        pool = multiprocessing.Pool(2)
        pool.apply_async(subprocess.check_call,(shlex.split(base_cmd + output + "_f.fasta " + output +"_f"),))
        pool.apply_async(subprocess.check_call,(shlex.split(base_cmd + output + "_r.fasta " + output+ "_r"),))
        pool.close()
        pool.join()
    else:
        subprocess.check_call(shlex.split(base_cmd + output + "_f.fasta " + output +"_f"))
        subprocess.check_call(shlex.split(base_cmd + output + "_r.fasta " + output+ "_r"))
    subprocess.check_call(["rm",output + "_f.fasta",output + "_r.fasta"])
    return 0
    
def convert_reads(inputf,output,buffer_line_number=100000):
    """
    This function takes a fastq file as input and converts all the cytosines in reads to thymines for
    mapping to bisulfite converted genomes. This function also stores an encoding of where the cytosines
    were located in the header of each fastq read. See encode_c_positions for more detail.
    
    input is a fastq file for conversion
    
    output is the name of the file you'd like to put the converted reads in
    """
    f = open(inputf,'r')
    g = open(output,'w')
    header = f.readline().rstrip()
    header = header.replace(" ","!")
    seq = f.readline()
    header2 = f.readline()
    qual = f.readline()
    encoding = encode_c_positions(seq)

    line_counts = 0
    out = ""
    while header:
        out += header+"!"+encoding+"\n"
        converted_seq = seq.replace("C","T")
        out += converted_seq
        out += header2
        out += qual
        line_counts += 4
        # output
        if line_counts > buffer_line_number:
            g.write(out)
            line_counts = 0
            out = ""
        # update
        header = f.readline().rstrip()
        header = header.replace(" ","!")
        seq = f.readline()
        header2 = f.readline()
        qual = f.readline()
        encoding = encode_c_positions(seq)
    # output
    if line_counts > 0:
        g.write(out)
        line_counts = 0
        out = ""

    f.close()
    g.close()
    
def encode_c_positions(seq,is_read2=False):
    """
    This function creates an encoding of where cytosine nucleotides are located in a converted read.
    The encoding uses ascii characters (minus an offset) to indicate an offset into the read. 
    For example, the ascii character # has an integer value of 36 and indicates that a C is located 
    2 bases from the previous position (36 - 34). The offsets build off of one another so if the first
    offset is 2 and the second offset is 5 the second C is located in the 9th position (since python indexing
    starts at 0). In other words, next_c_index = prev_c_index + offset + 1.
    
    seq is a string of nucleotides you'd like to encode.
    """
    indexes = ""
    prev_index = 0
    if is_read2==False:        
        index = seq.find("C",prev_index)
        offset = index + 34
        while True:
            if index < 0:
                break
            while offset >= 255:
                indexes += chr(255)
                offset -= 255
            if offset < 34:
                offset += 34
            indexes += chr(offset)
            prev_index = index + 1
            index = seq.find("C",prev_index)
            offset = index - prev_index + 34
    else:
        index = seq.find("G",prev_index)
        offset = index + 34
        while True:
            if index < 0:
                break
            while offset >= 255:
                indexes += chr(255)
                offset -= 255
            if offset < 34:
                offset += 34
            indexes += chr(offset)                
            prev_index = index + 1
            index = seq.find("G",prev_index)
            offset = index - prev_index + 34
    return indexes

def decode_c_positions(seq,indexes,strand,is_read2=False):
    """
    This function takes the encodings generated by encode_c_position and replaces the appropriate
    positions with C nucleotides.
    
    seq is a string of nucleotides to have Cs or Gs replaced.
    
    indexes is a string of characters indicating the offsets for the positions of the Cs or Gs.
    
    strand is the DNA strand (+ or -) that seq mapped to. This is important because
        sequences in sam files are always represented on the forward strand

    is_read2 indicates whether the read to be deconverted is a read2.
    """

    prev_index = 0
    new_seq=""
    index = 0
    saturated = False
    if is_read2 == False:
        if strand == "-":
            seq = seq[::-1]
        for char in indexes:
            offset = ord(char)
            if offset == 255:
                index += 255
                saturated = True
                continue
            if saturated and (offset - 34) < 34:
                offset -= 68
            else:
                offset -= 34
            index += offset
            if strand == "+":
                new_seq += seq[prev_index:index]+"C"
            elif strand == "-":
                new_seq += seq[prev_index:index]+"G"
            prev_index = index + 1
            index = prev_index
            saturated = False
    else:
        if strand == "-":
            seq = seq[::-1]
        for char in indexes:
            offset = ord(char)
            if offset == 255:
                index += 255
                saturated = True
                continue
            if saturated and (offset - 34) < 34:
                offset -= 68
            else:
                offset -= 34
            index += offset
            if strand == "+":
                new_seq += seq[prev_index:index]+"G"
            elif strand == "-":
                new_seq += seq[prev_index:index]+"C"
            prev_index = index + 1
            index = prev_index
            saturated = False
    new_seq += seq[prev_index:]
    if strand == "-":
        new_seq = new_seq[::-1]
    return new_seq

def run_alignment(current_library,library_read_files,
                  sample,
                  forward_reference,reverse_reference,reference_fasta,
                  path_to_output="",
                  path_to_samtools="",
                  aligner="bowtie2",
                  path_to_aligner="",aligner_options=None,
                  merge_by_max_mapq=False,min_mapq=30,
                  num_procs=1,keep_temp_files=False,sort_mem="500M"):
    """
    This function runs bowtie on the forward and reverse converted bisulfite references 
    (generated by build_ref). It removes any read that maps to both the forward and reverse
    strands.
    
    files is a list of file paths to be mapped
    
    forward_reference is a string indicating the path to the forward strand reference created by
        build_ref
    
    reverse_reference is a string indicating the path to the reverse strand reference created by
        build_ref
    
    prefix is a string that you would like prepended to the output files (e.g., the sample name)
    
    options is a list of strings indicating options you'd like passed to bowtie 
        (e.g., ["-k 1","-l 2"]
    
    path_to_aligner is a string indicating the path to the folder in which bowtie resides. Bowtie
        is assumed to be in your path if this option isn't used
    
    num_procs is an integer indicating the number of processors you'd like used for removing multi
        mapping reads and for bowtie mapping
    
    keep_temp_files is a boolean indicating that you'd like to keep the intermediate files generated
        by this function. This can be useful for debugging, but in general should be left False.
    
    bowtie2 specifies whether to use the bowtie2 aligner instead of bowtie
    
    sort_mem is the parameter to pass to unix sort with -S/--buffer-size command
    """

    if not sort_mem:
        sort_option = ""
    else:
        sort_option = " -S "+sort_mem

    if len(path_to_aligner) !=0:
        path_to_aligner+="/"

    if len(path_to_output) !=0:
        path_to_output+="/"

    prefix = path_to_output+sample+"_"+str(current_library)

    #Default bowtie option
    if aligner_options is None:
        if aligner.lower() == "minimap2":
            aligner_options = ["-ax", "sr","--secondary=no"]
        elif aligner.lower() == "bowtie":
            aligner_options = ["-S", "-k 1", "-m 1", "--chunkmbs 3072",
                               "--best", "--strata", "-o 4", "-e 80",
                               "-l 20", "-n 0"]
        else: #bowtie 2
            aligner_options = []

    options = aligner_options

    input_read_file = ",".join(library_read_files)
    input_read_file = prefix+"_converted_reads.fastq"
    subprocess.check_call(["mv",library_read_files[0],input_read_file])
    if len(library_read_files) > 1:
        with open(input_read_file,'a') as g:
            for library_read_file in library_read_files[1:]:
                with open(library_read_file,'r') as f:
                    g.write(f.read())
                subprocess.check_call(["rm",library_read_file])

    if aligner != "minimap2":
        if " ".join(options).find(" -p ") == -1:
            options.append("-p "+str(num_procs))
    else:
        if " ".join(options).find(" -t ") == -1:
            options.append("-t "+str(num_procs))

    if aligner.lower() == "minimap2":
        args = [path_to_aligner+"minimap2"]
        args.extend(options)
        args.append("--for-only")
        args.append(forward_reference)
        args.append(input_read_file)
    elif aligner.lower() == "bowtie":
        args = [path_to_aligner+"bowtie"]
        args.extend(options)
        args.append("--norc")
        args.append(forward_reference)
        args.append(input_read_file)
    else: # bowtie2
        args = [path_to_aligner+"bowtie2"]
        args.extend(options)
        args.append("--norc")
        args.append("-x "+forward_reference)
        args.append("-U "+input_read_file)
    ## run
    with open(prefix+"_forward_strand_hits.sam","w") as f:
        subprocess.check_call(shlex.split(" ".join(args)),stdout=f)
    print_checkpoint("Processing forward strand hits")
    find_multi_mappers(prefix+"_forward_strand_hits.sam",
                       prefix,
                       num_procs=num_procs,
                       min_mapq=min_mapq,
                       append=False,
                       keep_temp_files=keep_temp_files)

    if aligner.lower() == "minimap2":
        args = [path_to_aligner+"minimap2"]
        args.extend(options)
        args.append("--rev-only")
        args.append(reverse_reference)
        args.append(input_read_file)
    elif aligner.lower() == "bowtie":
        args = [path_to_aligner+"bowtie"]
        args.extend(options)
        args.append("--nofw")
        args.append(reverse_reference)
        args.append(input_read_file)
    else:
        args = [path_to_aligner+"bowtie2"]
        args.extend(options)
        args.append("--nofw")
        args.append("-x "+reverse_reference)
        args.append("-U "+input_read_file)
    ## run
    with open(prefix+"_reverse_strand_hits.sam","w") as f:
        subprocess.check_call(shlex.split(" ".join(args)),stdout=f)
    subprocess.check_call(["rm",input_read_file])        
    print_checkpoint("Processing reverse strand hits")
    sam_header = find_multi_mappers(prefix+"_reverse_strand_hits.sam",
                                    prefix,
                                    num_procs=num_procs,
                                    min_mapq=min_mapq,
                                    append=True,
                                    keep_temp_files=keep_temp_files)
    ## Clear temporary files
    if num_procs > 1:
        pool = multiprocessing.Pool(num_procs)
        for file_num in range(0,num_procs):
            pool.apply_async(subprocess.check_call,
                             (shlex.split(
                                 "env LC_COLLATE=C sort"+sort_option+ \
                                 " -t '\t' -k 1 -o "+prefix+"_sorted_"+str(file_num)+ \
                                 " "+prefix+"_sorted_"+str(file_num)),))
        pool.close()
        pool.join()
    else:
        for file_num in range(0,num_procs):
            subprocess.check_call(shlex.split(
                "env LC_COLLATE=C sort"+sort_option + " -t '\t' -k 1 -o "+ \
                prefix+"_sorted_"+str(file_num)+" "+prefix+"_sorted_"+str(file_num)))

    print_checkpoint("Finding multimappers")

    if merge_by_max_mapq:
        total_unique = merge_sorted_multimap_max_mapq(
            current_library,
            [prefix+"_sorted_"+str(file_num) for file_num in range(0,num_procs)],
            prefix,
            reference_fasta,
            path_to_samtools="")
    else:
        total_unique = merge_sorted_multimap(
            current_library,
            [prefix+"_sorted_"+str(file_num) for file_num in range(0,num_procs)],
            prefix,
            reference_fasta,
            path_to_samtools="")

    subprocess.check_call(shlex.split("rm "+" ".join([prefix+"_sorted_"+str(file_num)
                                                      for file_num in range(0,num_procs)])))

    output_bam_file = prefix+"_processed_reads.bam"
    if not sort_mem:
        sort_option = ""
    else:
        sort_option = " -m "+sort_mem

    try:
        subprocess.check_call(shlex.split(path_to_samtools+"samtools sort "+
                                          " -@ " + str(num_procs) +
                                          sort_option + " " +
                                          " -o "+output_bam_file + " " +
                                          output_bam_file ))
    except:
        subprocess.check_call(shlex.split(path_to_samtools+"samtools sort "+
                                          " -o "+output_bam_file + " " +
                                          output_bam_file ))
    return total_unique
 
def find_multi_mappers(inputf,output,num_procs=1,min_mapq=30,
                       keep_temp_files=False,append=False):
    """
    This function takes a sam file generated by bowtie and pulls out any mapped reads.
    It splits these mapped reads into num_procs number of files.
    
    inputf is a string of the path to a sam file from bowtie
    
    output is a string of the prefix you'd like prepended to the output files
        The output files will be named as <output>_sorted_<index num>
    
    num_procs is an integer indicating how many files the bowtie sam file should be split
        into
    
    keep_temp_files is a boolean indicating that you'd like to keep the intermediate files generated
        by this function. This can be useful for debugging, but in general should be left False.
    
    append is a boolean that should be False for the first bowtie sam file you process (i.e., for the forward
        mapped reads) and True for the second. This option is mainly for safety. It ensures that files from
        previous runs are erased.
    """
    min_mapq = max(3,min_mapq)
    sam_header = []
    file_handles = {}
    f = open(inputf,'r')
    cycle = itertools.cycle(list(range(0,num_procs)))
    for file_num in range(0,num_procs):
        if append == False:
            file_handles[file_num]=open(output+"_sorted_"+str(file_num),'w')
        else:
            file_handles[file_num]=open(output+"_sorted_"+str(file_num),'a')
    for line in f:
        #To deal with the way chromosomes were named in some of our older references
        if line[0] == "@":
            continue

        fields = line.split("\t")
        flag = int(fields[1])
        # minimum QC
        if fields[2] == "*" or int(fields[4]) < min_mapq or (flag & 2048 == 2048):
            continue
        header = fields[0].split("!")
        #BIG ASSUMPTION!! NO TABS IN FASTQ HEADER LINES EXCEPT THE ONES I ADD!
        if (flag & 16) == 16:
            strand = "-"
        elif (flag & 16) == 0:
            strand = "+"
        try:
            seq = decode_c_positions(fields[9],header[-1],strand)
            file_handles[next(cycle)].write(" ".join(header[:-1])+"\t"+"\t".join(fields[1:9])
                                            +"\t"+seq+"\t"+"\t".join(fields[10:]))
        except:
            print_warning("  Failed to recover unconverted sequence for:\n"+line+"\n")
            print_warning(header[-1]+"\n")
    f.close()

    subprocess.check_call(shlex.split("rm "+inputf))
    for file_num in range(0,num_procs):
        file_handles[file_num].close()
    
def merge_sorted_multimap(current_library,files,prefix,reference_fasta,path_to_samtools=""):
    """
    This function takes the files from find_multi_mappers and outputs the uniquely mapping reads.
    
    files is a list of filenames containing the output of find_multi_mappers
    
    output is a prefix you'd like prepended to the bam file containing the uniquely mapping reads
        This file will be named as <output>+"_no_multimap_"+<index_num>
    """

    output_sam_file = prefix+"_processed_reads.sam"
    output_bam_file = prefix+"_processed_reads.bam"
    output_handle = open(output_sam_file,'w')

    #output_pipe = subprocess.Popen(
    #    shlex.split(path_to_samtools+"samtools view -S -b -"),
    #    stdin=subprocess.PIPE,stdout=output_handle)

    try:
        f = open(reference_fasta+".fai",'r')
    except:
        print("Reference fasta not indexed. Indexing.")
        try:
            subprocess.check_call(shlex.split(path_to_samtools+"samtools faidx "+reference_fasta))
            f = open(reference_fasta+".fai",'r')
        except:
            sys.exit("Reference fasta wasn't indexed, and couldn't be indexed. Please try indexing it manually and running methylpy again.")
    #Create sam header based on reference genome
    output_handle.write("@HD\tVN:1.0\tSO:unsorted\n")
    for line in f:
        fields = line.split("\t")
        output_handle.write("@SQ\tSN:"+fields[0]+"\tLN:"+fields[1]+"\n")
    f.close()

    ## Merging alignment results of both strands    
    lines = {}
    fields = {}
    file_handles = {}
    total_unique = 0
    count= 0    
    for index,filen in enumerate(files):
        file_handles[filen]=open(filen,'r')
        lines[filen]=file_handles[filen].readline()
        fields[filen] = lines[filen].split("\t")[0]
    while True:
        all_fields = [field for field in list(fields.values()) if field != ""]
        if len(all_fields) == 0:
            break
        min_field = min(all_fields)
        count = 0
        current_line = ""
        current_field = ""
        for key in fields:
            while fields[key] == min_field:
                count += 1
                current_line = lines[key]
                lines[key]=file_handles[key].readline()
                fields[key]=lines[key].split("\t")[0]
        if count == 1:
            output_handle.write(current_line)
            total_unique += 1
            
    #output_pipe.stdin.close()
    output_handle.close()
    
    for index,filen in enumerate(files):
        file_handles[filen].close()

    f = open(output_bam_file,'w')
    subprocess.check_call(shlex.split(path_to_samtools+"samtools view -S -b -h "+output_sam_file),stdout=f)
    f.close()
    subprocess.check_call(shlex.split("rm "+output_sam_file))
    return total_unique

def merge_sorted_multimap_max_mapq(current_library,files,prefix,reference_fasta,path_to_samtools=""):
    """
    This function takes the files from find_multi_mappers and outputs the uniquely mapping reads.
    
    files is a list of filenames containing the output of find_multi_mappers
    
    output is a prefix you'd like prepended to the bam file containing the uniquely mapping reads
        This file will be named as <output>+"_no_multimap_"+<index_num>
    """

    output_sam_file = prefix+"_processed_reads.sam"
    output_bam_file = prefix+"_processed_reads.bam"
    output_handle = open(output_sam_file,'w')

    #output_pipe = subprocess.Popen(
    #    shlex.split(path_to_samtools+"samtools view -S -b -"),
    #    stdin=subprocess.PIPE,stdout=output_handle)

    try:
        f = open(reference_fasta+".fai",'r')
    except:
        print("Reference fasta not indexed. Indexing.")
        try:
            subprocess.check_call(shlex.split(path_to_samtools+"samtools faidx "+reference_fasta))
            f = open(reference_fasta+".fai",'r')
        except:
            sys.exit("Reference fasta wasn't indexed, and couldn't be indexed. Please try indexing it manually and running methylpy again.")
    #Create sam header based on reference genome
    output_handle.write("@HD\tVN:1.0\tSO:unsorted\n")
    for line in f:
        fields = line.split("\t")
        output_handle.write("@SQ\tSN:"+fields[0]+"\tLN:"+fields[1]+"\n")
    f.close()

    ## Merging alignment results of both strands    
    lines = {}
    fields = {}
    file_handles = {}
    total_unique = 0
    count= 0    
    for index,filen in enumerate(files):
        file_handles[filen]=open(filen,'r')
        lines[filen]=file_handles[filen].readline()
        fields[filen] = lines[filen].split("\t")[0]
    while True:
        all_fields = [field for field in list(fields.values()) if field != ""]
        if len(all_fields) == 0:
            break
        min_field = min(all_fields)
        count, count_diff_mapq = 0, 0
        current_line = ""
        current_field = ""
        max_mapq = -100
        for key in fields:
            while fields[key] == min_field:
                mapq = int(lines[key].split("\t")[4])
                count += 1
                if mapq > max_mapq:
                    count_diff_mapq += 1
                    max_mapq = mapq
                    current_line = lines[key]
                lines[key]=file_handles[key].readline()
                fields[key]=lines[key].split("\t")[0]
        if count == 1 or count_diff_mapq > 1:
            output_handle.write(current_line)
            total_unique += 1

    #output_pipe.stdin.close()
    output_handle.close()
    
    for index,filen in enumerate(files):
        file_handles[filen].close()

    f = open(output_bam_file,'w')
    subprocess.check_call(shlex.split(path_to_samtools+"samtools view -S -b -h "+output_sam_file),stdout=f)
    f.close()

    subprocess.check_call(shlex.split("rm "+output_sam_file))
    return total_unique

def quality_trim(inputf, output = None, quality_base = None, min_qual_score = None, min_read_len = None, 
                 adapter_seq = "AGATCGGAAGAGCACACGTCTG", num_procs = 1, input_format = None, error_rate = None, 
                 max_adapter_removal = None, overlap_length = None, zero_cap = False, path_to_cutadapt = ""):
    """
    Information from cutadapt documentation:
    input_format:
        Input file format; can be either 'fasta', 'fastq' or 'sra-fastq'. Ignored when reading csfasta/qual files
        (default: auto-detect from file name extension).

    adapter_seq:
        Sequence of an adapter that was ligated to the 3' end. The adapter itself and anything that follows is
        trimmed.
    
    error_rate:
        Maximum allowed error rate (no. of errors divided by the length of the matching region) (default: 0.1)

    max_adapter_removal:
        Try to remove adapters at most COUNT times. Useful when an adapter gets appended multiple times.
        
    overlap_length:
        Minimum overlap length. If the overlap between the read and the adapter is shorter than LENGTH, the read
        is not modified.This reduces the no. of bases trimmed purely due to short random adapter matches.

    min_read_len:
        Discard trimmed reads that are shorter than LENGTH. Reads that are too short even before adapter removal
        are also discarded. In colorspace, an initial primer is not counted.

    output:
        Write the modified sequences to this file instead of standard output and send the summary report to
        standard output. The format is FASTQ if qualities are available, FASTA otherwise.

    min_qual_score:
        Trim low-quality ends from reads before adapter removal. The algorithm is the same as the one used by
        BWA (Subtract CUTOFF from all qualities; compute partial sums from all indices to the end of the
        sequence; cut sequence at the index at which the sum is minimal).
        
    quality_base:
        Assume that quality values are encoded as ascii(quality + QUALITY_BASE). The default (33) is
        usually correct, except for reads produced by some versions of the Illumina pipeline, where this should
        be set to 64.
    
    zero_cap:
        Change negative quality values to zero (workaround to avoid segmentation faults in BWA).
    
    path_to_cutadapt:
        Path to the folder where cutadapt executable exists. If none, assumes it can be run from current directory
        
    input:
        list of filenames
    """
    if path_to_cutadapt:  #see if cutadapt is installed
        if path_to_cutadapt[-1]!="/":
            path_to_cutadapt += "/"
    path_to_cutadapt += "cutadapt"
    try:
        devnull = open('/dev/null', 'w')
        subprocess.check_call([path_to_cutadapt], stdout=devnull, stderr=devnull)
    except OSError:
        sys.exit("Cutadapt must be installed to run quality_trim")
    except:
        devnull.close()
                 
    if not isinstance(inputf, list):
        if isinstance(inputf, str):
            inputf = [inputf]
        else:
            sys.exit("input must be a list of strings")
    if not isinstance(output, list):
        if isinstance(output, str):
            output = [output]
        else:
            sys.exit("output must be a list of strings")
            
    if len(output) != len(inputf):
        sys.exit("Must provide an equal number of input and output files")
    base_cmd = path_to_cutadapt
    options = " --quiet "
    if zero_cap:
        zero = "-z "
    else:
        zero = ""  
    
    if input_format:
        options += " -f " + input_format
    if error_rate:
        options += " -e " + str(error_rate)
    if max_adapter_removal:
        options += " -n " + str(max_adapter_removal)
    if overlap_length:
        options += " -O " + str(overlap_length)
    if min_read_len:
        options += " -m " + str(min_read_len)
    if min_qual_score:
        options += " -q " + str(min_qual_score)
    if quality_base:
        options += " --quality-base=" + str(quality_base)
    options += " -a " + adapter_seq
    options += " " + zero
    if num_procs > 1:
        pool = multiprocessing.Pool(num_procs)
        #adapter trimming
        for current_input,current_output in zip(inputf,output):
            if output:
                options += " -o " + current_output + " "
            pool.apply_async(subprocess.check_call,(base_cmd + options + current_input,),{"shell":True})
        pool.close()
        pool.join()
    else:
        for current_input,current_output in zip(inputf,output):
            if output:
                options += " -o " + current_output + " "
            subprocess.check_call(base_cmd + options + current_input, shell=True)
        

def remove_clonal_bam(input_bam,output_bam,metric,is_pe=False,path_to_picard="",java_options="-Xmx20g"):
    """
    Running picard to remove clonal reads in input_bam and output non-clonal reads
    to output_bam.
    """
    
    subprocess.check_call(
        shlex.split(
            " ".join(["java",java_options,
                      "-jar",
                      path_to_picard+"/picard.jar MarkDuplicates",
                      "INPUT="+input_bam,
                      "OUTPUT="+output_bam,
                      "ASSUME_SORTED=true",
                      "REMOVE_DUPLICATES=true",
                      "METRICS_FILE="+metric,
                      "VALIDATION_STRINGENCY=LENIENT",
                      "QUIET=true"])
        )
    )
    total_clonal = 0
    with open(metric,'r') as f:
        while True:
            line = f.readline()
            if line[0] != "#" and len(line) != 1:
                break
        line = f.readline()
        fields = line.split("\t")
        if is_pe:
            total_clonal = fields[6]
        else:
            total_clonal = fields[5]
    return int(total_clonal)

    
def fasta_iter(fasta_name,query_chrom):
    """
    given a fasta file. yield tuples of header, sequence
    """
    fh = open(fasta_name)
    # ditch the boolean (x[0]) and just keep the header or sequence since
    # we know they alternate.
    from itertools import groupby
    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
    seq = None
    for header in faiter:
        # drop the ">"
        header = header.next()[1:].strip()
        header = header.split(" ")[0]
        # join all sequence lines to one.
        if header != query_chrom:
            for s in next(faiter):
                pass
            continue
        seq = "".join(s.strip() for s in next(faiter))
        return seq

def get_chromosome_sequence(fasta_name,query_chrom):
    chrom_pointer = None
    with open(fasta_name+".fai",'r') as f:
        for line in f:
            fields = line.split("\t")
            if fields[0] == query_chrom:
                chrom_pointer = int(fields[2])
    if chrom_pointer is None: return(None)

    seq = ""
    with open(fasta_name,'r') as f:
        f.seek(chrom_pointer)
        for line in f:
            if line[0] == ">": break
            seq += line.rstrip("\n")
    return(seq)
                    
def call_methylated_sites(inputf, sample, reference_fasta,
                          unmethylated_control=None,
                          sig_cutoff=.01,num_procs = 1,
                          num_upstr_bases=0,num_downstr_bases=2,
                          generate_mpileup_file=True,
                          compress_output=True,
                          bgzip=False,
                          path_to_bgzip="",
                          path_to_tabix="",
                          buffer_line_number = 100000,
                          min_mapq=30,
                          min_cov=1,binom_test=True,
                          path_to_samtools="",
                          remove_chr_prefix=True,
                          sort_mem="500M",
                          add_snp_info=False,
                          path_to_files="",min_base_quality=1,
                          keep_temp_files=False):

    """
    inputf is the path to a bam file that contains mapped bisulfite sequencing reads
    
    sample is the name you'd like for the allc files. The files will be named like so:
        allc_<sample>_<chrom>.tsv
    
    reference is the path to a samtools indexed fasta file
    
    control is the name of the chromosome/region that you want to use to estimate the non-conversion rate of your 
        sample, or the non-conversion rate you'd like to use. Consequently, control is either a string, or a decimal
        If control is a string then it should be in the following format: "chrom:start-end". 
        If you'd like to specify an entire chromosome simply use "chrom:"
    
    sig_cutoff is a float indicating the adjusted p-value cutoff you wish to use for determining whether or not
        a site is methylated
    
    num_procs is an integer indicating how many num_procs you'd like to run this function over
    
    min_cov is an integer indicating the minimum number of reads for a site to be tested.
    
    path_to_files is a string indicating the path for the output and the input bam, mpileup, or allc files
        for methylation calling.
    min_base_quality is an integer indicating the minimum PHRED quality score for a base to be included in the
        mpileup file (and subsequently to be considered for methylation calling)
    """

    if add_snp_info:
        return call_methylated_sites_with_SNP_info(inputf, sample, reference_fasta,
                                                   unmethylated_control=unmethylated_control,
                                                   sig_cutoff=sig_cutoff,
                                                   num_procs=num_procs,
                                                   num_upstr_bases=num_upstr_bases,
                                                   num_downstr_bases=num_downstr_bases,
                                                   generate_mpileup_file=generate_mpileup_file,
                                                   compress_output=compress_output,
                                                   buffer_line_number=buffer_line_number,
                                                   min_mapq=min_mapq,
                                                   min_cov=min_cov,
                                                   binom_test=binom_test,
                                                   path_to_samtools=path_to_samtools,
                                                   remove_chr_prefix=remove_chr_prefix,
                                                   sort_mem=sort_mem,
                                                   path_to_files=path_to_files,
                                                   min_base_quality=min_base_quality,
                                                   keep_temp_files=keep_temp_files)

    
    if binom_test and unmethylated_control is None:
        print_error("Please specify unmethylated_control if you would like to do binomial test!\n")

    #Figure out all the correct quality options based on the offset or CASAVA version given
    # quality_version >= 1.8:
    quality_base = 33
    if len(path_to_files)!=0:
        path_to_files+="/"
    if len(path_to_samtools)!=0:
        path_to_samtools+="/"

    try:
        num_procs = int(num_procs)
    except:
        sys.exit("num_procs must be an integer")
        
    try:
        #make sure bam file is indexed
        open(inputf+".bai",'r')
    except:
        print_checkpoint("Input not indexed. Indexing...")
        subprocess.check_call(shlex.split(path_to_samtools+"samtools index "+inputf))

    ## Check fasta index
    try:
        f = open(reference_fasta+".fai",'r')
    except:
        print("Reference fasta not indexed. Indexing.")
        try:
            subprocess.check_call(shlex.split(path_to_samtools+"samtools faidx "+reference_fasta))
            f = open(reference_fasta+".fai",'r')
        except:
            sys.exit("Reference fasta wasn't indexed, and couldn't be indexed. "
                     +"Please try indexing it manually and running methylpy again.")

    ## Input
    if not generate_mpileup_file:
        cmd = path_to_samtools+"samtools mpileup -Q "+str(min_base_quality)+\
              +" -q "+str(min_mapq)+" -B -f "+reference_fasta+" "+inputf
        pipes = subprocess.Popen(shlex.split(cmd),
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE,
                                 universal_newlines=True)
        fhandle = pipes.stdout
    else:
        with open(path_to_files+sample+"_mpileup_output.tsv",'w') as f:
            subprocess.check_call(
                shlex.split(
                    path_to_samtools+"samtools mpileup -Q "+str(min_base_quality)
                    +" -q "+str(min_mapq)
                    +" -B -f "+reference_fasta+" "+inputf),
                stdout=f)
        fhandle = open(path_to_files+sample+"_mpileup_output.tsv" ,'r')

    ## Output
    if compress_output:
        output_filehandler = gzip.open(path_to_files+"allc_"+sample+".tsv.gz",'wt')
        output_file = path_to_files+"allc_"+sample+".tsv.gz"
    else:
        output_filehandler = open(path_to_files+"allc_"+sample+".tsv",'w')
        output_file = path_to_files+"allc_"+sample+".tsv"
        
    complement = {"A":"T","C":"G","G":"C","T":"A","N":"N"}
    context_len = num_upstr_bases+1+num_downstr_bases
    cur_chrom = ""
    #cur_chrom_nochr = ""
    line_counts = 0
    out = ""
    for line in fhandle:
        fields = line.split("\t")
        if fields[0] != cur_chrom:
            cur_chrom = fields[0]
            cur_chrom_nochr = cur_chrom
            if remove_chr_prefix and cur_chrom.startswith("chr"):
                cur_chrom_nochr = cur_chrom_nochr[3:]
            seq = get_chromosome_sequence(reference_fasta,cur_chrom)
            if seq != None:
                seq = seq.upper()

        if seq == None:
            continue
        if (not fields[2] == "C") and (not fields[2] == "G"):
            continue

        # indels
        read_bases = fields[4]
        incons_basecalls = read_bases.count("+") + read_bases.count("-")
        if incons_basecalls > 0:
            read_bases_no_indel = ""
            index = 0
            prev_index = 0
            while index < len(read_bases):
                if read_bases[index] == "+" or read_bases[index] == "-":
                    # get insert size
                    indel_size = ""
                    ind = index+1
                    while True:
                        try:
                            int(read_bases[ind])
                            indel_size += read_bases[ind]
                            ind += 1
                        except:                    
                            break
                    try:
                        # sometimes +/- does not follow by a number and
                        # it should be ignored
                        indel_size = int(indel_size)
                    except:
                        index += 1
                        continue
                    read_bases_no_indel += read_bases[prev_index:index]
                    index = ind + indel_size
                    prev_index = index
                else:
                    index += 1
            read_bases_no_indel += read_bases[prev_index:index]
            fields[4] = read_bases_no_indel

        # count converted and unconverted bases
        if fields[2] == "C":            
            pos = int(fields[1])-1
            try:
                context = seq[(pos-num_upstr_bases):(pos+num_downstr_bases+1)]
            except: # complete context is not available, skip
                continue
            unconverted_c = fields[4].count(".")
            converted_c = fields[4].count("T")
            cov = unconverted_c+converted_c
            if cov > 0 and len(context) == context_len:
                line_counts += 1
                out += "\t".join([cur_chrom_nochr,str(pos+1),"+",context,
                                  str(unconverted_c),str(cov),"1"])+"\n"
        elif fields[2] == "G":
            pos = int(fields[1])-1
            try:
                context = "".join([complement[base]
                                   for base in reversed(
                                           seq[(pos-num_downstr_bases):(pos+num_upstr_bases+1)]
                                   )]
                )
            except: # complete context is not available, skip
                continue
            unconverted_c = fields[4].count(",")
            converted_c = fields[4].count("a")
            cov = unconverted_c+converted_c
            if cov > 0 and len(context) == context_len:
                line_counts += 1
                out += "\t".join([cur_chrom_nochr,str(pos+1),"-",context,
                                  str(unconverted_c),str(cov),"1"])+"\n"
        if line_counts > buffer_line_number:
            output_filehandler.write(out)
            line_counts = 0
            out = ""

    if line_counts > 0:
        output_filehandler.write(out)
        line_counts = 0
        out = ""
    fhandle.close()
    output_filehandler.close()

    if generate_mpileup_file and not keep_temp_files:
        subprocess.check_call(shlex.split("rm -f "+path_to_files+sample+"_mpileup_output.tsv"))

    if binom_test:
        print_checkpoint('Perform binomial test')
        perform_binomial_test(allc_file=output_file,
                              sample=sample,
                              path_to_output=path_to_files,
                              unmethylated_control=unmethylated_control,
                              min_cov=min_cov,
                              sig_cutoff=sig_cutoff,
                              num_procs=num_procs,
                              sort_mem=sort_mem,
                              compress_output=compress_output,
                              buffer_line_number=buffer_line_number,
                              remove_chr_prefix=remove_chr_prefix)
    elif not unmethylated_control is None:
        non_conversion = calculate_non_conversion_rate(unmethylated_control,
                                                       output_file,
                                                       chrom_pointer=None,
                                                       remove_chr_prefix=remove_chr_prefix)
    index_allc_file(output_file)
    if bgzip:
        bgzip_allc_file(output_file,path_to_bgzip,path_to_tabix,buffer_line_number)
    return(0)


def analyze_read_basecalls(ref,read_bases):
    # indels
    incons_basecalls = read_bases.count("+") + read_bases.count("-")
    if incons_basecalls > 0:
        read_bases_no_indel = ""
        index = 0
        prev_index = 0
        while index < len(read_bases):
            if read_bases[index] == "+" or read_bases[index] == "-":
                # get insert size
                indel_size = ""
                ind = index+1
                while True:
                    try:
                        int(read_bases[ind])
                        indel_size += read_bases[ind]
                        ind += 1
                    except:                    
                        break
                try:
                    indel_size = int(indel_size)
                except:
                    indel_size = 0
                read_bases_no_indel += read_bases[prev_index:index]
                index = ind + indel_size
                prev_index = index
            else:
                index += 1
        read_bases_no_indel += read_bases[prev_index:index]
        read_bases = read_bases_no_indel
    # counting matches and mismatches
    if ref == "C":
        incons_basecalls += read_bases.count('a') + \
                            read_bases.count('g') + \
                            read_bases.count('t') + \
                            read_bases.count('G') + \
                            read_bases.count('A')
        unconverted_c = read_bases.count(".")
        converted_c = read_bases.count("T")
        cons_basecalls = read_bases.count(',') + unconverted_c               
        return (str(cons_basecalls),str(incons_basecalls),unconverted_c,converted_c)
    elif ref == "G":
        incons_basecalls += read_bases.count('A') + \
                            read_bases.count('C') + \
                            read_bases.count('T') + \
                            read_bases.count('c') + \
                            read_bases.count('t')
        unconverted_c = read_bases.count(",")
        converted_c = read_bases.count("a")
        cons_basecalls = read_bases.count('.') + unconverted_c
        return (str(cons_basecalls),str(incons_basecalls),unconverted_c,converted_c)
    elif ref == "T":
        incons_basecalls += read_bases.count('a') + \
                            read_bases.count('c') + \
                            read_bases.count('g') + \
                            read_bases.count('A') + \
                            read_bases.count('C') + \
                            read_bases.count('G')
        cons_basecalls = read_bases.count(',')
        return (str(cons_basecalls),str(incons_basecalls))
    elif ref == "A":
        incons_basecalls += read_bases.count('T') + \
                            read_bases.count('C') + \
                            read_bases.count('G') + \
                            read_bases.count('t')+ \
                            read_bases.count('c')+ \
                            read_bases.count('g')
        cons_basecalls = read_bases.count('.')
        return (str(cons_basecalls),str(incons_basecalls))
    else:
        return ('0',str(incons_basecalls))

def call_methylated_sites_with_SNP_info(inputf, sample, reference_fasta,
                                        unmethylated_control=None,
                                        sig_cutoff=.01,num_procs = 1,
                                        num_upstr_bases=0,num_downstr_bases=2,
                                        generate_mpileup_file=True,
                                        compress_output=True,
                                        bgzip=False,
                                        path_to_bgzip="",
                                        path_to_tabix="",
                                        buffer_line_number=100000,
                                        min_mapq=30,
                                        min_cov=1,binom_test=True,
                                        path_to_samtools="",
                                        remove_chr_prefix=True,
                                        sort_mem="500M",
                                        add_snp_info=False,
                                        path_to_files="",min_base_quality=1,
                                        keep_temp_files=False):

    """
    inputf is the path to a bam file that contains mapped bisulfite sequencing reads
    
    sample is the name you'd like for the allc files. The files will be named like so:
        allc_<sample>_<chrom>.tsv
    
    reference is the path to a samtools indexed fasta file
    
    control is the name of the chromosome/region that you want to use to estimate the non-conversion rate of your 
        sample, or the non-conversion rate you'd like to use. Consequently, control is either a string, or a decimal
        If control is a string then it should be in the following format: "chrom:start-end". 
        If you'd like to specify an entire chromosome simply use "chrom:"
    
    sig_cutoff is a float indicating the adjusted p-value cutoff you wish to use for determining whether or not
        a site is methylated
    
    num_procs is an integer indicating how many num_procs you'd like to run this function over
    
    min_cov is an integer indicating the minimum number of reads for a site to be tested.
    
    path_to_files is a string indicating the path for the output and the input bam, mpileup, or allc files
        for methylation calling.
    min_base_quality is an integer indicating the minimum PHRED quality score for a base to be included in the
        mpileup file (and subsequently to be considered for methylation calling)
    """

    if binom_test and unmethylated_control is None:
        print_error("Please specify unmethylated_control if you would like to do binomial test!\n")

    #Figure out all the correct quality options based on the offset or CASAVA version given
    # quality_version >= 1.8:
    quality_base = 33
    if len(path_to_files)!=0:
        path_to_files+="/"
    if len(path_to_samtools)!=0:
        path_to_samtools+="/"

    try:
        num_procs = int(num_procs)
    except:
        sys.exit("num_procs must be an integer")
        
    try:
        #make sure bam file is indexed
        open(inputf+".bai",'r')
    except:
        print_checkpoint("Input not indexed. Indexing...")
        subprocess.check_call(shlex.split(path_to_samtools+"samtools index "+inputf))

    ## Check fasta index
    try:
        f = open(reference_fasta+".fai",'r')
    except:
        print("Reference fasta not indexed. Indexing.")
        try:
            subprocess.check_call(shlex.split(path_to_samtools+"samtools faidx "+reference_fasta))
            f = open(reference_fasta+".fai",'r')
        except:
            sys.exit("Reference fasta wasn't indexed, and couldn't be indexed. "
                     +"Please try indexing it manually and running methylpy again.")

    ## Input
    if not generate_mpileup_file:
        cmd = path_to_samtools+"samtools mpileup -Q "+str(min_base_quality)+\
              +" -q "+str(min_mapq)+" -B -f "+reference_fasta+" "+inputf
        pipes = subprocess.Popen(shlex.split(cmd),
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE,
                                 universal_newlines=True)
        fhandle = pipes.stdout
    else:
        with open(path_to_files+sample+"_mpileup_output.tsv",'w') as f:
            subprocess.check_call(
                shlex.split(
                    path_to_samtools+"samtools mpileup -Q "+str(min_base_quality)
                    +" -q "+str(min_mapq)
                    +" -B -f "+reference_fasta+" "+inputf),
                stdout=f)
        fhandle = open(path_to_files+sample+"_mpileup_output.tsv" ,'r')

    ## Output
    if compress_output:
        output_filehandler = gzip.open(path_to_files+"allc_"+sample+".tsv.gz",'wt')
        output_file = path_to_files+"allc_"+sample+".tsv.gz"
    else:
        output_filehandler = open(path_to_files+"allc_"+sample+".tsv",'w')
        output_file = path_to_files+"allc_"+sample+".tsv"
        
    complement = {"A":"T","C":"G","G":"C","T":"A","N":"N"}
    context_len = num_upstr_bases+1+num_downstr_bases
    cur_chrom = ""
    #cur_chrom_nochr = ""
    line_counts = 0
    out = ""
    SNP_info = {}
    SNP_info_end = 0
    to_chrom_end = False
    #for line in fhandle:
    line = True
    while line:
        line = fhandle.readline()
        fields = line.split("\t")
        # get reference genome information
        if fields[0] != cur_chrom:
            cur_chrom = fields[0]
            cur_chrom_nochr = cur_chrom
            if remove_chr_prefix and cur_chrom.startswith("chr"):
                cur_chrom_nochr = cur_chrom_nochr[3:]
            seq = get_chromosome_sequence(reference_fasta,cur_chrom)
            if seq != None:
                seq = seq.upper()
                to_chrom_end = False
                SNP_info_end = 0
        if seq == None:
            continue
        # get SNP information
        pos = int(fields[1])-1
        if pos > SNP_info_end-context_len-1 and not to_chrom_end:
            anchor = fhandle.tell()
            SNP_info_end = pos + 100000
            new_SNP_info = {}
            for tmp_pos in range(pos-context_len,pos+1):
                new_SNP_info[tmp_pos] = SNP_info.get(tmp_pos,('0','0'))
            SNP_info = new_SNP_info
            # make sure the current position is corrected analyzed
            # Otherwise, error happens if there is a region that has zero coverage
            SNP_info[pos] = analyze_read_basecalls(fields[2],fields[4])
            for line in fhandle:
                tmp_fields = line.split("\t")
                tmp_pos = int(tmp_fields[1])-1
                if tmp_pos > SNP_info_end:
                    break
                if tmp_fields[0] != cur_chrom:
                    to_chrom_end = True
                    break
                SNP_info[tmp_pos] = analyze_read_basecalls(tmp_fields[2],tmp_fields[4])
            fhandle.seek(anchor)
        # generate allc line
        if fields[2] == "C":
            strand = "+"
            try:
                context = seq[(pos-num_upstr_bases):(pos+num_downstr_bases+1)]
            except: # complete context is not available, skip
                continue
            incons_bases = ",".join(
                    [SNP_info.get(tmp_pos,('0','0'))[0]
                     for tmp_pos in range(pos-num_upstr_bases,pos+num_downstr_bases+1)]
                )
            incons_bases_cov = ",".join(
                    [SNP_info.get(tmp_pos,('0','0'))[1]
                     for tmp_pos in range(pos-num_upstr_bases,pos+num_downstr_bases+1)]
                )
            incons_base, incons_base_cov, unconverted_c, converted_c  = SNP_info[pos]
        elif fields[2] == "G":
            strand = "-"
            try:
                context = "".join([complement[base]
                                   for base in reversed(
                                           seq[(pos-num_downstr_bases):(pos+num_upstr_bases+1)]
                                   )]
                )
            except: # complete context is not available, skip
                continue
            incons_bases = ",".join([SNP_info.get(tmp_pos,('0','0'))[0]
                                     for tmp_pos in reversed(
                                             range(pos-num_downstr_bases,pos+num_upstr_bases+1)
                                     )]
            )
            incons_bases_cov = ",".join([SNP_info.get(tmp_pos,('0','0'))[1]
                                         for tmp_pos in reversed(
                                                 range(pos-num_downstr_bases,pos+num_upstr_bases+1)
                                         )]
            )
            incons_base, incons_base_cov, unconverted_c, converted_c  = SNP_info[pos]
        else:
            continue

        cov = unconverted_c+converted_c    
        if cov > 0 and len(context) == context_len:
            line_counts += 1
            out += "\t".join([cur_chrom_nochr,str(pos+1),strand,context,
                              str(unconverted_c),str(cov),"1",
                              incons_bases,incons_bases_cov])+"\n"
        if line_counts > buffer_line_number:
            output_filehandler.write(out)
            line_counts = 0
            out = ""

    if line_counts > 0:
        output_filehandler.write(out)
        line_counts = 0
        out = ""
    fhandle.close()
    output_filehandler.close()

    if generate_mpileup_file and not keep_temp_files:
        subprocess.check_call(shlex.split("rm -f "+path_to_files+sample+"_mpileup_output.tsv"))

    if binom_test:
        print_checkpoint('Perform binomial test')
        perform_binomial_test(allc_file=output_file,
                              sample=sample,
                              path_to_output=path_to_files,
                              unmethylated_control=unmethylated_control,
                              min_cov=min_cov,
                              sig_cutoff=sig_cutoff,
                              num_procs=num_procs,
                              sort_mem=sort_mem,
                              compress_output=compress_output,
                              buffer_line_number=buffer_line_number,
                              remove_chr_prefix=remove_chr_prefix)
    elif not unmethylated_control is None:
        non_conversion = calculate_non_conversion_rate(unmethylated_control,
                                                       output_file,
                                                       chrom_pointer=None,
                                                       remove_chr_prefix=remove_chr_prefix)
    index_allc_file(output_file)
    if bgzip:
        bgzip_allc_file(output_file,path_to_bgzip,path_to_tabix,buffer_line_number)
    return 0

def do_split_allc_file(allc_file,
                       sample,
                       path_to_output = "",
                       compress_output=True,
                       buffer_line_number=100000):
    """
    """
    if len(path_to_output)!=0:
        path_to_output+="/"

    fhandle = open_allc_file(allc_file)

    output_files = []
    cur_chrom = ""
    out = ""
    line_counts = 0
    for line in fhandle:
        fields = line.split("\t")
        if fields[0] != cur_chrom:
            # output
            if line_counts > 0:
                output_handle.write(out)
                line_counts = 0
                out = ""
                output_handle.close()
            cur_chrom = fields[0]
            if compress_output:
                output_handle = gzip.open(path_to_output+"allc_"+sample+"_"+cur_chrom+".tsv.gz",'wt')
                output_files.append(path_to_output+"allc_"+sample+"_"+cur_chrom+".tsv")
            else:
                output_handle = open(path_to_output+"allc_"+sample+"_"+cur_chrom+".tsv",'w')
                output_files.append(path_to_output+"allc_"+sample+"_"+cur_chrom+".tsv")
        # data
        line_counts += 1
        out += line
        if line_counts >= buffer_line_number:
            output_handle.write(out)
            line_counts = 0
            out = ""

    if line_counts > 0:
        output_handle.write(out)
        line_counts = 0
        out = ""
    output_handle.close()
    return(output_files)

def do_split_allc_file_chunk(allc_file,
                             sample,
                             num_chunks,
                             path_to_output = "",
                             compress_output=True):
    """
    """
    if len(path_to_output)!=0:
        path_to_output+="/"

    f = open_allc_file(allc_file)
    chrom_line_counts = {}
    total_line_count = 0
    for line in f:
        fields = line.split("\t")
        chrom_line_counts[fields[0]] = chrom_line_counts.get(fields[0],0)+1
        total_line_count += 1

    if total_line_count == 0:
        return 0

    chunk_size = math.ceil(float(total_line_count)/float(num_chunks))    

    index, cur_line_count = 0, 0
    chrom2index = {}
    for chrom in sorted(chrom_line_counts.keys(),key=lambda x:chrom_line_counts[x]):
        if cur_line_count >= chunk_size:
            cur_line_count = 0
            index += 1
        cur_line_count += chrom_line_counts[chrom]
        chrom2index[chrom] = index
        
    output_files = []
    output_handles = {}
    for index in range(index+1):
        if compress_output:
            output_handle = gzip.open(path_to_output+"allc_"+sample+"_"+str(index)+".tsv.gz",'wt')
            output_files.append(path_to_output+"allc_"+sample+"_"+str(index)+".tsv")
        else:
            output_handle = open(path_to_output+"allc_"+sample+"_"+str(index)+".tsv",'w')
            output_files.append(path_to_output+"allc_"+sample+"_"+str(index)+".tsv")
        output_handles[index] = output_handle

    f.seek(0)
    for line in f:
        fields = line.split("\t")
        output_handles[chrom2index[fields[0]]].write(line)

    for index in output_handles:
        output_handles[index].close()

    return(output_files)

def perform_binomial_test(allc_file,
                          sample,
                          path_to_output,
                          unmethylated_control,
                          min_cov=2,
                          sig_cutoff=0.01,
                          num_procs=1,
                          sort_mem="500M",
                          compress_output=True,
                          buffer_line_number=100000,
                          remove_chr_prefix=True):
    """
    """
    if len(path_to_output)!=0:
        path_to_output+="/"
    # calculate non-conversion rate
    chrom_pointer = read_allc_index(allc_file)
    non_conversion = calculate_non_conversion_rate(unmethylated_control,
                                                   allc_file,
                                                   chrom_pointer)    
    # binomial test
    if num_procs > 1:
        if len(chrom_pointer.keys()) > 100: # with too many files open can cause problem
            input_files = do_split_allc_file_chunk(allc_file,
                                                   sample,
                                                   min(100,num_procs),
                                                   path_to_output,
                                                   compress_output=False)
        else:
            # split allc file by chromosome
            input_files = do_split_allc_file(allc_file,
                                             sample,
                                             path_to_output,
                                             compress_output=False,
                                             buffer_line_number=buffer_line_number)
        output_files = [input_file+"_binom_results.tsv" for input_file in input_files]
        pool=multiprocessing.Pool(num_procs)
        results = []
        for input_file,output_file in zip(input_files,output_files):
            results.append(pool.apply_async(allc_run_binom_tests,
                                            (input_file,output_file,non_conversion),
                                            {"min_cov":min_cov,"sort_mem":sort_mem}))
        mc_class_counts = {}
        for result in results:
            result_mc_class_counts = result.get()
            for mc_class in result_mc_class_counts:
                mc_class_counts[mc_class] = mc_class_counts.get(mc_class,0) + result_mc_class_counts[mc_class]
        pool.close()
        pool.join()
        subprocess.check_call(shlex.split("rm "+" ".join(input_files)))
    else:
        output_files = [path_to_output+"allc_"+sample+".tsv_binom_results.tsv"]
        output_file = output_files[0]
        mc_class_counts = allc_run_binom_tests(filen=allc_file,
                                               output_file=output_file,
                                               non_conversion=non_conversion,
                                               min_cov=min_cov,
                                               sort_mem=sort_mem)

    # FDR correction
    p_value_cutoff = benjamini_hochberg_correction_call_methylated_sites(
        files=output_files,
        mc_class_counts=mc_class_counts,
        sig_cutoff=sig_cutoff)

    output_file = path_to_output+"allc_"+sample+".tsv"
    if compress_output:
        output_file += ".gz"

    filter_files_by_pvalue_combined(input_files=output_files,
                                    output_file=output_file,
                                    best_pvalues=p_value_cutoff,
                                    num_procs=num_procs,
                                    sort_mem=sort_mem,
                                    compress_output=compress_output)
    # remove _binom_results.tsv files
    subprocess.check_call(shlex.split("rm "+" ".join(output_files)))

def calculate_non_conversion_rate(unmethylated_control,
                                  allc_file,
                                  chrom_pointer=None,
                                  remove_chr_prefix=True):
    """
    """
    # Parse unmethylated_control
    try:
        non_conversion = float(unmethylated_control)
        if non_conversion < 0 or non_conversion > 1:
            print_error("Invalid unmethylated_control! "
                        +"It should be either a string, or a decimal between 0 and 1!\n")
        else:
            print_checkpoint("The non-conversion rate is "+str(non_conversion*100)+"%")
            return(non_conversion)
    except:
        if isinstance(unmethylated_control,str):
            fields = [field for field in
                      re.split("[\:\-]",unmethylated_control)
                      if len(field) > 0]
            um_chrom,um_start,um_end = None,None,None
            if len(fields) == 0:
                print_error("Invalid unmethylated_control! "
                            +"It should be either a string, or a decimal between 0 and 1!\n")
            # decode
            fields[0] = fields[0]
            if remove_chr_prefix and fields[0].startswith("chr"):
                fields[0] = fields[0][3:]
            if len(fields) == 1: # chrom only
                um_chrom = fields[0]
            elif len(fields) == 2: # chrom and start
                um_chrom,um_start = fields
            else:
                um_chrom,um_start,um_end = fields[:3]
            # further parsing
            try:
                if not (um_start is None):
                    um_start = int(um_start)
                    if not (um_end is None):
                        um_end = int(um_end)
            except:
                print_error("Invalid unmethylated_control! "
                            +"It should be either a string, or a decimal between 0 and 1!\n")

    # scan allc file to set up a table for fast look-up of lines belong
    # to different chromosomes
    chrom_pointer = read_allc_index(allc_file)
    f = open_allc_file(allc_file)
    if um_chrom not in chrom_pointer:
        print_error("The chromosome specified in unmethylated_control is not in the output allc file!\n")
    f.seek(chrom_pointer[um_chrom])
    
    # calculate non-conversion rate
    mc, h = 0, 0
    for line in f:
        line = line.rstrip("\n")
        fields = line.split("\t")
        if fields[0] != um_chrom:
            break
        if not (um_end is None) and int(fields[1]) > um_end:
            break
        if not (um_start is None) and int(fields[1]) < um_start:
            continue
        mc += int(fields[4])
        h += int(fields[5])

    non_conversion = None
    if h > 0:
        non_conversion = float(mc) / float(h)
        print_checkpoint("The non-conversion rate is "+str(non_conversion*100)+"%")
        return(non_conversion)
    else:
        print_error("The chromosome and range specified in unmethylated_control "
                    +"is not in the output allc file!\n")

def benjamini_hochberg_correction_call_methylated_sites(files,mc_class_counts,sig_cutoff):
    """
    This function is similar to the one defined here:
    http://stats.stackexchange.com/questions/870/multiple-hypothesis-testing-correction-with-benjamini-hochberg-p-values-or-q-va
    But takes advantage of the fact that the elements provided to it are in a sorted file.
    This way, it doesn't have to load much into memory.
    This link:
    http://brainder.org/2011/09/05/fdr-corrected-fdr-adjusted-p-values/
    was also helpful as the monotonicity correction from stats.stackexchange is not correct.
    
    file is a string indicating the path to an allc file (generated by run_binom_tests).
    
    mc_class_counts is a dictionary indicating the total number of statistical tests performed for each mc context
    
    sig_cutoff is the FDR cutoff you'd like to use to indicate if a site is significant.
    """
    #A dict of file_names to file handles for the benjamini hochberg correction step
    input_files = {}
    input_lines = {}
    input_fields ={}
    input_pvalues={}
    test_num={}
    prev_bh_value = {}
    best_fdr = {}
    best_pvalue = {}
    
    output_files = {}
    for filen in files:
        input_files[filen]=open(filen,'r')
        input_lines[filen] = input_files[filen].readline().rstrip()
        input_fields[filen] = input_lines[filen].split("\t")
        try:
            input_pvalues[filen] = float(input_fields[filen][6])
        except:
            #Dummy value that will never be the minimum
            input_pvalues[filen] = 2.0
    min_pvalue = min(input_pvalues,key=input_pvalues.get)
    #pdb.set_trace()
    while [i for i in input_pvalues if input_pvalues[i]!=2.0]:
        fields = input_fields[min_pvalue]
        bh_value = float(fields[6]) * mc_class_counts[fields[3]] / (test_num.get(fields[3],1) + 1)
        # Sometimes this correction can give values greater than 1,
        # so we set those values at 1
        bh_value = min(bh_value, 1.0)
        prev_bh_value[fields[3]] = bh_value
        #if bh_value <= sig_cutoff and bh_value >= best_fdr:
        if bh_value <= sig_cutoff:
            best_fdr[fields[3]] = bh_value
            best_pvalue[fields[3]] = float(fields[6])
        
        test_num[fields[3]] = test_num.get(fields[3],1) + 1
        input_lines[min_pvalue]=input_files[min_pvalue].readline().rstrip()
        input_fields[min_pvalue]=input_lines[min_pvalue].split("\t")
        try:
            input_pvalues[min_pvalue]=float(input_fields[min_pvalue][6])
        except:
            #Dummy value that will never be the minimum
            input_pvalues[min_pvalue]=2.0
        min_pvalue = min(input_pvalues,key=input_pvalues.get)
    for mc_class in mc_class_counts:
        best_pvalue[mc_class] = best_pvalue.get(mc_class,0)
        print("The closest p-value cutoff for "
              +mc_class+" at your desired FDR is "+
              str(best_pvalue[mc_class])+
              " which corresponds to an FDR of "+
              str(best_fdr.get(mc_class,1)))
    for filen in files:
        input_files[filen].close()
    return best_pvalue

def allc_run_binom_tests(filen,output_file,non_conversion,min_cov=1,sort_mem="500M"):
    """
    This function is used to recall methylated sites. This is faster than
    going through the original mpileup files.
    file is a string containing the path to an mpileup file
    
    non_conversion is a float indicating the estimated non-conversion rate and sequencing
        error
        
    min_cov is the minimum number of reads a site must have to be tested
    
    sort_mem is the parameter to pass to unix sort with -S/--buffer-size command
    """
    if sort_mem:
        sort_option = " -S " + sort_mem
    else:
        sort_option = ""

    mc_class_counts = {}
    obs_pvalues = {}
    reverse_complement = {"A":"T","C":"G","G":"C","T":"A","N":"N"}
    f = open_allc_file(filen)
    g = open(output_file,'w')
    for line in f:
        line = line.rstrip()
        fields = line.split("\t")
        mc_class = fields[3]
        unconverted_c = int(fields[4])
        converted_c = int(fields[5]) - unconverted_c
        total = int(fields[5])
        if total >= min_cov and unconverted_c != 0:
            try:
                p_value = obs_pvalues[(unconverted_c,total)]
            except:
                p_value = sci.binom.sf(unconverted_c-1,total,non_conversion)
                obs_pvalues[(unconverted_c,total)] = p_value
            g.write("\t".join(fields[:6])+"\t"+str(p_value)+"\n")
            mc_class_counts[mc_class] = mc_class_counts.get(mc_class,0) + 1
        elif total != 0:
            #a dummy value that will always sort to the bottom of the BH correction and be interpreted as
            #a unmethylated site
            p_value = 2.0
            g.write("\t".join(fields[:6])+"\t"+str(p_value)+"\n")

    f.close()
    g.close()
    subprocess.check_call(shlex.split("sort" + sort_option + " -k 7g,7g -o "+output_file+" "+output_file))
    return mc_class_counts

def filter_files_by_pvalue_combined(input_files,output_file,
                                    best_pvalues,num_procs,
                                    compress_output=True,
                                    sort_mem="500M"):
    """
    sort_mem is the parameter to pass to unix sort with -S/--buffer-size command
    """
    if sort_mem:
        sort_option = " -S " + sort_mem
    else:
        sort_option = ""
        
    print_checkpoint("Begin sorting file by position")
    # sort input files
    if num_procs > 1:
        pool=multiprocessing.Pool(num_procs)
        for input_file in input_files:
            pool.apply_async(subprocess.check_call,
                             (shlex.split(
                                 "sort" + sort_option + " -k 1,1 -k 2,2g -o "+input_file+" "+input_file),
                             ))
        pool.close()
        pool.join()
    else:
        for input_file in input_files:
            subprocess.check_call(
                shlex.split("sort" + sort_option + " -k 1,1 -k 2,2g -o "+input_file+" "+input_file))
    # output file
    if compress_output:
        if output_file[-3:] != ".gz":
            output_file += ".gz"
        g = gzip.open(output_file,'wt')
    else:
        g = open(output_file,'w')
    # write to output file
    for input_file in input_files:
        f = open(input_file,'r')

        for line in f:
            line = line.rstrip()
            fields = line.split("\t")
            if fields[6] != "2.0" and float(fields[6]) <= best_pvalues[fields[3]]:
                g.write("\t".join(fields[:6])+"\t1\n")
            else:
                g.write("\t".join(fields[:6])+"\t0\n")
        f.close()
    g.close()

def bam_quality_mch_filter(inputf,
                           outputf,
                           reference_fasta,
                           min_mapq=30,
                           min_ch=3,
                           max_mch_level=0.7,
                           buffer_line_number=100000,
                           path_to_samtools=""):
    """
    
    """

    ## Check fasta index
    try:
        f = open(reference_fasta+".fai",'r')
    except:
        print("Reference fasta not indexed. Indexing.")
        try:
            subprocess.check_call(shlex.split(path_to_samtools+"samtools faidx "+reference_fasta))
            f = open(reference_fasta+".fai",'r')
        except:
            sys.exit("Reference fasta wasn't indexed, and couldn't be indexed. "
                     +"Please try indexing it manually and running methylpy again.")

    min_ch = int(min_ch)
    max_mch_level = float(max_mch_level)
    
    # quality filter
    cmd = path_to_samtools+"samtools view -q"+str(min_mapq)+" "+inputf
    pipes = subprocess.Popen(shlex.split(cmd),
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE,
                             universal_newlines=True)
    fhandle = pipes.stdout

    # open up output file
    out_handle = open(outputf+".sam",'w')
    # print header
    subprocess.check_call(shlex.split(path_to_samtools+"samtools view -H "+inputf),
                          stdout=out_handle)

    line_counts = 0
    out = ""
    cur_chrom = ""
    for line in fhandle:
        fields = line.split("\t")
        if fields[2] != cur_chrom:
            cur_chrom = fields[2]
            seq = get_chromosome_sequence(reference_fasta,cur_chrom)
            if seq != None:
                seq = seq.upper()

        if seq == None:
            continue

        # get mc state
        read_pos = int(fields[3]) - 1
        uch = 0
        mch = 0
        num_ch = 0
        if int(fields[1]) & 16 == 0: # + strand
            for ind,base in enumerate(fields[9]):
                pos = ind + read_pos
                try:
                    if seq[pos] != "C":
                        continue
                    if base == "C":
                        if seq[pos+1] == "G":
                            continue
                        else:
                            mch += 1
                    elif base == "T":
                        if seq[pos+1] == "G":
                            continue
                        else:
                            uch += 1
                except: # pos + 1 exceed reference boundary
                    pass
        else: # - strand
            for ind,base in enumerate(fields[9]):
                pos = ind + read_pos
                try:
                    if seq[pos] != "G":
                        continue
                    if base == "G":
                        if seq[pos-1] == "C":
                            continue
                        else:
                            mch += 1
                    elif base == "A":
                        if seq[pos-1] == "C":
                            continue
                        else:
                            uch += 1
                except: # pos - 1 exceed reference boundary
                    pass
        # apply filter
        tot_ch = float(mch+uch)
        if tot_ch >= min_ch and float(mch)/float(tot_ch) > max_mch_level:
            continue
        
        line_counts += 1
        out += line
        if line_counts > buffer_line_number:
            out_handle.write(out)
            out = ""
            line_counts = 0
        
    if line_counts > 0:
        out_handle.write(out)
        out = ""
    out_handle.close()
    
    # sam to bam
    f = open(outputf,'w')
    subprocess.check_call(shlex.split(path_to_samtools+"samtools view -Sb "+outputf+".sam"),stdout=f)
    f.close()
    
    # remove sam
    subprocess.check_call(shlex.split("rm " +outputf+".sam"))
