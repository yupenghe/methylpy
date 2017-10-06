import sys
import multiprocessing
import subprocess
import scipy.stats as sci
from scipy.stats.mstats import mquantiles
from methylpy.utilities import print_checkpoint, print_error, split_fastq_file
import pdb
import shlex
import itertools
import re
import glob
import cStringIO as cStr
import bisect
from methylpy.call_mc_se import call_methylated_sites, remove_clonal_bam
try:
    from argparse import ArgumentParser
except Exception, e:
    exc_type, exc_obj, exc_tb = sys.exc_info()
    print(exc_type, exc_tb.tb_lineno)
    print(e)
    exit("methylpy.call_mc_pe requires ArgumentParser from the argparse module")
# bz2
try:
    import bz2
except Exception, e:
    exc_type, exc_obj, exc_tb = sys.exc_info()
    print(exc_type, exc_tb.tb_lineno)
    print e
    sys.exit("methylpy.call_mc_pe requires the bz2 module")
# gzip
try:
    import gzip
except Exception, e:
    exc_type, exc_obj, exc_tb = sys.exc_info()
    print(exc_type, exc_tb.tb_lineno)
    print e
    sys.exit("methylpy.call_mc_pe requires the gzip module")

def run_methylation_pipeline_pe(read1_files, read2_files, libraries, sample,
                                forward_reference, reverse_reference, reference_fasta,
                                unmethylated_control="chrL:",
                                path_to_output="", sig_cutoff=0.01,
                                num_procs=1, sort_mem="500M",
                                num_upstr_bases=0,
                                num_downstr_bases=2,
                                generate_mpileup_file=True,
                                compress_output=True,
                                split_allc_files=False,
                                binom_test=True, bh=True, min_cov=2,
                                trim_reads=True, path_to_cutadapt="",
                                bowtie2=True, path_to_aligner="", aligner_options=[],
                                pbat=False,
                                remove_clonal=True, keep_clonal_stats=False,
                                path_to_picard="",java_options="-Xmx20g",
                                path_to_samtools="",
                                adapter_seq_read1="AGATCGGAAGAGCACACGTCTGAAC",
                                adapter_seq_read2="AGATCGGAAGAGCGTCGTGTAGGGA",
                                max_adapter_removal=None,
                                overlap_length=None, zero_cap=None,
                                error_rate=None, min_qual_score=10,
                                min_read_len=30,
                                keep_temp_files=False,
                                min_base_quality=1):

    """
    This function

    read1_files and read2_files are lists of fastq files of the forward and reverse reads
        from paired-end bisulfite sequencing data, which you'd like to run through the pipeline.
        The length of read1_files and read2_files should be the same. Also, Files in there two
        lists should be ordered such that the forward reads for a particular read set are in the
        same position of read1_files as the reverse reads are in the read2_files.
        (i.e. the elements in each of these lists are paired)
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
        Default: "chrL:"

    path_to_output is the path to a directory where you would like the output to be stored.
        The default is the same directory as the input fastqs.

    num_procs is an integer indicating how many num_procs you'd like to run this function over

    sort_mem is the parameter to pass to unix sort with -S/--buffer-size command. Default: "500M"

    sig_cutoff is a float indicating the adjusted p-value cutoff you wish to use for
        determining whether or not a site is methylated.

    pbat indicates whether the input data is from PBAT protocol instead of MethylC-seq protocol.
        Default: False. To be implemented.

    binom_tests indicates that you'd like to use a binomial test, rather than the alternative
        method outlined here: https://bitbucket.org/schultzmattd/methylpy/wiki/Methylation%20Calling

    min_cov is an integer indicating the minimum number of reads for a site to be tested.

    trim_reads is a boolean indicating that you want to have reads trimmed by cutadapt.

    path_to_cutadapt is the path to the cutadapt execuatable. Otherwise this is assumed to be in
        your path.

    bowtie2 specifies whether to use the bowtie2 aligner instead of bowtie. Default: True

    path_to_aligner is a string indicating the path to the folder in which aligner resides. Aligner
        is assumed to be in your path if this option isn't used

    aligner_options is a list of strings indicating options you'd like passed to aligner.
        (default for bowtie2: "-X 1000 -k 2 --no-mixed --no-discordant")
        (default for bowtie: "-X 1000 -S -k 1 -m 1 --best --strata --chunkmbs 64 -n 1")

    remove_clonal is a boolean indicating that you want to remove clonal reads (PCR duplicates).
        If true, executable picard should be available in folder specified in path_to_picard.

    path_to_picard is a string of the path to "picard.jar". "picard.jar" is assumed to be
        in your path if this option isn't used.

    path_to_samtools is a string indicating the path to the directory containing your
        installation of samtools. Samtools is assumed to be in your path if this is not
        provided.

    adapter_seq_read1:
        Sequence of an adapter that was ligated to the 3' end of read 1. The adapter itself
        and anything that follows is trimmed.

    adapter_seq_read2:
        Sequence of an adapter that was ligated to the 3' end of read 2. The adapter itself
        and anything that follows is trimmed.

    max_adapter_removal indicates the maximum number of times to try to remove adapters. Useful
        when an adapter gets appended multiple times.

    overlap_length is the minimum overlap length. If the overlap between the read and the
        adapter is shorter than LENGTH, the read is not modified. This reduces the no. of bases
        trimmed purely due to short random adapter matches.

    zero_cap causes negative quality values to be set to zero (workaround to avoid segmentation
        faults in BWA).

    error_rate is the maximum allowed error rate (no. of errors divided by the length of the
        matching region). Default: 0.1

    min_qual_score allows you to trim low-quality ends from reads before adapter removal.
        The algorithm is the same as the one used by BWA (Subtract CUTOFF from all qualities;
        compute partial sums from all indices to the end of the sequence; cut sequence at the
        index at which the sum is minimal).

    min_read_len indicates the minimum length a read must be to be kept. Reads that are too short even
        before adapter removal are also discarded. In colorspace, an initial primer is not counted.
        It is not recommended to change this value in paired-end processing because it may results in
        the situation that one of the two reads in a pair is discarded. And this will lead to error.

    keep_temp_files is a boolean indicating that you'd like to keep the intermediate files generated
        by this function. This can be useful for debugging, but in general should be left False.

    min_base_quality is an integer indicating the minimum PHRED quality score for a base to be included
        in the mpileup file (and subsequently to be considered for methylation calling)    
    """

    #Default bowtie option
    if len(aligner_options) == 0:
        if bowtie2:
            aligner_options = ["-X 1000", "-k 2", "--no-discordant", "--no-mixed"]
        else:
            aligner_options = ["-X 1000", "-S", "-k 1", "-m 1", "--best", "--strata",
                               "--chunkmbs 3072", "-n 1", "-e 100"]

    # CASAVA >= 1.8
    aligner_options.append("--phred33-quals")
    quality_base = 33

    #just to avoid any paths missing a slash. It's ok to add an extra slash if
    #there's already one at the end
    if len(path_to_samtools) != 0:
        path_to_samtools += "/"
    if len(path_to_aligner) != 0:
        path_to_aligner += "/"
    if len(path_to_output) != 0:
        path_to_output += "/"

    if sort_mem:
        if sort_mem.find("-S") == -1:
            sort_mem = " -S " + sort_mem
    else:
        sort_mem = ""

    # This code allows the user to supply paths with "*" in them rather than listing
    # out every single file
    total_input = 0
    total_unique = 0
    total_clonal = 0

    # Get expanded file list
    if pbat:
        expanded_read1_file_list, expanded_library_list = expand_input_files(read2_files,libraries)
        expanded_read2_file_list, expanded_library_list = expand_input_files(read1_files,libraries)
    else:
        expanded_read1_file_list, expanded_library_list = expand_input_files(read1_files,libraries)
        expanded_read2_file_list, expanded_library_list = expand_input_files(read2_files,libraries)

    #Processing
    for current_library in set(libraries):
        library_read1_files = [filen for filen, library
                               in zip(expanded_read1_file_list, expanded_library_list)
                               if library == current_library]
        library_read2_files = [filen for filen, library
                               in zip(expanded_read2_file_list, expanded_library_list)
                               if library == current_library]

        #deal with actual filename rather than path to file
        lib_input, lib_unique = run_mapping_pe(
            current_library, library_read1_files, library_read2_files, sample,
            forward_reference, reverse_reference, reference_fasta,
            path_to_output=path_to_output,
            path_to_samtools=path_to_samtools, path_to_aligner=path_to_aligner,
            aligner_options=aligner_options, num_procs=num_procs,
            trim_reads=trim_reads, path_to_cutadapt=path_to_cutadapt,
            adapter_seq_read1=adapter_seq_read1, adapter_seq_read2=adapter_seq_read2,
            max_adapter_removal=max_adapter_removal, overlap_length=overlap_length,
            zero_cap=zero_cap, quality_base=quality_base, error_rate=error_rate,
            min_qual_score=min_qual_score, min_read_len=min_read_len,
            keep_temp_files=keep_temp_files,
            bowtie2=bowtie2,
            sort_mem=sort_mem)

        total_input += lib_input
        total_unique += lib_unique

        ## Remove clonal reads
        if remove_clonal == True:
            lib_clonal = remove_clonal_bam(input_bam = path_to_output+sample+"_"+str(current_library)
                                           +"_processed_reads.bam",

                                           output_bam = path_to_output+sample+"_"+str(current_library)
                                           +"_processed_reads_no_clonal.bam",                                           

                                           metric = path_to_output+sample+"_"+
                                           str(current_library)+".metric",

                                           is_pe = True,
                                           path_to_picard=path_to_picard,
                                           java_options=java_options)

            subprocess.check_call(shlex.split("rm "+path_to_output+sample+"_"+
                                              str(current_library)+"_processed_reads.bam"))
            if not keep_clonal_stats:
                subprocess.check_call(shlex.split("rm "+" "+path_to_output+sample+"_"+
                                                  str(current_library)+".metric"))
            total_clonal += lib_clonal

    print_checkpoint("There are " + str(total_input) + " total input reads")
    print_checkpoint("There are " + str(total_unique) + " uniquely mapping read pairs, " +
                     str(float(total_unique) / total_input*100) + " percent remaining")

    if remove_clonal == True:
        total_non_clonal = total_unique - total_clonal
        print_checkpoint("There are "+str(total_non_clonal)+" non-clonal reads, "+
                         str(float(total_non_clonal) / total_input*100)+" percent remaining")
        ## Merge bam files to get final bam file
        library_files = [path_to_output+sample+"_"+str(library)+"_processed_reads_no_clonal.bam"
                         for library in set(libraries)]
        if len(library_files) > 1:
            merge_bam_files(library_files, path_to_output+sample+
                            "_processed_reads_no_clonal.bam", path_to_samtools)
            subprocess.check_call(shlex.split("rm "+" ".join(library_files)))
        else:
            subprocess.check_call(
                shlex.split("mv "+library_files[0]+" "+
                            path_to_output+sample+"_processed_reads_no_clonal.bam")
            )
    ## If not removing clonal reads
    else:
        library_files = [path_to_output+sample+"_"+str(library)+"_processed_reads.bam"
                         for library in set(libraries)]
        if len(library_files) > 1:
            merge_bam_files(library_files,path_to_output+sample+"_processed_reads.bam",path_to_samtools)
            subprocess.check_call(shlex.split("rm "+" ".join(library_files)))
        else:
            subprocess.check_call(shlex.split("mv "+library_files[0]+" "+path_to_output+sample+"_processed_reads.bam"))

    #Calling methylated sites
    print_checkpoint("Begin calling mCs")
    if remove_clonal == True:
        output_bam_file = path_to_output+sample+"_processed_reads_no_clonal.bam"
    else:
        output_bam_file = path_to_output+sample+"_processed_reads.bam"

    call_methylated_sites_pe(output_bam_file,
                             sample,
                             reference_fasta,
                             unmethylated_control,
                             sig_cutoff=sig_cutoff,
                             num_procs=num_procs,
                             num_upstr_bases=num_upstr_bases,
                             num_downstr_bases=num_downstr_bases,
                             generate_mpileup_file=generate_mpileup_file,
                             compress_output=compress_output,
                             split_allc_files=split_allc_files,
                             min_cov=min_cov,
                             binom_test=binom_test,
                             bh=bh,
                             sort_mem=sort_mem,
                             path_to_files=path_to_output,
                             path_to_samtools=path_to_samtools,
                             min_base_quality=min_base_quality)
    print_checkpoint("Done")

def run_mapping_pe(current_library, library_read1_files, library_read2_files,
                   sample, forward_reference, reverse_reference, reference_fasta,
                   path_to_output="",
                   path_to_samtools="", path_to_aligner="",
                   aligner_options=[],
                   num_procs=1, trim_reads=True, path_to_cutadapt="",
                   adapter_seq_read1="AGATCGGAAGAGCACACGTCTGAAC",
                   adapter_seq_read2="AGATCGGAAGAGCGTCGTGTAGGGA",
                   max_adapter_removal=None, overlap_length=None, zero_cap=None,
                   quality_base=None, error_rate=None, min_qual_score=10,
                   min_read_len=30, keep_temp_files=False,
                   bowtie2=True, sort_mem="500M"):
    """
    This function runs the mapping portion of the methylation calling pipeline.
    For Paired-end data processing.

    current_library is the ID that you'd like to run mapping on.

    library_read1_files is a list of library IDs (in the same order as the files list) indiciating
    which libraries each set of fastq files belong to. If you use a glob, you only need to 
    indicate the library ID for those fastqs once (i.e., the length of files and libraries 
    should be the same)

    library_read2_files is a list of library IDs (in the same order as the files list) indiciating
    which libraries each set of fastq files belong to. If you use a glob, you only need to 
    indicate the library ID for those fastqs once (i.e., the length of files and libraries 
    should be the same)
    
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
        is assumed to be in your path if this option isn't used
            
    aligner_options is a list of strings indicating options you'd like passed to bowtie2 (or bowtie)
    
    num_procs is an integer indicating how many num_procs you'd like to run this function over
    
    trim_reads is a boolean indicating that you want to have reads trimmed by cutadapt

    path_to_cutadapt is the path to the cutadapt execuatable. Otherwise this is assumed to be in your
        path.

    adapter_seq_read1:
        Sequence of an adapter that was ligated to the 3' end of read 1. The adapter itself and anything that follows is
        trimmed.

    adapter_seq_read2:
        Sequence of an adapter that was ligated to the 3' end of read 2. The adapter itself and anything that follows is
        trimmed.

    max_adapter_removal indicates the maximum number of times to try to remove adapters. Useful when an adapter 
        gets appended multiple times.
        
    overlap_length is the minimum overlap length. If the overlap between the read and the adapter is shorter than 
        LENGTH, the read is not modified. This reduces the no. of bases trimmed purely due to short random adapter matches.

    zero_cap causes negative quality values to be set to zero (workaround to avoid segmentation faults in BWA).

    quality_base is the offset for quality scores. In other words, assume that quality values are encoded as ascii(quality + QUALITY_BASE). 
        The default (33) is usually correct, except for reads produced by some versions of the Illumina pipeline, where this should
        be set to 64.
        
    error_rate is the maximum allowed error rate (no. of errors divided by the length of the matching region) 
        (default: 0.1)
    
    min_qual_score allows you to trim low-quality ends from reads before adapter removal. The algorithm is the same as the 
        one used by BWA (Subtract CUTOFF from all qualities; compute partial sums from all indices to the end of the
        sequence; cut sequence at the index at which the sum is minimal).
        
    min_read_len indicates the minimum length a read must be to be kept. Reads that are too short even before adapter removal
        are also discarded. In colorspace, an initial primer is not counted.

    keep_temp_files is a boolean indicating that you'd like to keep the intermediate files generated
        by this function. This can be useful for debugging, but in general should be left False.
        
    bowtie2 specifies whether to use the bowtie2 aligner instead of bowtie
    
    sort_mem is the parameter to pass to unix sort with -S/--buffer-size command
    """
    if len(path_to_output) !=0:
        path_to_output+="/"
        
    total_input = 0
    total_unique = 0
    prefix = path_to_output+sample+"_"+str(current_library)

    #Split files
    print_checkpoint("Begin splitting reads for "+str(current_library))
    total_input_read1 = split_fastq_file(num_procs,library_read1_files,prefix+"_read1_split_")
    total_input_read2 = split_fastq_file(num_procs,library_read2_files,prefix+"_read2_split_")

    #Check if there are same number of reads in read 1 and read 2
    if total_input_read1 != total_input_read2:
        print_error("There are different numbers of read 1 and read 2 "+
                    "for library \"" + str(current_library) + "\" !")
    else:
        total_input = total_input_read1

    if trim_reads:
        #Trimming
        print_checkpoint("Begin trimming reads for "+str(current_library))  
        quality_trim_pe(
            inputf_read1=[prefix+"_read1_split_"+str(i) for i in xrange(0,num_procs)],
            outputf_read1=[prefix+"_read1_split_trimmed_"+str(i) for i in xrange(0,num_procs)],
            inputf_read2=[prefix+"_read2_split_"+str(i) for i in xrange(0,num_procs)],
            outputf_read2=[prefix+"_read2_split_trimmed_"+str(i) for i in xrange(0,num_procs)],
            adapter_seq_read1=adapter_seq_read1,
            adapter_seq_read2=adapter_seq_read2,
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
        
        subprocess.check_call(shlex.split("rm "+" ".join([prefix+"_read1_split_"+str(i) for i in xrange(0,num_procs)])))
        subprocess.check_call(shlex.split("rm "+" ".join([prefix+"_read2_split_"+str(i) for i in xrange(0,num_procs)])))
        
        #Conversion
        print_checkpoint("Begin converting reads for "+str(current_library))
        if num_procs > 1:
            pool = multiprocessing.Pool(num_procs)#read1
            for inputf,output in zip([prefix+"_read1_split_trimmed_"+str(i) for i in xrange(0,num_procs)],
                                     [prefix+"_read1_split_trimmed_converted_"+str(i) for i in xrange(0,num_procs)]):
                pool.apply_async(convert_reads_pe,(inputf,output))
            for inputf,output in zip([prefix+"_read2_split_trimmed_"+str(i) for i in xrange(0,num_procs)],
                                     [prefix+"_read2_split_trimmed_converted_"+str(i) for i in xrange(0,num_procs)]):
                pool.apply_async(convert_reads_pe,(inputf,output,True))
            pool.close()
            pool.join()
        else:
            for inputf,output in zip([prefix+"_read1_split_trimmed_"+str(i) for i in xrange(0,num_procs)],
                                     [prefix+"_read1_split_trimmed_converted_"+str(i) for i in xrange(0,num_procs)]):
                convert_reads_pe(inputf,output)
            for inputf,output in zip([prefix+"_read2_split_trimmed_"+str(i) for i in xrange(0,num_procs)],
                                     [prefix+"_read2_split_trimmed_converted_"+str(i) for i in xrange(0,num_procs)]):
                convert_reads_pe(inputf,output,True)
                
        subprocess.check_call(
            shlex.split("rm "+" ".join([prefix+"_read1_split_trimmed_"+str(i) for i in xrange(0,num_procs)]))
        )
        subprocess.check_call(
            shlex.split("rm "+" ".join([prefix+"_read2_split_trimmed_"+str(i) for i in xrange(0,num_procs)]))
        )        
        #Run bowtie
        input_fastq_read1 = [prefix+"_read1_split_trimmed_converted_"+str(i) for i in xrange(0,num_procs)]
        input_fastq_read2 = [prefix+"_read2_split_trimmed_converted_"+str(i) for i in xrange(0,num_procs)]
    else:
        print_checkpoint("No trimming applied on reads")  
        #Conversion
        print_checkpoint("Begin converting reads for "+str(current_library))
        if num_procs > 1:
            pool = multiprocessing.Pool(num_procs)#read1
            for inputf,output in zip([prefix+"_read1_split_"+str(i) for i in xrange(0,num_procs)],
                                     [prefix+"_read1_split_converted_"+str(i) for i in xrange(0,num_procs)]):
                pool.apply_async(convert_reads_pe,(inputf,output))
            for inputf,output in zip([prefix+"_read2_split_"+str(i) for i in xrange(0,num_procs)],
                                     [prefix+"_read2_split_converted_"+str(i) for i in xrange(0,num_procs)]):
                pool.apply_async(convert_reads_pe,(inputf,output,True))
            pool.close()
            pool.join()
        else:
            for inputf,output in zip([prefix+"_read1_split_"+str(i) for i in xrange(0,num_procs)],
                                     [prefix+"_read1_split_converted_"+str(i) for i in xrange(0,num_procs)]):
                convert_reads_pe(inputf,output)
            for inputf,output in zip([prefix+"_read2_split_"+str(i) for i in xrange(0,num_procs)],
                                     [prefix+"_read2_split_converted_"+str(i) for i in xrange(0,num_procs)]):
                convert_reads_pe(inputf,output,True)
        subprocess.check_call(shlex.split("rm "+" ".join([prefix+"_read1_split_"+str(i) for i in xrange(0,num_procs)])))
        subprocess.check_call(shlex.split("rm "+" ".join([prefix+"_read2_split_"+str(i) for i in xrange(0,num_procs)])))
        input_fastq_read1 = [prefix+"_read1_split_converted_"+str(i) for i in xrange(0,num_procs)]
        input_fastq_read2 = [prefix+"_read2_split_converted_"+str(i) for i in xrange(0,num_procs)]
                    
    #Run bowtie
    if bowtie2:
        print_checkpoint("Begin Running Bowtie2 for "+current_library)
    else:
        print_checkpoint("Begin Running Bowtie for "+current_library)
    total_unique = run_bowtie_pe(current_library,
                                  input_fastq_read1,
                                  input_fastq_read2,
                                  sample,
                                  forward_reference,reverse_reference,reference_fasta,
                                  path_to_output=path_to_output,
                                  path_to_samtools=path_to_samtools,
                                  aligner_options=aligner_options,
                                  path_to_aligner=path_to_aligner,num_procs=num_procs,
                                  keep_temp_files=keep_temp_files,
                                  bowtie2=bowtie2, sort_mem=sort_mem)
    
    return total_input,total_unique

def run_bowtie_pe(current_library,library_read1_files,library_read2_files,
                  sample,
                  forward_reference,reverse_reference,reference_fasta,
                  path_to_output="",                  
                  path_to_samtools="",
                  aligner_options="",path_to_aligner="",
                  num_procs=1,keep_temp_files=False, bowtie2=True, sort_mem="500M"):
    """
    This function runs bowtie on the forward and reverse converted bisulfite references 
    (generated by build_ref). The function is for processing paired-end data. It removes 
    any read that maps to both the forward and reverse strands. In addition, any inproperly 
    paired reads are removed. 
    
    library_read1_files is a list of fastq file paths of the first reads to be mapped

    library_read2_files is a list of fastq file paths of the second reads to be mapped
    
    forward_reference is a string indicating the path to the forward strand reference created by
        build_ref
    
    reverse_reference is a string indicating the path to the reverse strand reference created by
        build_ref
    
    prefix is a string that you would like prepended to the output files (e.g., the sample name)
    
    aligner_options is a list of strings indicating options you'd like passed to the aligner
        (bowtie2 or bowtie)
    
    path_to_aligner is a string indicating the path to the folder in which bowtie resides. Bowtie
        is assumed to be in your path if this option isn't used
    
    num_procs is an integer indicating the number of processors you'd like used for removing multi
        mapping reads and for bowtie mapping
    
    keep_temp_files is a boolean indicating that you'd like to keep the intermediate files generated
        by this function. This can be useful for debugging, but in general should be left False.
    
    bowtie2 specifies whether to use the bowtie2 aligner instead of bowtie
    
    sort_mem is the parameter to pass to unix sort with -S/--buffer-size command
    """
    options = aligner_options
    if " ".join(options).find(" -p ") == -1:
        options.append("-p "+str(num_procs))

    if len(path_to_output) !=0:
        path_to_output+="/"
        
    prefix = path_to_output+sample+"_"+str(current_library)
        
    ## Forward
    if bowtie2:
        args = [path_to_aligner+"bowtie2"]
        args.extend(options)
        args.append("--norc")
        args.append("-x "+forward_reference)
        args.append("-1 "+",".join(library_read1_files))
        args.append("-2 "+",".join(library_read2_files))
        args.append("-S "+prefix+"_forward_strand_hits.sam")
    else:
        args = [path_to_aligner+"bowtie"]
        args.extend(options)
        args.append("--norc")
        args.append(forward_reference)
        args.append("-1 "+",".join(library_read1_files))
        args.append("-2 "+",".join(library_read2_files))
        args.append(prefix+"_forward_strand_hits.sam")
    subprocess.check_call(shlex.split(" ".join(args)))
    print_checkpoint("Processing forward strand hits")
    find_multi_mappers_pe(prefix+"_forward_strand_hits.sam",prefix,
                          num_procs=num_procs,keep_temp_files=keep_temp_files)

    ## Reverse    
    if bowtie2:
        args = [path_to_aligner+"bowtie2"]
        args.extend(options)
        args.append("--nofw")
        args.append("-x "+reverse_reference)
        args.append("-1 "+",".join(library_read1_files))
        args.append("-2 "+",".join(library_read2_files))
        args.append("-S "+prefix+"_reverse_strand_hits.sam")
    else:
        args = [path_to_aligner+"bowtie"]
        args.extend(options)
        args.append("--nofw")
        args.append(reverse_reference)
        args.append("-1 "+",".join(library_read1_files))
        args.append("-2 "+",".join(library_read2_files))
        args.append(prefix+"_reverse_strand_hits.sam")
    subprocess.check_call(shlex.split(" ".join(args)))
    print_checkpoint("Processing reverse strand hits")
    find_multi_mappers_pe(prefix+"_reverse_strand_hits.sam",prefix,
                          num_procs=num_procs,append=True,keep_temp_files=keep_temp_files)
    
    ## Clear temporary files
    if keep_temp_files==False:
        subprocess.check_call(shlex.split("rm "+" ".join(library_read1_files+library_read2_files)))
    if num_procs > 1:
        pool = multiprocessing.Pool(num_procs)
        for file_num in xrange(0,num_procs):
            pool.apply_async(subprocess.check_call,(
                shlex.split("env LC_COLLATE=C sort" + sort_mem +
                            " -t '\t' -k 1 -o "+prefix+"_sorted_"+str(file_num)+" "+
                            prefix+"_sorted_"+str(file_num)),))
        pool.close()
        pool.join()
    else:
        for file_num in xrange(0,num_procs):
            subprocess.check_call(shlex.split("env LC_COLLATE=C sort" + sort_mem +
                                              " -t '\t' -k 1 -o "+
                                              prefix+"_sorted_"+str(file_num)+" "+
                                              prefix+"_sorted_"+str(file_num)))
    print_checkpoint("Finding multimappers")

    total_unique = merge_sorted_multimap_pe(current_library,
                                            [prefix+"_sorted_"+str(file_num) for file_num in xrange(0,num_procs)],
                                            prefix,
                                            reference_fasta,
                                            path_to_samtools="")
    subprocess.check_call(shlex.split("rm "+" ".join([prefix+"_sorted_"+str(file_num) for file_num in xrange(0,num_procs)])))
    return total_unique


def find_multi_mappers_pe(inputf,output,num_procs=1,keep_temp_files=False,append=False):
    """
    This function takes a sam file generated by bowtie/bowtie2 and pulls out any mapped reads.
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
    sam_header = []
    file_handles = {}
    f = open(inputf,'r')
    cycle = itertools.cycle(range(0,num_procs))
    for file_num in xrange(0,num_procs):
        if append == False:
            file_handles[file_num]=open(output+"_sorted_"+str(file_num),'w')
        else:
            file_handles[file_num]=open(output+"_sorted_"+str(file_num),'a')
    for line in f:
        if line[0] == "@":
            continue
        
        fields = line.split("\t")

        ## Check if it is proper pair
        if int(fields[1]) & 2 == 0:
            continue;
        
        header = fields[0].split("!")
        if (int(fields[1]) & 16) == 16:
            strand = "-"
        else:
            strand = "+"
        if (int(fields[1]) & 128) == 128:
            is_read2 = True
        else:
            is_read2 = False
        seq = decode_converted_positions(fields[9],header[-1],strand,is_read2)
        file_handles[cycle.next()].write(" ".join(header[:-1])+"\t"+"\t".join(fields[1:9])+"\t"+seq+"\t"+"\t".join(fields[10:]))
            #file_handles[cycle.next()].write("\t".join(fields[0:9])+"\t"+seq+"\t"+"\t".join(fields[10:]))
    f.close()
    if keep_temp_files == False:
        subprocess.check_call(shlex.split("rm "+inputf))
        pass
    for file_num in xrange(0,num_procs):
        file_handles[file_num].close()

def merge_sorted_multimap_pe(current_library,files,prefix,reference_fasta,path_to_samtools=""):
    """
    This function takes the files from find_multi_mappers and outputs the uniquely mapping reads
    
    files is a list of filenames containing the output of find_multi_mappers
    
    output is a prefix you'd like prepended to the file containing the uniquely mapping reads
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
        print "Reference fasta not indexed. Indexing."
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
        fields[filen] = lines[filen].split("\t")[0]#Read ID
    while True:
        all_fields = [field for field in fields.values() if field != ""]
        if len(all_fields) == 0:
            break
        min_field = min(all_fields)
        count_1 = 0
        count_2 = 0
        current_line_1 = ""
        current_line_2 = ""
        for key in fields:
            while fields[key] == min_field:
                if(int(lines[key].split("\t")[1]) & 64 == 64): #First in pair
                    count_1 += 1
                    current_line_1 = lines[key]
                else:
                    count_2 += 1
                    current_line_2 = lines[key]
                lines[key]=file_handles[key].readline()
                fields[key]=lines[key].split("\t")[0]
        #Check if there is only one valid alignment
        if count_1 == 1:
            output_handle.write(current_line_1)
            output_handle.write(current_line_2)
            total_unique += 1

    #output_pipe.stdin.close()
    output_handle.close()
    
    for index,filen in enumerate(files):
        file_handles[filen].close()

    f = open(output_bam_file,'w')
    subprocess.check_call(shlex.split(path_to_samtools+"samtools view -S -b -h "+output_sam_file),stdout=f)
    f.close()

    subprocess.check_call(shlex.split("rm "+output_sam_file))
    subprocess.check_call(shlex.split(path_to_samtools+"samtools sort "+output_bam_file+
                                      " -o "+output_bam_file))
    
    return total_unique

def convert_reads_pe(inputf,output,is_read2=False):
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
    encoding = encode_converted_positions(seq,is_read2=is_read2)
    
    if is_read2 == False:
        while header:
            g.write(header+"!"+encoding+"\n")
            converted_seq = seq.replace("C","T")
            g.write(converted_seq)
            g.write(header2)
            g.write(qual)            
            header = f.readline().rstrip()
            header = header.replace(" ","!")
            seq = f.readline()
            header2 = f.readline()
            qual = f.readline()
            encoding = encode_converted_positions(seq,is_read2=is_read2)
    else:
        while header:
            g.write(header+"!"+encoding+"\n")
            converted_seq = seq.replace("G","A")
            g.write(converted_seq)
            g.write(header2)
            g.write(qual)                    
            header = f.readline().rstrip()
            header = header.replace(" ","!")
            seq = f.readline()
            header2 = f.readline()
            qual = f.readline()
            encoding = encode_converted_positions(seq,is_read2=is_read2)
    f.close()
    g.close()

def quality_trim_pe(inputf_read1, outputf_read1,inputf_read2, outputf_read2,quality_base = None, min_qual_score = 10,
                    min_read_len = 30,adapter_seq_read1 = "AGATCGGAAGAGCACACGTCTGAAC",
                    adapter_seq_read2 = "AGATCGGAAGAGCGTCGTGTAGGGA",num_procs = 1, input_format = None,
                    error_rate = None, max_adapter_removal = None,overlap_length = None, zero_cap = False,
                    path_to_cutadapt = ""):
    """
    Information from cutadapt documentation:
    input_format:
        Input file format; can be either 'fasta', 'fastq' or 'sra-fastq'. Ignored when reading csfasta/qual files
        (default: auto-detect from file name extension).

    inputf_read1,inputf_read2:
        list of filenames for read 1 and read 2 respectively

    outputf_read1,outputf_read2:
        Write the modified sequences to these files instead of standard output and send the summary report to
        standard output. The format is FASTQ if qualities are available, FASTA otherwise. outputf_read1 and outputf_read2
        specify the output filenames of read 1 and read 2 respectively.

    adapter_seq_read1:
        Sequence of an adapter that was ligated to the 3' end of read 1. The adapter itself and anything that follows is
        trimmed.

    adapter_seq_read2:
        Sequence of an adapter that was ligated to the 3' end of read 2. The adapter itself and anything that follows is
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
                 
    if not isinstance(inputf_read1, list):
        if isinstance(inputf_read1, basestring):
            inputf = [inputf_read1]
        else:
            sys.exit("inputf_read1 must be a list of strings")
    if not isinstance(inputf_read2, list):
        if isinstance(inputf_read2, basestring):
            inputf = [inputf_read2]
        else:
            sys.exit("inputf_read2 must be a list of strings")

    if not isinstance(outputf_read1, list):
        if isinstance(outputf_read1, basestring):
            output = [outputf_read1]
        else:
            sys.exit("outputf_read1 must be a list of strings")
    if not isinstance(outputf_read2, list):
        if isinstance(outputf_read2, basestring):
            output = [outputf_read2]
        else:
            sys.exit("outputf_read2 must be a list of strings")            
            
    if len(outputf_read1) != len(inputf_read2) or len(outputf_read1) != len(outputf_read1) or len(outputf_read1) != len(outputf_read2):
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
    options += " -a " + adapter_seq_read1
    options += " -A " + adapter_seq_read2
    options += " " + zero
    if num_procs > 1:
        pool = multiprocessing.Pool(num_procs)
        #adapter trimming
        for current_input_read1,current_output_read1,current_input_read2,current_output_read2 in zip(inputf_read1,outputf_read1,inputf_read2,outputf_read2):
            options += " -o " + current_output_read1 + " " + " -p " + current_output_read2 + " "
            pool.apply_async(subprocess.check_call,(base_cmd + options + current_input_read1 + " " + current_input_read2,),{"shell":True})
        pool.close()
        pool.join()
    else:
        for current_input_read1,current_output_read1,current_input_read2,current_output_read2 in zip(inputf_read1,outputf_read1,inputf_read2,outputf_read2):
            options += " -o " + current_output_read1 + " " + " -p " + current_output_read2 + " "
        subprocess.check_call(base_cmd + options + current_input_read1 + " " + current_input_read2,shell = True)

def flip_read2_strand(input_file,output_file,path_to_samtools=""):
    """
    This function flips the strand of all read2's of mapped paired-end
        reads in input bam file
    
    input_file:
        Input bam file storing the mapped paired-end reads

    output_file:
        Output bam file storing the paired-end reads with strand of read 2 flipped

    path_to_samtools:
        A string of the directory where samtools executive is available. By default,
        samtools is assumed to be included in your path (PATH environmental vairable).
    """
        
    if len(path_to_samtools) > 0:
        path_to_samtools += "/"
    # Input initialization
    input_pipe = subprocess.Popen(
        shlex.split(path_to_samtools+"samtools view -h "+input_file),
        stdout=subprocess.PIPE)
    # Output initialization
    output_handle = open(output_file,'w')
    output_pipe = subprocess.Popen(
        shlex.split(path_to_samtools+"samtools view -S -b -"),
        stdin=subprocess.PIPE,stdout=output_handle)

    for line in input_pipe.stdout:
        # header
        if line[0] == "@":
            output_pipe.stdin.write(line)
            continue

        fields = line.split("\t")
        flag = int(fields[1])
        # Check if it is read 2
        if( (flag & 128) == 0 ): #Not read 2
            output_pipe.stdin.write(line)
            continue
        # flip the strand of read 2
        if( (flag & 16) == 0):
            flag += 16
        else:
            flag -= 16
        fields[1] = str(flag)
        # Write
        output_pipe.stdin.write("\t".join(fields))
        
    # End
    output_handle.close()
    output_pipe.stdin.close()

    
def call_methylated_sites_pe(inputf, sample, reference_fasta, control,sig_cutoff=.01,num_procs = 1,
                             num_upstr_bases=0,num_downstr_bases=2,
                             generate_mpileup_file=True,compress_output=True,
                             split_allc_files=False,
                             min_cov=1,binom_test=True,min_mc=0,path_to_samtools="",sort_mem="500M",bh=True,
                             path_to_files="",min_base_quality=1):

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
    
    sort_mem is the parameter to pass to unix sort with -S/--buffer-size command
    
    bh is a True/False flag indicating whether or not you'd like to use the benjamini-hochberg FDR
        instead of an FDR calculated from the control reference
    
    path_to_files is a string indicating the path for the output and the input bam, mpileup, or allc files
        for methylation calling.
    min_base_quality is an integer indicating the minimum PHRED quality score for a base to be included in the
        mpileup file (and subsequently to be considered for methylation calling)
    """

    #Flip the strand of read 2 and create a new bam file
    ##input: sample+"_processed_reads_no_clonal.bam"
    ##output: sample+"_processed_reads_no_clonal_flipped.bam"

    print_checkpoint("Begin flipping the strand of read 2")
    flip_read2_strand(input_file = inputf,
                   output_file = inputf+".read2flipped.bam",
                   path_to_samtools=path_to_samtools)

        
    #Call methylated sites
    call_methylated_sites(inputf = inputf+".read2flipped.bam",
                          sample = sample,
                          reference_fasta = reference_fasta,
                          control = control,
                          sig_cutoff = sig_cutoff,
                          num_procs = num_procs,
                          num_upstr_bases=num_upstr_bases,
                          num_downstr_bases=num_downstr_bases,
                          generate_mpileup_file=generate_mpileup_file,
                          compress_output=compress_output,
                          split_allc_files=split_allc_files,
                          min_cov = min_cov,
                          binom_test = binom_test,
                          min_mc = min_mc,
                          path_to_samtools = path_to_samtools,
                          sort_mem = sort_mem,
                          bh = bh,
                          path_to_files = path_to_files,
                          min_base_quality = min_base_quality)

    #Remove intermediate bam file
    try:
        subprocess.check_call(shlex.split("rm -f "+inputf+".read2flipped.bam"+
                                          " "+inputf+".read2flipped.bam.bai"))
    except:
        pass

def expand_input_files(read_files,libraries):
    expanded_read_file_list = []
    expanded_library_list = []
    for path,library in zip(read_files,libraries):
        # Assumption: matched read 1 and read2 are stored in files with similar filenames
        # Need to add sort because glob returns files in arbitrary order
        glob_list = glob.glob(path)
        glob_list.sort()
        for filen in glob_list:
            expanded_read_file_list.append(filen)
            expanded_library_list.append(library)
    return(expanded_read_file_list,expanded_library_list)

def merge_bam_files(input_files,output,path_to_samtools=""):
    """
    This function will merge several bam files and create the correct header.
    
    input_files is a list of files produced by collapse_clonal reads. In other words,
        they're assumed to be named like <sample>_<processed_reads>_<lib_id>_no_clonal.bam
    
    output is the name of the merged bam file
    
    path_to_samtools is a string indicating the path to the directory containing your 
        installation of samtools. Samtools is assumed to be in your path if this is not
        provided.
    """
    
    f=open("header.sam",'w')
    subprocess.check_call(
        shlex.split(path_to_samtools+"samtools view -H "+input_files[0]),
        stdout=f
    )

    ## Header
    for filen in input_files:
        f.write("@RG\tID:" + filen[:filen.rindex(".bam")] + "\tLB:" +
                filen +
                "\tSM:NA" + "\n")
    f.close()
    
    subprocess.check_call(
        shlex.split(path_to_samtools+"samtools merge -r -h header.sam "+
                    output +" "+" ".join(input_files))
    )
    subprocess.check_call(["rm", "header.sam"])


def encode_c_positions(seq):
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
    index = seq.find("C",prev_index)
    offset = index + 34
    while True:
        if index < 0:
            break
        while offset > 255:
            indexes += chr(255)
            offset -= 255
        indexes += chr(offset)
        
        prev_index = index + 1
        index = seq.find("C",prev_index)
        offset = index - prev_index + 34
    return indexes
def decode_c_positions(seq,indexes,strand):
    """
    This function takes the encodings generated by encode_c_position and replaces the appropriate
    positions with C nucleotides.
    
    seq is a string of nucleotides to have Cs replaced.
    
    indexes is a string of characters indicating the offsets for the positions of the Cs.
    
    strand is the DNA strand (+ or -) that seq mapped to. This is important because
        sequences in sam files are always represented on the forward strand
    """
    
    prev_index = 0
    new_seq=""
    index = 0
    if strand == "-":
        seq = seq[::-1]
    for char in indexes:
        offset = ord(char)-34
        while offset == 255:
            index+=offset
            offset=ord(char)-34
        index += offset
        if strand == "+":
            new_seq += seq[prev_index:index]+"C"
        elif strand == "-":
            new_seq += seq[prev_index:index]+"G"
        prev_index = index + 1
        index = prev_index
    new_seq += seq[prev_index:]
    if strand == "-":
        new_seq = new_seq[::-1]
    return new_seq

def encode_converted_positions(seq,is_read2=False):
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
            while offset > 255:
                indexes += chr(255)
                offset -= 255
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
            while offset > 255:
                indexes += chr(255)
                offset -= 255
            indexes += chr(offset)                
            prev_index = index + 1
            index = seq.find("G",prev_index)
            offset = index - prev_index + 34
    return indexes

def decode_converted_positions(seq,indexes,strand,is_read2=False):
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
    if is_read2 == False:
        if strand == "-":
            seq = seq[::-1]
        for char in indexes:
            offset = ord(char)-34
            while offset == 255:
                index+=offset
                offset=ord(char)-34
            index += offset
            if strand == "+":
                new_seq += seq[prev_index:index]+"C"
            elif strand == "-":
                new_seq += seq[prev_index:index]+"G"
            prev_index = index + 1
            index = prev_index
    else:
        if strand == "-":
            seq = seq[::-1]
        for char in indexes:
            offset = ord(char)-34
            while offset == 255:
                index+=offset
                offset=ord(char)-34
            index += offset
            if strand == "+":
                new_seq += seq[prev_index:index]+"G"
            elif strand == "-":
                new_seq += seq[prev_index:index]+"C"
            prev_index = index + 1
            index = prev_index
    new_seq += seq[prev_index:]
    if strand == "-":
        new_seq = new_seq[::-1]
    return new_seq
 
def find_multi_mappers(inputf,output,num_procs=1,keep_temp_files=False,append=False):
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
    sam_header = []
    file_handles = {}
    f = open(inputf,'r')
    cycle = itertools.cycle(range(0,num_procs))
    for file_num in xrange(0,num_procs):
        if append == False:
            file_handles[file_num]=open(output+"_sorted_"+str(file_num),'w')
        else:
            file_handles[file_num]=open(output+"_sorted_"+str(file_num),'a')
    for line in f:
        #To deal with the way chromosomes were named in some of our older references
        if line[0] == "@":
            continue

        fields = line.split("\t")
        fields[2] = fields[2].replace("_f","")
        fields[2] = fields[2].replace("_r","")
        if fields[2] != "*":
            header = fields[0].split("!")
            #BIG ASSUMPTION!! NO TABS IN FASTQ HEADER LINES EXCEPT THE ONES I ADD!
            if (int(fields[1]) & 16) == 16:
                strand = "-"
            elif (int(fields[1]) & 16) == 0:
                strand = "+"
            seq = decode_c_positions(fields[9],header[-1],strand)
            file_handles[cycle.next()].write(" ".join(header[:-1])+"\t"+"\t".join(fields[1:9])+"\t"+seq+"\t"+"\t".join(fields[10:]))
    f.close()
    if keep_temp_files == False:
        subprocess.check_call(shlex.split("rm "+inputf))
    for file_num in xrange(0,num_procs):
        file_handles[file_num].close()
               
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
    for first in ["A","T","C","G","N"]:
        for second in ["A","T","C","G","N"]:
            test_num["C"+first+second] = 1
            prev_bh_value["C"+first+second] = 0
            best_fdr["C"+first+second] = 0
            best_pvalue["C"+first+second] = 1
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
        bh_value = float(fields[6]) * mc_class_counts[fields[3]] / (test_num[fields[3]] + 1)
        # Sometimes this correction can give values greater than 1,
        # so we set those values at 1
        bh_value = min(bh_value, 1)
        prev_bh_value[fields[3]] = bh_value
        #if bh_value <= sig_cutoff and bh_value >= best_fdr:
        if bh_value <= sig_cutoff:
            best_fdr[fields[3]] = bh_value
            best_pvalue[fields[3]] = float(fields[6])
        
        test_num[fields[3]] += 1
        input_lines[min_pvalue]=input_files[min_pvalue].readline().rstrip()
        input_fields[min_pvalue]=input_lines[min_pvalue].split("\t")
        try:
            input_pvalues[min_pvalue]=float(input_fields[min_pvalue][6])
        except:
            #Dummy value that will never be the minimum
            input_pvalues[min_pvalue]=2.0
        min_pvalue = min(input_pvalues,key=input_pvalues.get)
    for mc_class in best_pvalue:
        print "The closest p-value cutoff for "+mc_class+" at your desired FDR is "+str(best_pvalue[mc_class])+" which corresponds to an FDR of "+str(best_fdr[mc_class])
    for filen in files:
        input_files[filen].close() 
    return best_pvalue
                
def allc_run_binom_tests(filen,non_conversion,min_cov=1,sort_mem="500M"):
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
        if sort_mem.find("-S") == -1:
            sort_mem = " -S " + sort_mem
    else:
        sort_mem = ""
    mc_class_counts = {}
    for first in ["A","T","C","G","N"]:
        for second in ["A","T","C","G","N"]:
            mc_class_counts["C"+first+second]=0
    obs_pvalues = {}
    reverse_complement = {"A":"T","C":"G","G":"C","T":"A","N":"N"}
    f = open(filen,'r')
    g = open(filen+"_binom_results.tsv",'w')
    for line in f:
        line = line.rstrip()
        
        fields = line.split("\t")
        
        mc_class = fields[3]
        try:
            unconverted_c = int(fields[4])
            converted_c = int(fields[5]) - unconverted_c
        except:
            continue
        total = int(fields[5])
        if total >= min_cov and unconverted_c != 0:
            try:
                p_value = obs_pvalues[(unconverted_c,total)]
            except:
                p_value = sci.binom.sf(unconverted_c-1,total,non_conversion)
                obs_pvalues[(unconverted_c,total)] = p_value
            g.write("\t".join(fields[:6])+"\t"+str(p_value)+"\n")
            mc_class_counts[mc_class]+=1
        elif total != 0:
            #a dummy value that will always sort to the bottom of the BH correction and be interpreted as
            #a unmethylated site
            p_value = 2.0
            g.write("\t".join(fields[:6])+"\t"+str(p_value)+"\n")
        
    g.close()
    subprocess.check_call(shlex.split("sort" + sort_mem + " -k 7g,7g -o "+filen+"_binom_results.tsv "+filen+"_binom_results.tsv"))
    return mc_class_counts

def parse_args():
     # create the top-level parser
     parser = ArgumentParser(prog='PROG')
     subparsers = parser.add_subparsers(help='Process all commands', dest='command')

     # create the parser for the "call_mc" command
     parser_pipeline = subparsers.add_parser('run_methylation_pipeline', help='Use to run the methylation pipeline')
     parser_pipeline.add_argument('--files', type=str, nargs="+", required=True, help="list of all the fastq files you'd like to run \
        through the pipeline. Note that globbing is supported here (i.e., you can use * in your paths)")
     parser_pipeline.add_argument('--libraries', type=str, nargs="+", required=True, help="list of library IDs (in the same order as \
        the files list) indiciating which libraries each set of fastq files belong to. If you use a glob, you only need to indicate \
        the library ID for those fastqs once (i.e., the length of files and libraries should be the same)")
     parser_pipeline.add_argument('--sample', type=str, required=True, help="String indicating the name of the sample you're processing. \
        It will be included in the output files.")
     parser_pipeline.add_argument('--forward_ref', type=str, required=True, help="string indicating the path to the forward strand \
        reference created by build_ref")
     parser_pipeline.add_argument('--reverse_ref', type=str, required=True, help="string indicating the path to the reverse strand \
        reference created by build_ref")
     parser_pipeline.add_argument('--ref_fasta', type=str, required=True, help="string indicating the path to a fasta file containing \
        the sequences you used for mapping")
     parser_pipeline.add_argument('--unmethylated_control', type=str, required=True, help="name of the chromosome/region that you want \
        to use to estimate the non-conversion rate of your sample, or the non-conversion rate you'd like to use. Consequently, control \
        is either a string, or a decimal. If control is a string then it should be in the following format: 'chrom:start-end'. \
        If you'd like to specify an entire chromosome simply use 'chrom:'")
     parser_pipeline.add_argument('--path_to_samtools', type=str, default="", help='Path to samtools installation (default is current dir)')     
     parser_pipeline.add_argument('--path_to_aligner', type=str, default="", help='Path to bowtie installation (default is current dir)')
     parser_pipeline.add_argument('--aligner_options', type=str, nargs='+', help="list of strings indicating options you'd like passed to bowtie \
        (e.g., ['-k 1','-l 2']")
     parser_pipeline.add_argument('--num_procs', type=int, default=1, help='Number of processors you wish to use to \
        parallelize this function')
     parser_pipeline.add_argument('--trim_reads', type=bool, default=True, help='Whether to trim reads using cutadapt (default is True)')
     parser_pipeline.add_argument('--path_to_cutadapt', type=str, default="", help='Path to cutadapt installation (default is current dir)')
     parser_pipeline.add_argument('--adapter_seq', type=str, default="AGATCGGAAGAGCACACGTCTG", help="sequence of an adapter that was ligated \
        to the 3' end. The adapter itself and anything that follows is trimmed.")
     parser_pipeline.add_argument('--max_adapter_removal', type=int, help="Indicates the maximum number of times to try to remove adapters. \
        Useful when an adapter gets appended multiple times.")
     parser_pipeline.add_argument('--overlap_length', type=int, help="Minimum overlap length. If the overlap between the read and the adapter \
        is shorter than LENGTH, the read is not modified. This reduces the no. of bases trimmed purely due to short random adapter matches.")
     parser_pipeline.add_argument('--zero_cap', type=bool, help="Flag that causes negative quality values to be set to zero (workaround to avoid \
        segmentation faults in BWA)")
     parser_pipeline.add_argument('--error_rate', type=float, help="maximum allowed error rate (no. of errors divided by the length \
        of the matching region)")
     parser_pipeline.add_argument('--min_qual_score', type=int, default=10, help="allows you to trim low-quality ends from reads before \
        adapter removal. The algorithm is the same as the one used by BWA (Subtract CUTOFF from all qualities; compute partial sums from \
        all indices to the end of the sequence; cut sequence at the index at which the sum is minimal).")
     parser_pipeline.add_argument('--min_read_len', type=int, default=30, help="indicates the minimum length a read must be to be kept. \
        Reads that are too short even before adapter removal are also discarded. In colorspace, an initial primer is not counted.")
     parser_pipeline.add_argument('--sig_cutoff', type=float, default=.01, help="float indicating the adjusted p-value cutoff you wish to \
        use for determining whether or not a site is methylated")
     parser_pipeline.add_argument('--min_cov', type=int, default=0, help="integer indicating the minimum number of reads for a site to be tested.")
     parser_pipeline.add_argument('--binom_test', type=bool, default=False, help="Indicates that you'd like to use a binomial test, rather than the \
        alternative method outlined here: https://bitbucket.org/schultzmattd/methylpy/wiki/Methylation%20Calling")
     parser_pipeline.add_argument('--keep_temp_files', type=bool, default=False, help="Boolean indicating that you'd like to keep the intermediate \
        files generated by this function. This can be useful for debugging, but in general should be left False.")
     parser_pipeline.add_argument('--save_space', type=bool, default=True, help="indicates whether or not you'd like to perform read collapsing \
        right after mapping or once all the libraries have been mapped. If you wait until after everything has been mapped, the collapsing can be \
        parallelized. Otherwise the collapsing will have to be done serially. The trade-off is that you must keep all the mapped files around, rather \
        than deleting them as they are processed, which can take up a considerable amount of space. It's safest to set this to True.")
     parser_pipeline.add_argument('--bowtie2', type=bool, default=False, help="Specifies whether to use the bowtie2 aligner instead of bowtie")
     parser_pipeline.add_argument('--sort_mem', type=str, help="Parameter to pass to unix sort with -S/--buffer-size command")
     parser_pipeline.add_argument('--path_to_output', type=str, default="", help="Path to a directory where you would like the output to be stored. \
        The default is the same directory as the input fastqs.")
     parser_pipeline.add_argument('--path_to_picard', type=str, default=False, help="The path to picard.jar. Default is false indicating that you don't want to use this jar for duplication removal")
     parser_pipeline.add_argument('--remove_clonal', type=bool, default=True, help="Remove clonal reads or not")

     
     #create the parser for the "call_methylated_sites" command
     parser_call = subparsers.add_parser('call_methylated_sites', help='Use to run the call_methylated_sites function')
     parser_call.add_argument('inputf', type=str, help='inputf is the path to a bam file that contains mapped bisulfite sequencing reads')
     parser_call.add_argument('sample', type=str, help="output is the name you'd like for the allc files. The files will be named like so: allc_<sample>_<chrom>.tsv")
     parser_call.add_argument('reference', type=str, help="reference is the path to a samtools indexed fasta file")
     parser_call.add_argument('control', type=str, help="control is the name of the chromosome/region that you want to use to \
        estimate the non-conversion rate of your sample, or the non-conversion rate you'd like to use. Consequently, control \
        is either a string, or a decimal. If control is a string then it should be in the following format: 'chrom:start-end'. \
        If you'd like to specify an entire chromosome simply use 'chrom:'")
     parser_call.add_argument('casava_version', type=float, help="casava_version is a float indicating which version of casava was used to generate the fastq files.")
     parser_call.add_argument('--sig_cutoff', type=float, default=0.01, help="sig_cutoff is a float indicating the adjusted \
        p-value cutoff you wish to use for determining whether or not a site is methylated")
     parser_call.add_argument('--num_procs', type=int, default=1, help="processers is an integer indicating how many processors you'd like to run this function over") 
     parser_call.add_argument('--min_cov', type=int, default=1, help="min_cov is an integer indicating the minimum number of reads for a site to be tested")
     parser_call.add_argument('--binom_test', type=bool, default=False, help="Boolean indicating if you want to run binomial tests")
     parser_call.add_argument('--min_mc', type=int, default=0, help="Minimum number of mCs that must be observed")
     parser_call.add_argument('--path_to_samtools', type=str, default="", help='Path to samtools installation (default is current dir)')
     parser_call.add_argument('--sort_mem', type=str, default=False, help="Parameter to pass to unix sort with -S/--buffer-size command")
     parser_call.add_argument('--bh', type=bool, default=False, help="Boolean flag indicating whether or not you'd like to use the benjamini-hochberg FDR \
        instead of an FDR calculated from the control reference")
     parser_call.add_argument('--path_to_files', type=str, default="", help="string indicating the path for the output and the input bam, mpileup, or allc files \
        for methylation calling.")                                                                                                                                                   
     
     args = parser.parse_args()

     if args.command == "run_methylation_pipeline":
         if not args.aligner_options:
             args.aligner_options = ["-S","-k 1","-m 1","--chunkmbs 3072","--best","--strata","-o 4","-e 80","-l 20","-n 0"]
        
         run_methylation_pipeline_pe(args.files,args.libraries,args.sample,args.forward_ref,args.reverse_ref,args.ref_fasta,
                                     args.unmethylated_control,args.path_to_samtools,args.path_to_aligner,
                                     args.aligner_options,args.num_procs,args.trim_reads,
                                     args.path_to_cutadapt,args.adapter_seq,args.max_adapter_removal,
                                     args.overlap_length,args.zero_cap,args.error_rate,args.min_qual_score,args.min_read_len,
                                     args.sig_cutoff,args.min_cov,args.binom_test,args.keep_temp_files,
                                     args.save_space,
                                     args.bowtie2,args.sort_mem,args.path_to_output)
                                  
     elif args.command == "call_methylated_sites":
         call_methylated_sites(args.inputf, args.sample, args.reference, args.control, args.casava_version, args.sig_cutoff,
                              args.num_procs, args.min_cov, args.binom_test, args.min_mc, args.path_to_samtools, args.sort_mem,
                              args.bh, args.path_to_files)
    
if __name__ == '__main__':
    parse_args()
