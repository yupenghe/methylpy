import sys
import os
import multiprocessing
import subprocess
import scipy.stats as sci
from scipy.stats.mstats import mquantiles
from methylpy.utilities import print_checkpoint, print_error, split_fastq_file, print_warning
import pdb
import shlex
import itertools
import re
import glob
import io as cStr
import bisect
from methylpy.call_mc_se import call_methylated_sites, remove_clonal_bam
from methylpy.call_mc_se import encode_c_positions,decode_c_positions
from methylpy.utilities import check_call_mc_dependencies

def run_methylation_pipeline_pe(read1_files, read2_files, sample,
                                forward_reference, reverse_reference, reference_fasta,
                                libraries = None,
                                unmethylated_control=None,
                                path_to_output="", sig_cutoff=0.01,
                                num_procs=1, sort_mem="500M",
                                num_upstr_bases=0,
                                num_downstr_bases=2,
                                generate_allc_file=True,
                                generate_mpileup_file=True,
                                compress_output=True,
                                bgzip=False,
                                path_to_bgzip="",
                                path_to_tabix="",
                                binom_test=True, min_cov=2,
                                trim_reads=True, path_to_cutadapt="",
                                aligner="bowtie2",
                                path_to_aligner="", aligner_options=None,
                                merge_by_max_mapq=False,min_mapq=30,
                                pbat=False,check_dependency=True,
                                remove_clonal=False, keep_clonal_stats=True,
                                path_to_picard="",java_options="-Xmx20g",
                                path_to_samtools="",
                                remove_chr_prefix=True,
                                add_snp_info=False,
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
    if check_dependency:
        check_call_mc_dependencies(path_to_samtools,
                                   trim_reads,
                                   path_to_cutadapt,
                                   aligner,
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

    if len(libraries) == 1 and len(read1_files) > 1:
        uniq_library = libraries[0]
        libraries = [uniq_library for ind in range(len(read1_files))]

    #Default bowtie option
    if aligner_options is None:
        if aligner.lower() == "minimap2":
            aligner_options = ["-ax","sr","--secondary=no"]
        elif aligner.lower() == "bowtie":
            aligner_options = ["-X 1000", "-S", "-k 1", "-m 1", "--best", "--strata",
                               "--chunkmbs 3072", "-n 1", "-e 100"]
            aligner_options.append("--phred33-quals")
        else: # bowtie2
            aligner_options = ["-X 1000", "--no-discordant", "--no-mixed"]
            aligner_options.append("--phred33-quals")

    # CASAVA >= 1.8
    quality_base = 33

    #just to avoid any paths missing a slash. It's ok to add an extra slash if
    #there's already one at the end
    if len(path_to_samtools) != 0:
        path_to_samtools += "/"

    if len(path_to_aligner) != 0:
        path_to_aligner += "/"
    if len(path_to_output) != 0:
        path_to_output += "/"
        if not os.path.exists(path_to_output):
            try:
                os.makedirs(path_to_output)
            except:
                print_error("  Failed to create output folder!")

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
            aligner=aligner,
            aligner_options=aligner_options,
            merge_by_max_mapq=merge_by_max_mapq,
            min_mapq=min_mapq,
            num_procs=num_procs,
            trim_reads=trim_reads, path_to_cutadapt=path_to_cutadapt,
            adapter_seq_read1=adapter_seq_read1, adapter_seq_read2=adapter_seq_read2,
            max_adapter_removal=max_adapter_removal, overlap_length=overlap_length,
            zero_cap=zero_cap, quality_base=quality_base, error_rate=error_rate,
            min_qual_score=min_qual_score, min_read_len=min_read_len,
            keep_temp_files=keep_temp_files,
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

    print_checkpoint("There are " + str(total_input) + " total input read pairs")
    print_checkpoint("There are " + str(total_unique) + " uniquely mapping read pairs, " +
                     str(float(total_unique) / total_input*100) + " percent remaining")

    if remove_clonal == True:
        total_non_clonal = total_unique - total_clonal
        print_checkpoint("There are "+str(total_non_clonal)+" non-clonal read pairs, "+
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
    if generate_allc_file:
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
                                 bgzip=bgzip,
                                 path_to_bgzip=path_to_bgzip,
                                 path_to_tabix=path_to_tabix,
                                 min_cov=min_cov,
                                 binom_test=binom_test,
                                 remove_chr_prefix=remove_chr_prefix,
                                 sort_mem=sort_mem,
                                 path_to_files=path_to_output,
                                 path_to_samtools=path_to_samtools,
                                 min_base_quality=min_base_quality,
                                 keep_temp_files=keep_temp_files)
    print_checkpoint("Done")

def run_mapping_pe(current_library, library_read1_files, library_read2_files,
                   sample, forward_reference, reverse_reference, reference_fasta,
                   path_to_output="",
                   path_to_samtools="", path_to_aligner="",
                   aligner="bowtie2",
                   aligner_options=None,
                   merge_by_max_mapq=False,
                   min_mapq=30,
                   num_procs=1, trim_reads=True, path_to_cutadapt="",
                   adapter_seq_read1="AGATCGGAAGAGCACACGTCTGAAC",
                   adapter_seq_read2="AGATCGGAAGAGCGTCGTGTAGGGA",
                   max_adapter_removal=None, overlap_length=None, zero_cap=None,
                   quality_base=None, error_rate=None, min_qual_score=10,
                   min_read_len=30, keep_temp_files=False,
                   sort_mem="500M"):
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
        
    sort_mem is the parameter to pass to unix sort with -S/--buffer-size command
    """

    #Default option
    if aligner_options is None:
        if aligner.lower() == "minimap2":
            aligner_options = ["-ax","sr","--secondary=no"]
        elif aligner.lower() == "bowtie":
            aligner_options = ["-X 1000", "-S", "-k 1", "-m 1", "--best", "--strata",
                               "--chunkmbs 3072", "-n 1", "-e 100"]
            aligner_options.append("--phred33-quals")
        else: # bowtie2
            aligner_options = []
            aligner_options = ["-X 1000", "--no-discordant", "--no-mixed"]
            aligner_options.append("--phred33-quals")

    # CASAVA >= 1.8
    quality_base = 33

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
            inputf_read1=[prefix+"_read1_split_"+str(i) for i in range(0,num_procs)],
            outputf_read1=[prefix+"_read1_split_trimmed_"+str(i) for i in range(0,num_procs)],
            inputf_read2=[prefix+"_read2_split_"+str(i) for i in range(0,num_procs)],
            outputf_read2=[prefix+"_read2_split_trimmed_"+str(i) for i in range(0,num_procs)],
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
        
        subprocess.check_call(shlex.split("rm "+" ".join([prefix+"_read1_split_"+str(i) for i in range(0,num_procs)])))
        subprocess.check_call(shlex.split("rm "+" ".join([prefix+"_read2_split_"+str(i) for i in range(0,num_procs)])))
        
        #Conversion
        print_checkpoint("Begin converting reads for "+str(current_library))
        if num_procs > 1:
            pool = multiprocessing.Pool(num_procs)#read1
            for inputf,output in zip([prefix+"_read1_split_trimmed_"+str(i) for i in range(0,num_procs)],
                                     [prefix+"_read1_split_trimmed_converted_"+str(i) for i in range(0,num_procs)]):
                pool.apply_async(convert_reads_pe,(inputf,output))
            for inputf,output in zip([prefix+"_read2_split_trimmed_"+str(i) for i in range(0,num_procs)],
                                     [prefix+"_read2_split_trimmed_converted_"+str(i) for i in range(0,num_procs)]):
                pool.apply_async(convert_reads_pe,(inputf,output,True))
            pool.close()
            pool.join()
        else:
            for inputf,output in zip([prefix+"_read1_split_trimmed_"+str(i) for i in range(0,num_procs)],
                                     [prefix+"_read1_split_trimmed_converted_"+str(i) for i in range(0,num_procs)]):
                convert_reads_pe(inputf,output)
            for inputf,output in zip([prefix+"_read2_split_trimmed_"+str(i) for i in range(0,num_procs)],
                                     [prefix+"_read2_split_trimmed_converted_"+str(i) for i in range(0,num_procs)]):
                convert_reads_pe(inputf,output,True)
                
        subprocess.check_call(
            shlex.split("rm "+" ".join([prefix+"_read1_split_trimmed_"+str(i) for i in range(0,num_procs)]))
        )
        subprocess.check_call(
            shlex.split("rm "+" ".join([prefix+"_read2_split_trimmed_"+str(i) for i in range(0,num_procs)]))
        )        
        #Run bowtie
        input_fastq_read1 = [prefix+"_read1_split_trimmed_converted_"+str(i) for i in range(0,num_procs)]
        input_fastq_read2 = [prefix+"_read2_split_trimmed_converted_"+str(i) for i in range(0,num_procs)]
    else:
        print_checkpoint("No trimming applied on reads")  
        #Conversion
        print_checkpoint("Begin converting reads for "+str(current_library))
        if num_procs > 1:
            pool = multiprocessing.Pool(num_procs)#read1
            for inputf,output in zip([prefix+"_read1_split_"+str(i) for i in range(0,num_procs)],
                                     [prefix+"_read1_split_converted_"+str(i) for i in range(0,num_procs)]):
                pool.apply_async(convert_reads_pe,(inputf,output))
            for inputf,output in zip([prefix+"_read2_split_"+str(i) for i in range(0,num_procs)],
                                     [prefix+"_read2_split_converted_"+str(i) for i in range(0,num_procs)]):
                pool.apply_async(convert_reads_pe,(inputf,output,True))
            pool.close()
            pool.join()
        else:
            for inputf,output in zip([prefix+"_read1_split_"+str(i) for i in range(0,num_procs)],
                                     [prefix+"_read1_split_converted_"+str(i) for i in range(0,num_procs)]):
                convert_reads_pe(inputf,output)
            for inputf,output in zip([prefix+"_read2_split_"+str(i) for i in range(0,num_procs)],
                                     [prefix+"_read2_split_converted_"+str(i) for i in range(0,num_procs)]):
                convert_reads_pe(inputf,output,True)
        subprocess.check_call(shlex.split("rm "+" ".join([prefix+"_read1_split_"+str(i) for i in range(0,num_procs)])))
        subprocess.check_call(shlex.split("rm "+" ".join([prefix+"_read2_split_"+str(i) for i in range(0,num_procs)])))
        input_fastq_read1 = [prefix+"_read1_split_converted_"+str(i) for i in range(0,num_procs)]
        input_fastq_read2 = [prefix+"_read2_split_converted_"+str(i) for i in range(0,num_procs)]
                    
    #Run bowtie
    if aligner.lower() == "minimap2":
        print_checkpoint("Begin Running minimap2 for "+current_library)
    elif aligner.lower() == "Bowtie":
        print_checkpoint("Begin Running Bowtie for "+current_library)
    else:
        print_checkpoint("Begin Running Bowtie2 for "+current_library)
    total_unique = run_alignment_pe(current_library,
                                    input_fastq_read1,
                                    input_fastq_read2,
                                    sample,
                                    forward_reference,reverse_reference,reference_fasta,
                                    path_to_output=path_to_output,
                                    path_to_samtools=path_to_samtools,
                                    aligner=aligner,
                                    aligner_options=aligner_options,
                                    merge_by_max_mapq=merge_by_max_mapq,
                                    min_mapq=min_mapq,
                                    path_to_aligner=path_to_aligner,num_procs=num_procs,
                                    keep_temp_files=keep_temp_files,
                                    sort_mem=sort_mem)
    
    return total_input,total_unique

def run_alignment_pe(current_library,library_read1_files,library_read2_files,
                     sample,
                     forward_reference,reverse_reference,reference_fasta,
                     path_to_output="",                  
                     path_to_samtools="",
                     aligner="bowtie2",
                     aligner_options=None,path_to_aligner="",
                     merge_by_max_mapq=False,min_mapq=30,
                     num_procs=1,keep_temp_files=False,sort_mem="500M"):
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

    if not sort_mem:
        sort_option = ""
    else:
        sort_option = " -S "+sort_mem

    if len(path_to_output) !=0:
        path_to_output+="/"

    #Default bowtie option
    if aligner_options is None:
        if aligner.lower() == "minimap2":
            aligner_options = ["-ax","sr","--secondary=no"]
        elif aligner.lower() == "bowtie":
            aligner_options = ["-X 1000", "-S", "-k 1", "-m 1", "--best", "--strata",
                               "--chunkmbs 3072", "-n 1", "-e 100"]
            aligner_options.append("--phred33-quals")
        else: # bowtie2
            aligner_options = ["-X 1000", "--no-discordant", "--no-mixed"]
            aligner_options.append("--phred33-quals")

    # CASAVA >= 1.8
    quality_base = 33

    options = aligner_options
    
    if aligner != "minimap2":
        if " ".join(options).find(" -p ") == -1:
            options.append("-p "+str(num_procs))
    else:
        if " ".join(options).find(" -t ") == -1:
            options.append("-t "+str(num_procs))

    prefix = path_to_output+sample+"_"+str(current_library)


    input_read1_file = ",".join(library_read1_files)
    input_read2_file = ",".join(library_read2_files)
    input_read1_file = prefix+"_converted_read1s.fastq"
    input_read2_file = prefix+"_converted_read2s.fastq"
    subprocess.check_call(["mv",library_read1_files[0],input_read1_file])
    subprocess.check_call(["mv",library_read2_files[0],input_read2_file])
    if len(library_read1_files) > 1:
        with open(input_read1_file,'a') as g:
            for library_read1_file in library_read1_files[1:]:
                with open(library_read1_file,'r') as f:
                    g.write(f.read())
                subprocess.check_call(["rm",library_read1_file])
    if len(library_read2_files) > 1:
        with open(input_read2_file,'a') as g:
            for library_read2_file in library_read2_files[1:]:
                with open(library_read2_file,'r') as f:
                    g.write(f.read())
                subprocess.check_call(["rm",library_read2_file])

                    
    ## Forward
    if aligner.lower() == "minimap2":
        args = [path_to_aligner+"minimap2"]
        args.extend(options)
        args.append("--for-only")
        args.append(forward_reference)
        args.append(input_read1_file)
        args.append(input_read2_file)
    elif aligner.lower() == "bowtie":
        args = [path_to_aligner+"bowtie"]
        args.extend(options)
        args.append("--norc")
        args.append(forward_reference)
        args.append("-1 "+input_read1_file)
        args.append("-2 "+input_read2_file)
    else: # bowtie2
        args = [path_to_aligner+"bowtie2"]
        args.extend(options)
        args.append("--norc")
        args.append("-x "+forward_reference)
        args.append("-1 "+input_read1_file)
        args.append("-2 "+input_read2_file)
    ## run
    with open(prefix+"_forward_strand_hits.sam","w") as f:
        subprocess.check_call(shlex.split(" ".join(args)),stdout=f)
    print_checkpoint("Processing forward strand hits")
    find_multi_mappers_pe(prefix+"_forward_strand_hits.sam",
                          prefix,
                          num_procs=num_procs,
                          min_mapq=min_mapq,
                          append=False,
                          keep_temp_files=keep_temp_files)
    
    ## Reverse    
    if aligner.lower() == "minimap2":
        args = [path_to_aligner+"minimap2"]
        args.extend(options)
        args.append("--rev-only")
        args.append(reverse_reference)
        args.append(input_read1_file)
        args.append(input_read2_file)
    elif aligner.lower() == "bowtie":
        args = [path_to_aligner+"bowtie"]
        args.extend(options)
        args.append("--nofw")
        args.append(reverse_reference)
        args.append("-1 "+input_read1_file)
        args.append("-2 "+input_read2_file)
    else: # bowtie2
        args = [path_to_aligner+"bowtie2"]
        args.extend(options)
        args.append("--nofw")
        args.append("-x "+reverse_reference)
        args.append("-1 "+input_read1_file)
        args.append("-2 "+input_read2_file)
    ## run
    with open(prefix+"_reverse_strand_hits.sam","w") as f:
        subprocess.check_call(shlex.split(" ".join(args)),stdout=f)
    subprocess.check_call(["rm",input_read1_file,input_read2_file])
    print_checkpoint("Processing reverse strand hits")
    find_multi_mappers_pe(prefix+"_reverse_strand_hits.sam",
                          prefix,
                          num_procs=num_procs,
                          min_mapq=min_mapq,
                          append=True,
                          keep_temp_files=keep_temp_files)

    if num_procs > 1:
        pool = multiprocessing.Pool(num_procs)
        for file_num in range(0,num_procs):
            pool.apply_async(subprocess.check_call,(
                shlex.split("env LC_COLLATE=C sort" + sort_option +
                            " -t '\t' -k 1 -o "+prefix+"_sorted_"+str(file_num)+" "+
                            prefix+"_sorted_"+str(file_num)),))
        pool.close()
        pool.join()
    else:
        for file_num in range(0,num_procs):
            subprocess.check_call(shlex.split("env LC_COLLATE=C sort" + sort_option +
                                              " -t '\t' -k 1 -o "+
                                              prefix+"_sorted_"+str(file_num)+" "+
                                              prefix+"_sorted_"+str(file_num)))
    print_checkpoint("Finding multimappers")

    if merge_by_max_mapq:
        total_unique = merge_sorted_multimap_pe_max_mapq(
            current_library,
            [prefix+"_sorted_"+str(file_num) for file_num in range(0,num_procs)],
            prefix,
            reference_fasta,
            path_to_samtools="")
    else:
        total_unique = merge_sorted_multimap_pe(
            current_library,
            [prefix+"_sorted_"+str(file_num) for file_num in range(0,num_procs)],
            prefix,
            reference_fasta,
            path_to_samtools="")

    subprocess.check_call(shlex.split("rm "+" ".join([prefix+"_sorted_"+str(file_num)
                                                      for file_num in range(0,num_procs)])))
    
    if not sort_mem:
        sort_option = ""
    else:
        sort_option = " -m "+sort_mem
        
    output_bam_file = prefix+"_processed_reads.bam"
    try:
        subprocess.check_call(shlex.split(path_to_samtools+"samtools sort "+
                                          " -@ " + str(num_procs) +
                                          sort_option + " " +
                                          " -o "+output_bam_file + " " + 
                                          output_bam_file))
    except:
        subprocess.check_call(shlex.split(path_to_samtools+"samtools sort "+
                                          " -o "+output_bam_file + " " +
                                          output_bam_file ))
    return total_unique

def find_multi_mappers_pe(inputf,output,num_procs=1,
                          min_mapq=30,
                          keep_temp_files=False,append=False):
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
        if line[0] == "@":
            continue
        
        fields = line.split("\t")

        ## Check if it is proper pair
        flag = int(fields[1])
        if (flag & 2) == 0 or int(fields[4]) < min_mapq or (flag & 2048) == 2048:
            continue;
        
        header = fields[0].split("!")
        if (flag & 16) == 16:
            strand = "-"
        else:
            strand = "+"
        if (flag & 128) == 128:
            is_read2 = True
        else:
            is_read2 = False
        try:
            seq = decode_c_positions(fields[9],header[-1],strand,is_read2)
            file_handles[next(cycle)].write(" ".join(header[:-1])+"\t"+"\t".join(fields[1:9])+
                                            "\t"+seq+"\t"+"\t".join(fields[10:]))
        except:
            print_warning("  Failed to recover unconverted sequence for:\n"+line+"\n")
            print_warning(header[-1]+"\n")
    f.close()
    if keep_temp_files == False:
        subprocess.check_call(shlex.split("rm "+inputf))
        pass
    for file_num in range(0,num_procs):
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
    for index,filen in enumerate(files):
        file_handles[filen]=open(filen,'r')
        lines[filen]=file_handles[filen].readline()
        fields[filen] = lines[filen].split("\t")[0]#Read ID
    while True:
        all_fields = [field for field in list(fields.values()) if field != ""]
        if len(all_fields) == 0:
            break
        min_field = min(all_fields)
        count_1, count_2 = 0, 0
        current_line_1, current_line_2 = "", ""
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
    return total_unique

def merge_sorted_multimap_pe_max_mapq(current_library,files,prefix,reference_fasta,path_to_samtools=""):
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
        print("Reference fasta not indexed. Indexing.")
        try:
            subprocess.check_call(shlex.split(path_to_samtools+"samtools faidx "+reference_fasta))
            f = open(reference_fasta+".fai",'r')
        except:
            sys.exit("Reference fasta wasn't indexed, and couldn't be indexed. "
                     +"Please try indexing it manually and running methylpy again.")
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
    for index,filen in enumerate(files):
        file_handles[filen]=open(filen,'r')
        lines[filen]=file_handles[filen].readline()
        fields[filen] = lines[filen].split("\t")[0]#Read ID
    while True:
        all_fields = [field for field in list(fields.values()) if field != ""]
        if len(all_fields) == 0:
            break
        min_field = min(all_fields)
        count_1, count_2 = 0, 0
        current_line_1, current_line_2 = "", ""
        count= 0
        max_mapq, min_mapq = -1000,1000
        for key in fields:
            while fields[key] == min_field:
                count += 1
                if int(lines[key].split("\t")[4]) >= max_mapq:
                    max_mapq = int(lines[key].split("\t")[4])
                    if(int(lines[key].split("\t")[1]) & 64 == 64): #First in pair
                        #count_1 += 1
                        current_line_1 = lines[key]
                    else:
                        #count_2 += 1
                        current_line_2 = lines[key]

                if int(lines[key].split("\t")[4]) < min_mapq:
                    min_mapq = int(lines[key].split("\t")[4])

                lines[key]=file_handles[key].readline()
                fields[key]=lines[key].split("\t")[0]
        #Check if there is only one valid alignment
        #if count_1 == 1:
        if count == 2 or max_mapq > min_mapq:
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
    return total_unique

def convert_reads_pe(inputf,output,is_read2=False,buffer_line_number=100000):
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
    encoding = encode_c_positions(seq,is_read2=is_read2)
    
    line_counts = 0
    out = ""
    if is_read2 == False:
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
            encoding = encode_c_positions(seq,is_read2=is_read2)
    else:
        while header:
            out += header+"!"+encoding+"\n"
            converted_seq = seq.replace("G","A")
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
            encoding = encode_c_positions(seq,is_read2=is_read2)

    # output
    if line_counts > 0:
        g.write(out)
        line_counts = 0
        out = ""
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
        if isinstance(inputf_read1, str):
            inputf = [inputf_read1]
        else:
            sys.exit("inputf_read1 must be a list of strings")
    if not isinstance(inputf_read2, list):
        if isinstance(inputf_read2, str):
            inputf = [inputf_read2]
        else:
            sys.exit("inputf_read2 must be a list of strings")

    if not isinstance(outputf_read1, list):
        if isinstance(outputf_read1, str):
            output = [outputf_read1]
        else:
            sys.exit("outputf_read1 must be a list of strings")
    if not isinstance(outputf_read2, list):
        if isinstance(outputf_read2, str):
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
        stdout=subprocess.PIPE,
        universal_newlines=True)
    # Output initialization
    output_handle = open(output_file,'w')
    output_pipe = subprocess.Popen(
        shlex.split(path_to_samtools+"samtools view -S -b -"),
        stdin=subprocess.PIPE,stdout=output_handle,
        universal_newlines=True)

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

    
def call_methylated_sites_pe(inputf, sample, reference_fasta,
                             unmethylated_control = None,
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
                             add_snp_info=False,
                             sort_mem="500M",
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
                          unmethylated_control = unmethylated_control,
                          sig_cutoff = sig_cutoff,
                          num_procs = num_procs,
                          num_upstr_bases=num_upstr_bases,
                          num_downstr_bases=num_downstr_bases,
                          generate_mpileup_file=generate_mpileup_file,
                          compress_output=compress_output,
                          bgzip=bgzip,
                          path_to_bgzip=path_to_bgzip,
                          path_to_tabix=path_to_tabix,
                          buffer_line_number = buffer_line_number,
                          min_mapq=min_mapq,
                          min_cov = min_cov,
                          binom_test = binom_test,
                          path_to_samtools = path_to_samtools,
                          sort_mem=sort_mem,
                          path_to_files = path_to_files,
                          min_base_quality = min_base_quality,
                          remove_chr_prefix = remove_chr_prefix,
                          add_snp_info  = add_snp_info,
                          keep_temp_files = keep_temp_files)

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
 
