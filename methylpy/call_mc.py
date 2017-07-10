import sys
import multiprocessing
import subprocess
import scipy.stats as sci
from scipy.stats.mstats import mquantiles
from methylpy.utilities import print_checkpoint,split_mpileup_file,split_fastq_file
import pdb
import shlex
import itertools
import re
import glob
import cStringIO as cStr
import bisect
try:
    from argparse import ArgumentParser
except Exception,e:
    exc_type, exc_obj, exc_tb = exc_info()
    print(exc_type, exc_tb.tb_lineno)
    print e
    exit("methylpy.call_mc requires ArgumentParser from the argparse module")
    
def run_methylation_pipeline(files, libraries, sample,
                             forward_reference, reverse_reference, reference_fasta, unmethylated_control,
                             quality_version="1.8", remove_clonal=True, path_to_picard="",
                             path_to_samtools="", path_to_bowtie="",
                             bowtie_options=["-S","-k 1","-m 1","--chunkmbs 3072",
                                             "--best","--strata","-o 4","-e 80","-l 20","-n 0"],
                             num_procs=1, trim_reads=True,
                             path_to_cutadapt="", adapter_seq = "AGATCGGAAGAGCACACGTCTG", max_adapter_removal=None,
                             overlap_length=None, zero_cap=None, error_rate=None, min_qual_score=10, min_read_len=30,
                             sig_cutoff=.01, min_cov=0, binom_test=True, bh=False,
                             keep_temp_files=False, num_reads=-1, save_space=True,
                             bowtie2=False, sort_mem="500M", path_to_output="", in_mem=False, min_base_quality=1):

    """
    files is a list of all the fastq files you'd like to run through the pipeline. 
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

    unmethylated_control is the name of the chromosome/region that you want to use to estimate the non-conversion rate of your 
        sample, or the non-conversion rate you'd like to use. Consequently, control is either a string, or a decimal
        If control is a string then it should be in the following format: "chrom:start-end". 
        If you'd like to specify an entire chromosome simply use "chrom:"
        
    quality_version is either an integer indicating the base offset for the quality scores or a float indicating
        which version of casava was used to generate the fastq files.

    remove_clonal is a boolean indicating that you want to remove clonal reads (PCR duplicates). If true,
        picard.jar should be available in folder specified in path_to_picard.
    
    path_to_picard is a string of the path to "picard.jar". "picard.jar" is assumed to be 
        in your path if this option isn't used
            
    path_to_samtools is a string indicating the path to the directory containing your 
        installation of samtools. Samtools is assumed to be in your path if this is not
        provided    

    path_to_bowtie is a string indicating the path to the folder in which bowtie resides. Bowtie
        is assumed to be in your path if this option isn't used
            
    bowtie_options is a list of strings indicating options you'd like passed to bowtie 
        (e.g., ["-k 1","-l 2"]
    
    num_procs is an integer indicating how many num_procs you'd like to run this function over
    
    trim_reads is a boolean indicating that you want to have reads trimmed by cutadapt. This option is currently
        not compatible with "in_mem" mode.

    path_to_cutadapt is the path to the cutadapt execuatable. Otherwise this is assumed to be in your
        path.

    adapter_seq is the sequence of an adapter that was ligated to the 3' end. The adapter itself and 
        anything that follows is trimmed.

    max_adapter_removal indicates the maximum number of times to try to remove adapters. Useful when an adapter 
        gets appended multiple times.
        
    overlap_length is the minimum overlap length. If the overlap between the read and the adapter is shorter than 
        LENGTH, the read is not modified. This reduces the no. of bases trimmed purely due to short random adapter matches.

    zero_cap causes negative quality values to be set to zero (workaround to avoid segmentation faults in BWA).
        
    error_rate is the maximum allowed error rate (no. of errors divided by the length of the matching region) 
        (default: 0.1)
    
    min_qual_score allows you to trim low-quality ends from reads before adapter removal. The algorithm is the same as the 
        one used by BWA (Subtract CUTOFF from all qualities; compute partial sums from all indices to the end of the
        sequence; cut sequence at the index at which the sum is minimal).
        
    min_read_len indicates the minimum length a read must be to be kept. Reads that are too short even before adapter removal
        are also discarded. In colorspace, an initial primer is not counted.
        
    sig_cutoff is a float indicating the adjusted p-value cutoff you wish to use for determining whether or not
        a site is methylated

    min_cov is an integer indicating the minimum number of reads for a site to be tested.

    binom_tests indicates that you'd like to use a binomial test, rather than the alternative method outlined here
        https://bitbucket.org/schultzmattd/methylpy/wiki/Methylation%20Calling

    keep_temp_files is a boolean indicating that you'd like to keep the intermediate files generated
        by this function. This can be useful for debugging, but in general should be left False.
        
    num_reads is an integer indicating how many reads you'd like to process at a time. Breaking up the processing this
        way can reduce the size of some of the intermediate files that are created. Use this if you're strapped for space.
        
    save_space indicates whether or not you'd like to perform read collapsing right after mapping or once all the 
        libraries have been mapped. If you wait until after everything has been mapped, the collapsing can be 
        parallelized. Otherwise the collapsing will have to be done serially. The trade-off is that you must keep
        all the mapped files around, rather than deleting them as they are processed, which can take up a 
        considerable amount of space. It's safest to set this to True.
    
    bowtie2 specifies whether to use the bowtie2 aligner instead of bowtie
    
    sort_mem is the parameter to pass to unix sort with -S/--buffer-size command
    
    in_mem: does not write files as part of run_mapping, keeps all input/output in memory
    
    path_to_output is the path to a directory where you would like the output to be stored. The default is the
        same directory as the input fastqs.

    min_base_quality is an integer indicating the minimum PHRED quality score for a base to be included in the
        mpileup file (and subsequently to be considered for methylation calling)        
    """
    #Figure out all the correct quality options based on the offset or CASAVA version given
    try:
        quality_version = float(quality_version)
        if quality_version >= 1.3 and quality_version<1.8:
            quality_base = 64
            bowtie_options.append("--phred64-quals")
        elif quality_version >= 1.8:
            quality_base = 33
            bowtie_options.append("--phred33-quals")
        elif quality_version < 1.3:
            quality_base = 64
            bowtie_options.append("--solexa-quals")
        else:
            sys.exit("Unrecognized quality_version. Either enter the CASAVA version used to generate the fastq files, or enter the base quality offset (either 33 or 64).")
    except:
        if int(quality_version) == 64:
            quality_base = 64
            bowtie_options.append("--phred64-quals")
        elif int(quality_version) == 33:
            quality_base = 33
            bowtie_options.append("--phred33-quals")
        else:
            sys.exit("Unrecognized quality_version. Either enter the CASAVA version used to generate the fastq files, or enter the base quality offset (either 33 or 64).")

    #This code allows the user to supply paths with "*" in them rather than listing out every single
    #file
    #just to avoid any paths missing a slash. It's ok to add an extra slash if
    #there's already one at the end
    if len(path_to_samtools)!=0:
        path_to_samtools +="/"
    if len(path_to_bowtie)!=0:
        path_to_bowtie+="/"
    if len(path_to_output) !=0:
        path_to_output+="/"
        
    if sort_mem:
        if sort_mem.find("-S") == -1:
            sort_mem = " -S " + sort_mem
    else:
        sort_mem = ""
        
    expanded_file_list = []
    expanded_library_list = []
    
    total_reads = 0
    total_unique = 0
    total_clonal = 0
    for path,library in zip(files,libraries):
        glob_list = glob.glob(path)
        for filen in glob_list:
            expanded_file_list.append(filen)
            expanded_library_list.append(library)
            total_reads += int(subprocess.check_output(['wc', '-l', filen]).split(' ')[0]) / 4
    print_checkpoint("There are " + str(total_reads) + " total input reads")
    for current_library in set(libraries):
        library_files = [filen for filen,library in zip(expanded_file_list,expanded_library_list) if library == current_library]
        #deal with actual filename rather than path to file
        if in_mem:
            total_unique += Mrun_mapping(current_library,library_files,sample,
                                         forward_reference,reverse_reference,reference_fasta,
                                         path_to_samtools=path_to_samtools,path_to_bowtie=path_to_bowtie,
                                         bowtie_options=bowtie_options,num_procs=num_procs,
                                         trim_reads=trim_reads,path_to_cutadapt=path_to_cutadapt,
                                         adapter_seq = adapter_seq,max_adapter_removal=max_adapter_removal,
                                         overlap_length=overlap_length,zero_cap=zero_cap,quality_base=quality_base,
                                         error_rate=error_rate,
                                         min_qual_score=min_qual_score,min_read_len=min_read_len,
                                         keep_temp_files=keep_temp_files,
                                         num_reads=num_reads, bowtie2=bowtie2, sort_mem=sort_mem,
                                         path_to_output=path_to_output)
        else:
            total_unique += run_mapping(current_library,library_files,sample,
                                        forward_reference,reverse_reference,reference_fasta,
                                        path_to_samtools=path_to_samtools,path_to_bowtie=path_to_bowtie,
                                        bowtie_options=bowtie_options,num_procs=num_procs,
                                        trim_reads=trim_reads,path_to_cutadapt=path_to_cutadapt,
                                        adapter_seq = adapter_seq,max_adapter_removal=max_adapter_removal,
                                        overlap_length=overlap_length,zero_cap=zero_cap,quality_base=quality_base,
                                        error_rate=error_rate,
                                        min_qual_score=min_qual_score,min_read_len=min_read_len,
                                        keep_temp_files=keep_temp_files,
                                        num_reads=num_reads, bowtie2=bowtie2, sort_mem=sort_mem,
                                        path_to_output=path_to_output,save_space=save_space)            
        if save_space == True:
            pool = multiprocessing.Pool(num_procs)
            processed_library_files = glob.glob(path_to_output+sample+"_"+str(current_library)+"_*_no_multimap_*")
            total_clonal += collapse_clonal_reads(reference_fasta,path_to_samtools,num_procs,sample,current_library,path_to_files=path_to_output,sort_mem=sort_mem,remove_clonal=remove_clonal,path_to_picard=path_to_picard)
    if save_space == True:
        print_checkpoint("There are " + str(total_unique) + " uniquely mapping reads, " + str(float(total_unique) / total_reads*100) + " percent remaining")
    else:
        print_checkpoint("Finding multimappers")
        total_unique = 0
        merge_results= []
        sorted_files = []
        for library in set(libraries):
            pool = multiprocessing.Pool(num_procs)
            library_sorted_files = glob.glob(sample+"_"+str(library)+"_[0-9]*_sorted_[0-9]*")
            sorted_files.extend(library_sorted_files)
            sort_results = []
            for filename in library_sorted_files:
                sort_results.append(pool.apply_async(subprocess.check_call,(shlex.split("env LC_COLLATE=C sort" + sort_mem + " -t '\t' -k 1 -o "+filename+" "+filename),)))
            #ensure the sorting has finished
            for result in sort_results:
                result.get()
            #Need to add the _1 for a bogus chunk ID
            merge_results.append(pool.apply_async(merge_sorted_multimap,(sorted_files,sample+"_"+str(library)+"_1")))
        pool.close()
        pool.join()
        subprocess.check_call(shlex.split("rm "+" ".join(sorted_files)))
        for result in merge_results:
            total_unique+=result.get()
        print_checkpoint("There are " + str(total_unique) + " uniquely mapping reads, " + str(float(total_unique) / total_reads*100) + " percent remaining")
        
        
        pool = multiprocessing.Pool(num_procs)
        #I avoid sorting by chromosome or strand because it's not strictly necessary and this speeds up 
        #the sorting. I had to add some special logic in the processing though
        processed_library_files = glob.glob(path_to_output+sample+"_"+str(current_library)+"_*_no_multimap_*")
        
        clonal_results = []
        pool = multiprocessing.Pool(num_procs)

        for library in set(libraries):
            clonal_results.append(pool.apply_async(collapse_clonal_reads,(reference_fasta,path_to_samtools,num_procs,sample,library), {"sort_mem":sort_mem,"path_to_files":path_to_output,"remove_clonal":remove_clonal,"path_to_picard":path_to_picard}))
        pool.close()
        pool.join()

        if remove_clonal==True:        
            for result in clonal_results:
                total_clonal += result.get()

    if remove_clonal==True:
        print_checkpoint("There are " + str(total_clonal) + " non-clonal reads, " + str(float(total_clonal) / total_reads*100) + " percent remaining")
        library_files = [path_to_output+sample+"_processed_reads_"+str(library)+"_no_clonal.bam" for library in set(libraries)]
        if len(library_files) > 1:
            merge_bam_files(library_files,path_to_output+sample+"_processed_reads_no_clonal.bam",path_to_samtools)
            subprocess.check_call(shlex.split("rm "+" ".join(library_files)))
        else:
            subprocess.check_call(shlex.split("mv "+library_files[0]+" "+path_to_output+sample+"_processed_reads_no_clonal.bam"))
    else:
        library_files = [path_to_output+sample+"_processed_reads_"+str(library)+".bam" for library in set(libraries)]
        if len(library_files) > 1:
            merge_bam_files(library_files,path_to_output+sample+"_processed_reads.bam",path_to_samtools)
            subprocess.check_call(shlex.split("rm "+" ".join(library_files)))
        else:
            subprocess.check_call(shlex.split("mv "+library_files[0]+" "+path_to_output+sample+"_processed_reads.bam"))

    print_checkpoint("Begin calling mCs")
    if remove_clonal==True:
        call_methylated_sites(sample+"_processed_reads_no_clonal.bam",sample,reference_fasta,unmethylated_control,quality_version,sig_cutoff=sig_cutoff,num_procs=num_procs,min_cov=min_cov,binom_test=binom_test,bh=bh,sort_mem=sort_mem,path_to_files=path_to_output,path_to_samtools=path_to_samtools,min_base_quality=min_base_quality)
    else:
        call_methylated_sites(sample+"_processed_reads.bam",sample,reference_fasta,unmethylated_control,quality_version,sig_cutoff=sig_cutoff,num_procs=num_procs,min_cov=min_cov,binom_test=binom_test,bh=bh,sort_mem=sort_mem,path_to_files=path_to_output,path_to_samtools=path_to_samtools,min_base_quality=min_base_quality)
    print_checkpoint("Done")
    
def run_mapping(current_library,library_files,sample,forward_reference,reverse_reference,reference_fasta,
                path_to_samtools="",path_to_bowtie="",
                bowtie_options=["-S","-k 1","-m 1","--chunkmbs 3072","--best","--strata","-o 4","-e 80","-l 20","-n 0"],
                num_procs=1,trim_reads=True,path_to_cutadapt="",adapter_seq = "AGATCGGAAGAGCACACGTCTG",
                max_adapter_removal=None,overlap_length=None,zero_cap=None,quality_base=None,error_rate=None,
                min_qual_score=10,min_read_len=30,keep_temp_files=False,num_reads=-1, bowtie2=False, sort_mem="500M",
                path_to_output="",save_space=True):
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
    
    path_to_bowtie is a string indicating the path to the folder in which bowtie resides. Bowtie
        is assumed to be in your path if this option isn't used
            
    bowtie_options is a list of strings indicating options you'd like passed to bowtie 
        (e.g., ["-k 1","-l 2"]
    
    num_procs is an integer indicating how many num_procs you'd like to run this function over
    
    trim_reads is a boolean indicating that you want to have reads trimmed by cutadapt

    path_to_cutadapt is the path to the cutadapt execuatable. Otherwise this is assumed to be in your
        path.

    adapter_seq is the sequence of an adapter that was ligated to the 3' end. The adapter itself and 
        anything that follows is trimmed.

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
        
    num_reads is an integer indicating how many reads you'd like to process at a time. Breaking up the processing this
        way can reduce the size of some of the intermediate files that are created. Use this if you're strapped for space.
    
    bowtie2 specifies whether to use the bowtie2 aligner instead of bowtie
    
    sort_mem is the parameter to pass to unix sort with -S/--buffer-size command
    """
    total_unique = 0
    chunk = 1
    file_name = sample+"_"+str(current_library)+"_"+str(chunk)
    file_path = path_to_output+file_name
    save_point = (False,0)
    while save_point != False:
        print_checkpoint("Begin chunk "+str(chunk)+" for library "+str(current_library))
        print_checkpoint("Begin splitting reads for "+file_name)
        save_point = split_fastq_file(num_procs,library_files,file_path+"_split_",num_reads=num_reads,current_file=save_point[0],current_pointer=save_point[1])
        if trim_reads:
            print_checkpoint("Begin trimming reads for "+file_name)
            quality_trim([file_path+"_split_"+str(i) for i in xrange(0,num_procs)],output=[file_path+"_split_trimmed_"+str(i) for i in xrange(0,num_procs)],adapter_seq=adapter_seq,error_rate=error_rate,quality_base = quality_base, min_qual_score=min_qual_score,min_read_len=min_read_len,format="fastq",num_procs=num_procs,max_adapter_removal=max_adapter_removal,overlap_length=overlap_length,zero_cap=zero_cap,path_to_cutadapt=path_to_cutadapt)
            subprocess.check_call(shlex.split("rm "+" ".join([file_path+"_split_"+str(i) for i in xrange(0,num_procs)])))
            print_checkpoint("Begin converting reads for "+file_name)
            pool = multiprocessing.Pool(num_procs)
            for inputf,output in zip([file_path+"_split_trimmed_"+str(i) for i in xrange(0,num_procs)],
                                     [file_path+"_split_trimmed_converted_"+str(i) for i in xrange(0,num_procs)]):
                pool.apply_async(convert_reads,(inputf,output))
            pool.close()
            pool.join()
            subprocess.check_call(shlex.split("rm "+
                                              " ".join([file_path+"_split_trimmed_"+str(i) for i in xrange(0,num_procs)])))
            print_checkpoint("Begin Running Bowtie for "+file_name)
            total_unique += run_bowtie([file_path+"_split_trimmed_converted_"+str(i) for i in xrange(0,num_procs)],
                                       forward_reference,reverse_reference,file_path,options=bowtie_options,
                                       path_to_bowtie=path_to_bowtie,num_procs=num_procs,keep_temp_files=keep_temp_files,
                                       bowtie2=bowtie2, sort_mem=sort_mem,save_space=save_space)
            if keep_temp_files==False:
                subprocess.check_call(shlex.split("rm "+" ".join(
                    [file_path+"_split_trimmed_converted_"+str(i) for i in xrange(0,num_procs)])))

        else:
            print_checkpoint("No trimming on reads")
            print_checkpoint("Begin converting reads for "+file_name)
            pool = multiprocessing.Pool(num_procs)
            for inputf,output in zip([file_path+"_split_"+str(i) for i in xrange(0,num_procs)],
                                     [file_path+"_split_converted_"+str(i) for i in xrange(0,num_procs)]):
                pool.apply_async(convert_reads,(inputf,output))
            pool.close()
            pool.join()
            print_checkpoint("Begin Running Bowtie for "+file_name)
            total_unique += run_bowtie([file_path+"_split_converted_"+str(i) for i in xrange(0,num_procs)],
                                       forward_reference,reverse_reference,file_path,options=bowtie_options,
                                       path_to_bowtie=path_to_bowtie,num_procs=num_procs,keep_temp_files=keep_temp_files,
                                       bowtie2=bowtie2, sort_mem=sort_mem,save_space=save_space)                
        chunk+=1

    return total_unique

def Mrun_mapping(current_library,library_files,sample,
                 forward_reference,reverse_reference,reference_fasta,path_to_samtools="",path_to_bowtie="",
                 bowtie_options=["-S","-k 1","-m 1","--chunkmbs 3072","--best","--strata","-o 4","-e 80","-l 20","-n 0"],
                 num_procs=1,trim_reads=True,path_to_cutadapt="",adapter_seq = "AGATCGGAAGAGCACACGTCTG",
                 max_adapter_removal=None,overlap_length=None,zero_cap=None,quality_base=None,error_rate=None,
                 min_qual_score=10,min_read_len=30,keep_temp_files=False,num_reads=-1,bowtie2=False,sort_mem="500M",
                 path_to_output=""):
    total_unique = 0
    file_name = sample+"_"+str(current_library)
    file_path = path_to_output+file_name
    save_point = (False,0)
    chunk = 1
    while save_point != False:
        print_checkpoint("Begin chunk "+str(chunk)+" for library "+str(current_library))
        print_checkpoint("Begin splitting reads for "+file_name)
        save_point, file_handles = Msplit_fastq_file(num_procs,library_files,num_reads=num_reads,
                                                     current_file=save_point[0],current_pointer=save_point[1])
        print_checkpoint("Begin trimming reads for "+file_name)
        trimmed = Mquality_trim(file_handles,adapter_seq=adapter_seq,error_rate=error_rate,quality_base = quality_base, min_qual_score=min_qual_score,min_read_len=min_read_len,format="fastq",num_procs=num_procs,max_adapter_removal=max_adapter_removal,overlap_length=overlap_length,zero_cap=zero_cap,path_to_cutadapt=path_to_cutadapt)
        print_checkpoint("Begin converting reads for "+file_name)
        converted = Mconvert_reads(trimmed, num_procs)
        print_checkpoint("Begin Running Bowtie for "+file_name)
        total_unique += Mrun_bowtie(converted,forward_reference,reverse_reference,file_path,options=bowtie_options,path_to_bowtie=path_to_bowtie,num_procs=num_procs,keep_temp_files=keep_temp_files, bowtie2=bowtie2, sort_mem=sort_mem)
        chunk+=1

    return total_unique

def Msplit_fastq_file(num_chunks, input_files, num_reads = -1,current_file=False,current_pointer=0):
    count = 0
    if not isinstance(input_files, list):
        if isinstance(input_files, basestring):
            input_files = [input_files]
        else:
            sys.exit("input_files must be a list of strings")
    if current_file != False:
        input_files = input_files[input_files.index(current_file):]
    file_handles = {}
    for index in xrange(0,num_chunks):
        file_handles[index]=cStr.StringIO()    
    cycle = itertools.cycle(range(0,num_chunks))
    for inputf in input_files:
        try:
            f = gzip.open(inputf,'r')
            f.readline()
        except:
            try:
                f = bz2.BZ2File(inputf,'r')
                f.readline()
            except:
                f = open(inputf,'r')
        finally:
            f.seek(current_pointer)

        line = True
        while line:
            current_file = cycle.next()
            for index in xrange(0,4):
                line = f.readline()
                file_handles[current_file].write(line)
            count += 1
            if num_reads != -1 and count >= num_reads:                
                save_point = ((inputf,f.tell()), file_handles)
                f.close()
                return save_point
        f.close()
        current_pointer=0

    return (False, file_handles)
 
def Mquality_trim(file_handles, output = None, quality_base = None, min_qual_score = None, min_read_len = None, 
                 adapter_seq = "AGATCGGAAGAGCACACGTCTG", num_procs = 1, format = None, error_rate = None, 
                 max_adapter_removal = None, overlap_length = None, zero_cap = False, path_to_cutadapt = ""):
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
            
    base_cmd = path_to_cutadapt
    options = ""
    if zero_cap:
        zero = "-z "
    else:
        zero = ""  
    
    if format:
        options += " -f " + format
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
    cmd = shlex.split(base_cmd + options + "-")
    pool = multiprocessing.Pool(num_procs)
    results = []
    #adapter trimming
    for key in file_handles:
        results.append(pool.apply_async(open_cmd, (cmd, file_handles[key].getvalue())))
        file_handles[key].close()
    pool.close()
    pool.join()
    trimmed = {}
    for key in file_handles:
        trimmed[key] = results[int(key)].get()
    
    return trimmed
        
def open_cmd(cmd, fstring):
    """
    Helper function for Mquality_trim
    """
    proc = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    out, err = proc.communicate(input=fstring)
    return out

def Mconvert_reads(file_handles, num_procs):
    pool = multiprocessing.Pool(num_procs)
    results = []
    #adapter trimming
    for key in file_handles:
        results.append(pool.apply_async(Mconvert_reads_helper, (file_handles[key],)))
    pool.close()
    pool.join()
    converted = {}
    for key in file_handles:
        converted[key] = results[int(key)].get()
    
    return converted

def Mconvert_reads_helper(readstr):
    g = cStr.StringIO()
    f = cStr.StringIO(readstr)
    header = f.readline().rstrip()
    header = header.replace(" ","!")
    seq = f.readline()
    header2 = f.readline()
    qual = f.readline()
    encoding = encode_c_positions(seq)
    
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
        encoding = encode_c_positions(seq)
    f.close()
    
    return g.getvalue()

def Mrun_bowtie(converted,forward_reference,reverse_reference,prefix,options="",path_to_bowtie="",
               num_procs=1,keep_temp_files=False, bowtie2=False, sort_mem="500M"):
    inreads = []
    for key in converted:
        inreads.append(converted[key])
    inreads = "".join(inreads)
    
    if sort_mem:
        if sort_mem.find("-S") == -1:
            sort_mem = " -S " + sort_mem
    else:
        sort_mem = ""
    if " ".join(options).find(" -p ") == -1:
        options.append("-p "+str(num_procs))
    if bowtie2:
        args = [path_to_bowtie+"bowtie2"]
        args.extend(options)
        args.append("--norc")
        args.append("-x "+forward_reference)
        args.append("-U "+",".join(files))
        args.append("-S "+prefix+"_forward_strand_hits.sam")

    else:
        args = [path_to_bowtie+"bowtie"]
        args.extend(options)
        args.append("--norc")
        args.append(forward_reference)
        args.append(",".join(files))
        args.append(prefix+"_forward_strand_hits.sam")
    proc = subprocess.Popen(shlex.split(" ".join(args)), stdin=subprocess.PIPE)
    proc.communicate(input=inreads)

    print_checkpoint("Processing forward strand hits")
    find_multi_mappers(prefix+"_forward_strand_hits.sam",prefix,num_procs=num_procs,keep_temp_files=keep_temp_files)
    if bowtie2:
        args = [path_to_bowtie+"bowtie2"]
        args.extend(options)
        args.append("--nofw")
        args.append("-x "+reverse_reference)
        args.append("-U "+",".join(files))
        args.append("-S "+prefix+"_reverse_strand_hits.sam")
    else:
        args = [path_to_bowtie+"bowtie"]
        args.extend(options)
        args.append("--nofw")
        args.append(reverse_reference)
        args.append(",".join(files))
        args.append(prefix+"_reverse_strand_hits.sam")
    proc = subprocess.Popen(shlex.split(" ".join(args)), stdin=subprocess.PIPE)
    proc.communicate(input=inreads)
        
    print_checkpoint("Processing reverse strand hits")
    sam_header = find_multi_mappers(prefix+"_reverse_strand_hits.sam",prefix,num_procs=num_procs,append=True,keep_temp_files=keep_temp_files)
    pool = multiprocessing.Pool(num_procs)
    for file_num in xrange(0,num_procs):
        pool.apply_async(subprocess.check_call,(shlex.split("env LC_COLLATE=C sort" + sort_mem + " -t '\t' -k 1 -o "+prefix+"_sorted_"+str(file_num)+" "+prefix+"_sorted_"+str(file_num)),))
    pool.close()
    pool.join()
    print_checkpoint("Finding multimappers")
    total_unique = merge_sorted_multimap([prefix+"_sorted_"+str(file_num) for file_num in xrange(0,num_procs)],prefix)
    subprocess.check_call(shlex.split("rm "+" ".join([prefix+"_sorted_"+str(file_num) for file_num in xrange(0,num_procs)])))
    #return "".join(sorted(set(sam_header)))
    return total_unique

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
        f.write("@RG\tID:" + filen[:filen.rindex(".bam")] + "\tLB:" + filen[filen.rindex("processed_reads_")+16:filen.rindex("_no_clonal.bam")] + "\tSM:NA" + "\n")
    f.close()
    subprocess.check_call(shlex.split(path_to_samtools+"samtools merge -r -h header.sam "+ output +" "+" ".join(input_files)))
    subprocess.check_call(["rm", "header.sam"])

def merge_unknown_bam(old_file, new_file, output, path_to_samtools=""):
    """
    This function merges two bam files. old_file can consist of many libraries, but new_file must only contain one library.
    If there are duplicates of reads from the same library in new_file and old_file, only one copy will be written after
    merging.
    """
    #find out if these come from the same library or not
    
    #find all libary names present in the old file
    old_header = subprocess.check_output(shlex.split(path_to_samtools+"samtools view -H "+old_file))
    rg1_index = [] #indices of '@RG' substrings in header for old_file
    lib1 = [] #list of (library name, id name) found in old_file
    start = 0
    #find all instances of "@RG" in the header of old_file
    while True: 
        start = old_header.find("@RG", start)
        if start==-1:
            break
        rg1_index.append(start)
        start += 3
    
    #get the ID and LB from each @RG value in the header
    for index in rg1_index:
        rg1_start = old_header.find("LB:", index) + 3
        rg1_end = old_header.find("\n", rg1_start) # (!) assumption that LB: is the last field on the @RG line    
        id1_start = old_header.find("ID:", index) + 3
        id1_end = old_header.find("\t", id1_start)
        lib1.append((old_header[rg1_start:rg1_end], old_header[id1_start:id1_end]))
    
    #find library name for this file, (!) assumption that there is only one library in new_file
    new_header = subprocess.check_output(shlex.split(path_to_samtools+"samtools view -H "+new_file)) 
    rg2_index = new_header.find("@RG")

    #If no @RG can be found in the file, try to use the library found in the filename
    id_filename = False #flag indicating if the ID is the current new_file filename to prevent overwriting
    if rg2_index == -1:
        try:
            lib2 = new_file[new_file.rindex("processed_reads_")+16:new_file.rindex("_no_clonal.bam")]
            id_filename = True
        except:
            sys.exit("Could not find library information for: " + new_file)
    else: #get the LB information from the @RG tag
        rg2_start = new_header.find("LB:", rg2_index) + 3
        rg2_end = new_header.find("\n", rg2_start)
        rg2_tab = new_header.find("\t", rg2_start) #in case the "LB:" is not the last field and there is \t instead of \n
        if rg2_end == -1 or (rg2_end!=-1 and rg2_tab!=-1 and rg2_tab < rg2_end):
            rg2_end = rg2_tab
        lib2 = new_header[rg2_start:rg2_end]
    
    #print "Library 2: " + lib2
    #print "Library 1: " + ", ".join([x[0] for x in lib1]) + "\n"
    
    if lib2 not in [x[0] for x in lib1]:
        print_checkpoint("New library not found in existing file, merging")
        #merge the two files together without worrying about duplicates
        f=open("header.sam",'w')
        subprocess.check_call(shlex.split(path_to_samtools+"samtools view -H "+old_file),stdout=f)
        f.write("@RG\tID:" + new_file[:new_file.rindex(".bam")] + "\tLB:" + lib2 + "\tSM:NA" + "\n")
        f.close()
        subprocess.check_call(shlex.split(path_to_samtools+"samtools merge -r -h header.sam "+ output +" "+new_file+" "+old_file))
        subprocess.check_call(["rm", "header.sam"])
    else:
        print_checkpoint("Finding duplicates and merging BAM files")
        try:
            g = open(old_file + ".bai", 'r') #check if old_file has a samtools index (needed to search by region)
            g.close()
        except:
            subprocess.check_call(shlex.split(path_to_samtools + "samtools index " + old_file))
        
        #Get the id of the matching library in old_file
        lib2_id = ""
        for lb, id in lib1:
            if lb==lib2:
                lib2_id = id
                break
        
        #check if the old file has the same filename as the new file that we are trying to create
        if old_file == lib2_id + ".bam":
            sys.exit("Error when merging bam files: the old bam file must not be the same filename as the one that it is merging with\n")
        
        #create the header from old_file because the @RG is the same in the new file
        g = open("header.sam", 'w')
        subprocess.check_call(shlex.split(path_to_samtools+"samtools view -H "+old_file),stdout=g)
        g.close()
        
        lines = subprocess.check_output(shlex.split(path_to_samtools+"samtools view "+new_file))
        
        if id_filename: #prevent overwriting the input file if the filename is going to be the same as the new merge file
            subprocess.check_call(["mv", new_file, new_file+"tmp"])
            
        f = open(lib2_id + ".sam",'w')
        #header is required in sam file, so put it in the new file to merge
        if id_filename:
            subprocess.check_call(shlex.split(path_to_samtools+"samtools view -H "+new_file+"tmp"),stdout=f)
        else:
            subprocess.check_call(shlex.split(path_to_samtools+"samtools view -H "+new_file),stdout=f)
        
        lines = lines.split('\n')
        written = False #flag used to indicate if any new reads were found (if not, the merge command will fail if attempted)
        for line in lines:
            if not line:
                continue
            line = line.split('\t')
            chr = line[2]
            start = line[3]
            strand = line[1]
            #check to see if this read exists in the old file
            results = subprocess.check_output(shlex.split(path_to_samtools+"samtools view "+old_file+" "+chr+":"+start+"-"+start))
            results = results.split("\n")
            valid = True #if any duplicate read is found, valid becomes false and the line is not written
            for read in results:
                if not read:
                    continue
                read = read.split("\t")
                if read[2]==chr and read[3]==start and read[1]==strand and read[-1]=="RG:Z:"+lib2_id:
                   valid = False #this read already exists in the old_file
                   break 
            if valid:
                f.write('\t'.join(line) + '\n')
                written = True
        f.close()
        
        if written:       
            with open(lib2_id + ".bam", 'w') as outfile:
                subprocess.check_call(shlex.split(path_to_samtools+"samtools view -bS " + lib2_id + ".sam"), stdout=outfile)
            subprocess.check_call(shlex.split(path_to_samtools+"samtools merge -r -h header.sam "+ output +" "+old_file+" "+lib2_id + ".bam"))
            subprocess.check_call(["rm", "header.sam", lib2_id + ".bam", lib2_id + ".sam"])
        else:
            print_checkpoint("No new reads were found")
            subprocess.check_call(["rm", "header.sam", lib2_id + ".sam"])
        
        if id_filename:
            subprocess.check_call(["mv", new_file+"tmp", new_file])
    
    print_checkpoint("Finished merging BAM files")

def collapse_clonal_reads(reference_fasta,path_to_samtools,num_procs,sample,
                          current_library, sort_mem="500M",path_to_files="",
                          remove_clonal=True,path_to_picard=""):
    """
    This function is a wrapper for collapsing clonal reads
    
    reference_fasta is a string indicating the path to a fasta file containing the sequences
        you used for mapping
        input is the path to a bam file that contains mapped bisulfite sequencing reads
    
    path_to_samtools is a string indicating the path to the directory containing your 
        installation of samtools. Samtools is assumed to be in your path if this is not
        provided    

    num_procs is an integer indicating how many num_procs you'd like to run this function over

    sample is a string indicating the name of the sample you're processing. It will be included
        in the output files.
    
    current_library is the library from which clonal reads should be removed. The expected file
        name is <sample>_processed_reads_<current_library>_no_clonal
        
    sort_mem is the parameter to pass to unix sort with -S/--buffer-size command
    
    path_to_files is a string indicating the path for the output files and input sam for clonal
        collapsing.

    remove_clonal is a boolean indicating that you want to remove clonal reads (PCR duplicates). If true,
        picard.jar should be available in folder specified in path_to_picard.
    
    path_to_picard is a string of the path to "picard.jar". "picard.jar" is assumed to be 
        in your path if this option isn't used
    """
    if remove_clonal == True:
        print_checkpoint("Collapsing clonal reads for library "+str(current_library))
    else:
        print_checkpoint("Processing uniquely mapped reads for library "+str(current_library))

    #Include the * between the library and _no_multimap to account for different chunks
    #Reads in each of these files are sorted by the position
    processed_library_files = glob.glob(path_to_files+sample+"_"+str(current_library)+"_*_no_multimap_*")
    merge_no_multimap(processed_library_files,
                      path_to_files+sample+"_processed_reads_"+str(current_library)+"",
                      reference_fasta,path_to_samtools=path_to_samtools,num_procs=num_procs,sort_mem=sort_mem)
    
    print_checkpoint("Converting to BAM")
    f = open(path_to_files+sample+"_processed_reads_"+str(current_library)+".bam",'w')
    subprocess.check_call(shlex.split(path_to_samtools+"samtools view -S -b -h "+path_to_files+sample+"_processed_reads_"+str(current_library)),stdout=f)
    subprocess.check_call(shlex.split(" ".join(["rm",
                                                path_to_files+sample+"_processed_reads_"+str(current_library)])))
    subprocess.check_call(shlex.split(path_to_samtools+"samtools sort "+path_to_files+sample+"_processed_reads_"+str(current_library)+".bam "+path_to_files+sample+"_processed_reads_"+str(current_library)))
    
    f.close()

    total_clonal = 0
    #Remove clonal reads
    if remove_clonal == True:
        print_checkpoint("Remove clonal reads for library "+str(current_library))
        subprocess.check_call(shlex.split(" ".join(["java","-Xmx2g","-jar",
                                                    path_to_picard+"/picard.jar MarkDuplicates",
                                                    "INPUT="+path_to_files+sample+"_processed_reads_"+str(current_library)+".bam",
                                                    "OUTPUT="+path_to_files+sample+"_processed_reads_"+str(current_library)+"_no_clonal.bam",
                                                    "ASSUME_SORTED=true",
                                                    "REMOVE_DUPLICATES=true",
                                                    "METRICS_FILE=/dev/null",
                                                    "VALIDATION_STRINGENCY=LENIENT"])))
        subprocess.check_call(shlex.split(" ".join(["rm",
                                                    path_to_files+sample+"_processed_reads_"+str(current_library)+".bam"])))
        total_clonal = subprocess.check_output(
            shlex.split(" ".join([path_to_samtools,
                                  "samtools view -c",
                                  path_to_files+sample+"_processed_reads_"+str(current_library)+"_no_clonal.bam"])
            ))
    
    return int(total_clonal)

def build_ref(input_files, output, buffsize=100, offrate=False,parallel=False,bowtie2=False):
    """
    Creates 2 reference files: one with all C's converted to T's, and one with all G's converted to A's
    
    input_files is a list of files to build a reference from
    
    output is the prefix of the two output reference files that will be created
    
    buffsize is the number of bytes that will be read in from the reference at once
    
    offrate refers to the Bowtie parameter, reference the bowtie manual to see more detail
    """
    if not isinstance(input_files, list):
        if isinstance(input_files, basestring):
            input_files = [input_files]
        else:
            sys.exit("input_files must be a list of strings")
    print input_files
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
    if bowtie2:
        base_cmd = "bowtie2-build -f "
    else:
        base_cmd = "bowtie-build -f "
    if offrate:
        base_cmd += "-o " + str(offrate) + " "
    
    if parallel == True:
        pool = multiprocessing.Pool(2)
        pool.apply_async(subprocess.check_call,(shlex.split(base_cmd + output + "_f.fasta " + output +"_f"),))
        pool.apply_async(subprocess.check_call,(shlex.split(base_cmd + output + "_r.fasta " + output+ "_r"),))
        pool.close()
        pool.join()
    else:
        subprocess.check_call(shlex.split(base_cmd + output + "_f.fasta " + output +"_f"))
        subprocess.check_call(shlex.split(base_cmd + output + "_r.fasta " + output+ "_r"))
    
def convert_reads(inputf,output):
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
        encoding = encode_c_positions(seq)
    f.close()
    g.close()
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

def encode_converted_positions(seq,is_R2=False):
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
    if is_R2==False:        
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

def decode_converted_positions(seq,indexes,strand,is_R2=False):
    """
    This function takes the encodings generated by encode_c_position and replaces the appropriate
    positions with C nucleotides.
    
    seq is a string of nucleotides to have Cs or Gs replaced.
    
    indexes is a string of characters indicating the offsets for the positions of the Cs or Gs.
    
    strand is the DNA strand (+ or -) that seq mapped to. This is important because
        sequences in sam files are always represented on the forward strand
    """

    prev_index = 0
    new_seq=""
    index = 0
    if is_R2 == False:
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

def run_bowtie(files,forward_reference,reverse_reference,prefix,options="",path_to_bowtie="",
               num_procs=1,keep_temp_files=False, bowtie2=False, sort_mem="500M",save_space=True):
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
    
    path_to_bowtie is a string indicating the path to the folder in which bowtie resides. Bowtie
        is assumed to be in your path if this option isn't used
    
    num_procs is an integer indicating the number of processors you'd like used for removing multi
        mapping reads and for bowtie mapping
    
    keep_temp_files is a boolean indicating that you'd like to keep the intermediate files generated
        by this function. This can be useful for debugging, but in general should be left False.
    
    bowtie2 specifies whether to use the bowtie2 aligner instead of bowtie
    
    sort_mem is the parameter to pass to unix sort with -S/--buffer-size command
    """
    if sort_mem:
        if sort_mem.find("-S") == -1:
            sort_mem = " -S " + sort_mem
    else:
        sort_mem = ""
    if " ".join(options).find(" -p ") == -1:
        options.append("-p "+str(num_procs))
    if bowtie2:
        args = [path_to_bowtie+"bowtie2"]
        args.extend(options)
        args.append("--norc")
        args.append("-x "+forward_reference)
        args.append("-U "+",".join(files))
        args.append("-S "+prefix+"_forward_strand_hits.sam")

    else:
        args = [path_to_bowtie+"bowtie"]
        args.extend(options)
        args.append("--norc")
        args.append(forward_reference)
        args.append(",".join(files))
        args.append(prefix+"_forward_strand_hits.sam")
    subprocess.check_call(shlex.split(" ".join(args)))
    print_checkpoint("Processing forward strand hits")
    find_multi_mappers(prefix+"_forward_strand_hits.sam",prefix,num_procs=num_procs,keep_temp_files=keep_temp_files)
    if bowtie2:
        args = [path_to_bowtie+"bowtie2"]
        args.extend(options)
        args.append("--nofw")
        args.append("-x "+reverse_reference)
        args.append("-U "+",".join(files))
        args.append("-S "+prefix+"_reverse_strand_hits.sam")
    else:
        args = [path_to_bowtie+"bowtie"]
        args.extend(options)
        args.append("--nofw")
        args.append(reverse_reference)
        args.append(",".join(files))
        args.append(prefix+"_reverse_strand_hits.sam")
    subprocess.check_call(shlex.split(" ".join(args)))
    
    print_checkpoint("Processing reverse strand hits")
    sam_header = find_multi_mappers(prefix+"_reverse_strand_hits.sam",prefix,num_procs=num_procs,append=True,keep_temp_files=keep_temp_files)
    if save_space == True:
        pool = multiprocessing.Pool(num_procs)
        for file_num in xrange(0,num_procs):
            pool.apply_async(subprocess.check_call,(shlex.split("env LC_COLLATE=C sort" + sort_mem + " -t '\t' -k 1 -o "+prefix+"_sorted_"+str(file_num)+" "+prefix+"_sorted_"+str(file_num)),))
        pool.close()
        pool.join()
        print_checkpoint("Finding multimappers")

        total_unique = merge_sorted_multimap([prefix+"_sorted_"+str(file_num) for file_num in xrange(0,num_procs)],prefix)
        subprocess.check_call(shlex.split("rm "+" ".join([prefix+"_sorted_"+str(file_num) for file_num in xrange(0,num_procs)])))
        return total_unique
    else:
        return 0
 
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
def merge_sorted_multimap(files, output):
    """
    This function takes the files from find_multi_mappers and outputs the uniquely mapping reads.
    
    files is a list of filenames containing the output of find_multi_mappers
    
    output is a prefix you'd like prepended to the bam file containing the uniquely mapping reads
        This file will be named as <output>+"_no_multimap_"+<index_num>
    """
    
    lines = {}
    fields = {}
    output_handles = {}
    file_handles = {}
    
    total_unique = 0
    count= 0
    cycle = itertools.cycle(range(0,len(files)))
    
    for index,filen in enumerate(files):
        output_handles[index] = open(output+"_no_multimap_"+str(index),'a')               
        file_handles[filen]=open(filen,'r')
        lines[filen]=file_handles[filen].readline()
        fields[filen] = lines[filen].split("\t")[0]
    while True:
        all_fields = [field for field in fields.values() if field != ""]
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
            index = cycle.next()
            output_handles[index].write(current_line)
            output_handles[index].flush()
            total_unique += 1
            
    for index,filen in enumerate(files):
        output_handles[index].close()                
        file_handles[filen].close()

    return total_unique
        
def merge_no_multimap(files,output,reference_fasta,num_procs=1,path_to_samtools="",sort_mem="500M"):
    """
    This function takes a list of *_no_multimap files and merge them into a sam file
    
    files is a list of files that you wish to have clonal reads removed from and reads in each of these files are sorted by position
    
    output is a string indicating the prefix you'd like prepended to the file containing the
        non-clonal output
    
    reference_fasta is a string indicating the path to a fasta file containing the sequences
        you used for mapping
    
    num_procs is an integer indicating the number of files you'd like to split the input into
    
    path_to_samtools is a string indicating the path to the directory containing your 
        installation of samtools. Samtools is assumed to be in your path if this is not
        provided
        
    sort_mem is the parameter to pass to unix sort with -S/--buffer-size command
    """

    g = open(output,'w')
    try:
        f = open(reference_fasta+".fai",'r')
    except:
        print "Reference fasta not indexed. Indexing."
        try:
            subprocess.check_call(shlex.split(path_to_samtools+"samtools faidx "+reference_fasta))
            f = open(reference_fasta+".fai",'r')
        except:
            sys.exit("Reference fasta wasn't indexed, and couldn't be indexed. Please try indexing it manually and running methylpy again.")
    #Write header for sam file
    g.write("@HD\tVN:1.0\tSO:unsorted\n")
    for line in f:
        fields = line.split("\t")
        g.write("@SQ\tSN:"+fields[0]+"\tLN:"+fields[1]+"\n")
    f.close()
    #Intialization
    for filen in files:
        file_handle=open(filen,'r')
        for line in file_handle:
            g.write(line)
        file_handle.close()
    g.close()
    subprocess.check_call(shlex.split("rm "+" ".join(files)))

def run_methylation_pipeline_pe(read1_files,read2_files,libraries,sample,
                                forward_reference,reverse_reference,reference_fasta,unmethylated_control,
                                quality_version="1.8",remove_clonal=True,path_to_picard="",
                                path_to_samtools="",path_to_bowtie="",
                                bowtie_options=[],
                                num_procs=1,trim_reads=True,path_to_cutadapt="",
                                adapter_seq_R1 = "AGATCGGAAGAGCACACGTCTGAAC",
                                adapter_seq_R2 = "AGATCGGAAGAGCGTCGTGTAGGGA",
                                max_adapter_removal=None,
                                overlap_length=None,zero_cap=None,error_rate=None,min_qual_score=10,
                                min_read_len=30,sig_cutoff=.001,min_cov=0,binom_test=True,bh=False,
                                keep_temp_files=False,num_reads=-1,save_space=True,
                                bowtie2=True,sort_mem="500M",path_to_output="",in_mem=False,min_base_quality=1):
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

    unmethylated_control is the name of the chromosome/region that you want to use to estimate the non-conversion rate of your 
        sample, or the non-conversion rate you'd like to use. Consequently, control is either a string, or a decimal
        If control is a string then it should be in the following format: "chrom:start-end". 
        If you'd like to specify an entire chromosome simply use "chrom:"
        
    quality_version is either an integer indicating the base offset for the quality scores or a float indicating
        which version of casava was used to generate the fastq files.

    remove_clonal is a boolean indicating that you want to remove clonal reads (PCR duplicates). If true,
        executable picard should be available in folder specified in path_to_picard.
    
    path_to_picard is a string of the path to "picard.jar". "picard.jar" is assumed to be 
        in your path if this option isn't used
            
    path_to_samtools is a string indicating the path to the directory containing your 
        installation of samtools. Samtools is assumed to be in your path if this is not
        provided    
    
    path_to_bowtie is a string indicating the path to the folder in which bowtie resides. Bowtie
        is assumed to be in your path if this option isn't used
            
    bowtie_options is a list of strings indicating options you'd like passed to bowtie2 (or bowtie)
    (default for bowtie2: "-X 1000 -k 2 --no-mixed --no-discordant")
    (default for bowtie: "-X 1000 -S -k 1 -m 1 --best --strata --chunkmbs 64 -n 1")

    
    num_procs is an integer indicating how many num_procs you'd like to run this function over
    
    trim_reads is a boolean indicating that you want to have reads trimmed by cutadapt. This option is currently
        not compatible with "in_mem" mode.

    path_to_cutadapt is the path to the cutadapt execuatable. Otherwise this is assumed to be in your
        path.

    adapter_seq_R1:
        Sequence of an adapter that was ligated to the 3' end of read 1. The adapter itself and anything that follows is
        trimmed.

    adapter_seq_R2:
        Sequence of an adapter that was ligated to the 3' end of read 2. The adapter itself and anything that follows is
        trimmed.

    max_adapter_removal indicates the maximum number of times to try to remove adapters. Useful when an adapter 
        gets appended multiple times.
        
    overlap_length is the minimum overlap length. If the overlap between the read and the adapter is shorter than 
        LENGTH, the read is not modified. This reduces the no. of bases trimmed purely due to short random adapter matches.

    zero_cap causes negative quality values to be set to zero (workaround to avoid segmentation faults in BWA).
        
    error_rate is the maximum allowed error rate (no. of errors divided by the length of the matching region) 
        (default: 0.1)
    
    min_qual_score allows you to trim low-quality ends from reads before adapter removal. The algorithm is the same as the 
        one used by BWA (Subtract CUTOFF from all qualities; compute partial sums from all indices to the end of the
        sequence; cut sequence at the index at which the sum is minimal).
        
    min_read_len indicates the minimum length a read must be to be kept. Reads that are too short even before adapter removal
        are also discarded. In colorspace, an initial primer is not counted. It is not recommended to change this value in paired-end processing because it may results in the situation that one of the two reads in a pair is discarded. And this will lead to error.
        
    sig_cutoff is a float indicating the adjusted p-value cutoff you wish to use for determining whether or not
        a site is methylated

    min_cov is an integer indicating the minimum number of reads for a site to be tested.

    binom_tests indicates that you'd like to use a binomial test, rather than the alternative method outlined here
        https://bitbucket.org/schultzmattd/methylpy/wiki/Methylation%20Calling

    keep_temp_files is a boolean indicating that you'd like to keep the intermediate files generated
        by this function. This can be useful for debugging, but in general should be left False.
        
    num_reads is an integer indicating how many reads you'd like to process at a time. Breaking up the processing this
        way can reduce the size of some of the intermediate files that are created. Use this if you're strapped for space.
        
    save_space indicates whether or not you'd like to perform read collapsing right after mapping or once all the libraries have
        been mapped. If you wait until after everything has been mapped, the collapsing can be parallelized. Otherwise the collapsing
        will have to be done serially. The trade-off is that you must keep all the mapped files around, rather than deleting them as they
        are processed, which can take up a considerable amount of space. It's safest to set this to True.
    
    bowtie2 specifies whether to use the bowtie2 aligner instead of bowtie
    
    sort_mem is the parameter to pass to unix sort with -S/--buffer-size command
    
    in_mem: does not write files as part of run_mapping, keeps all input/output in memory (This mode is unavailable now)
    
    path_to_output is the path to a directory where you would like the output to be stored. The default is the
        same directory as the input fastqs.
    min_base_quality is an integer indicating the minimum PHRED quality score for a base to be included in the
        mpileup file (and subsequently to be considered for methylation calling)    

    path_to_picard is the path to MarkDuplicates jar from picard. Default is false indicating that you
        don't want to use this jar for duplication removal
    """

    #Need to work out the in_mem mode later
    in_mem=False

    #Default bowtie option
    if len(bowtie_options) == 0:
        if bowtie2:
            bowtie_options = ["-X 1000","-k 2","--no-discordant","--no-mixed"]
        else:
            bowtie_options=["-X 1000","-S","-k 1","-m 1","--best","--strata",
                            "--chunkmbs 3072","-n 1","-e 100"]

    #Figure out all the correct quality options based on the offset or CASAVA version given
    try:
        quality_version = float(quality_version)
        if quality_version >= 1.3 and quality_version<1.8:
            quality_base = 64
            bowtie_options.append("--phred64-quals")
        elif quality_version >= 1.8:
            quality_base = 33
            bowtie_options.append("--phred33-quals")
        elif quality_version < 1.3:
            quality_base = 64
            bowtie_options.append("--solexa-quals")
        else:
            sys.exit("Unrecognized quality_version. Either enter the CASAVA version used to generate the fastq files, or enter the base quality offset (either 33 or 64).")
    except:
        if int(quality_version) == 64:
            quality_base = 64
            bowtie_options.append("--phred64-quals")
        elif int(quality_version) == 33:
            quality_base = 33
            bowtie_options.append("--phred33-quals")
        else:
            sys.exit("Unrecognized quality_version. Either enter the CASAVA version used to generate the fastq files, or enter the base quality offset (either 33 or 64).")

    #This code allows the user to supply paths with "*" in them rather than listing out every single
    #file
    #just to avoid any paths missing a slash. It's ok to add an extra slash if
    #there's already one at the end
    if len(path_to_samtools)!=0:
        path_to_samtools +="/"
    if len(path_to_bowtie)!=0:
        path_to_bowtie+="/"
    if len(path_to_output) !=0:
        path_to_output+="/"
        
    if sort_mem:
        if sort_mem.find("-S") == -1:
            sort_mem = " -S " + sort_mem
    else:
        sort_mem = ""
        
    expanded_file_list_R1 = []
    expanded_file_list_R2 = []
    expanded_library_list = []
    
    total_reads = 0
    total_unique = 0
    total_clonal = 0
    #Get expanded file list and the number of total reads
    total_reads_R1 = 0
    total_reads_R2 = 0
    #R1
    for path,library in zip(read1_files,libraries):
        glob_list = glob.glob(path)
        for filen in glob_list:            
            expanded_file_list_R1.append(filen)
            expanded_library_list.append(library)
            total_reads_R1 += int(subprocess.check_output(['wc', '-l', filen]).split(' ')[0]) / 4
    #R2
    for path,library in zip(read2_files,libraries):
        glob_list = glob.glob(path)
        for filen in glob_list:
            expanded_file_list_R2.append(filen)
            total_reads_R2 += int(subprocess.check_output(['wc', '-l', filen]).split(' ')[0]) / 4
    #Check if there are same number of reads in read 1 and read 2
    if total_reads_R1 != total_reads_R2:
        print "There are different numbers of read 1 and read 2!"
        exit()
    else:
        total_reads = total_reads_R1
    print_checkpoint("There are " + str(total_reads) + " total input read pairs")
    #print expanded_file_list_R1
    #print expanded_file_list_R2
    #print expanded_library_list

    #Processing
    for current_library in set(libraries):
        library_read1_files = [filen for filen,library in zip(expanded_file_list_R1,expanded_library_list) if library == current_library]
        library_read2_files = [filen for filen,library in zip(expanded_file_list_R2,expanded_library_list) if library == current_library]

        #deal with actual filename rather than path to file
        if in_mem:
            total_unique += Mrun_mapping_pe(
                current_library,library_read1_files,library_read2_files,
                sample,forward_reference,reverse_reference,reference_fasta,
                path_to_samtools=path_to_samtools,path_to_bowtie=path_to_bowtie,
                bowtie_options=bowtie_options,num_procs=num_procs,
                trim_reads=trim_reads,path_to_cutadapt=path_to_cutadapt,
                adapter_seq = adapter_seq,
                #adapter_seq_R1 = adapter_seq_R1, adapter_seq_R2 = adapter_seq_R2,
                max_adapter_removal=max_adapter_removal,
                overlap_length=overlap_length,zero_cap=zero_cap,quality_base=quality_base,
                error_rate=error_rate,min_qual_score=min_qual_score,min_read_len=min_read_len,
                keep_temp_files=keep_temp_files,num_reads=num_reads, bowtie2=bowtie2, 
                sort_mem=sort_mem, path_to_output=path_to_output)
        else:
            total_unique += run_mapping_pe(
                current_library,library_read1_files,library_read2_files,sample,
                forward_reference,reverse_reference,reference_fasta,
                path_to_samtools=path_to_samtools,path_to_bowtie=path_to_bowtie,
                bowtie_options=bowtie_options,num_procs=num_procs,
                trim_reads=trim_reads,path_to_cutadapt=path_to_cutadapt,
                adapter_seq_R1 = adapter_seq_R1, adapter_seq_R2 = adapter_seq_R2,
                max_adapter_removal=max_adapter_removal,overlap_length=overlap_length,
                zero_cap=zero_cap,quality_base=quality_base,error_rate=error_rate,
                min_qual_score=min_qual_score,min_read_len=min_read_len,
                keep_temp_files=keep_temp_files,num_reads=num_reads, bowtie2=bowtie2, 
                sort_mem=sort_mem, path_to_output=path_to_output,save_space=save_space)

        if save_space == True:
            pool = multiprocessing.Pool(num_procs)
            #I avoid sorting by chromosome or strand because it's not strictly necessary and this speeds up 
            #the sorting. I had to add some special logic in the processing though
            processed_library_files = glob.glob(path_to_output+sample+"_"+str(current_library)+"_*_no_multimap_*")
            total_clonal += collapse_clonal_reads(reference_fasta,path_to_samtools,num_procs,sample,current_library,path_to_files=path_to_output,sort_mem=sort_mem,remove_clonal=remove_clonal,path_to_picard=path_to_picard)
    if save_space == True:
        print_checkpoint("There are " + str(total_unique) + " uniquely mapping reads, " + str(float(total_unique) / total_reads*100) + " percent remaining")
    else:
        print_checkpoint("Finding multimappers")
        total_unique = 0
        merge_results= []
        sorted_files = []
        for library in set(libraries):
            pool = multiprocessing.Pool(num_procs)
            library_sorted_files = glob.glob(sample+"_"+str(library)+"_[0-9]*_sorted_[0-9]*")
            sorted_files.extend(library_sorted_files)
            sort_results = []
            for filename in library_sorted_files:
                sort_results.append(pool.apply_async(subprocess.check_call,(shlex.split("env LC_COLLATE=C sort" + sort_mem + " -t '\t' -k 1 -o "+filename+" "+filename),)))
            #ensure the sorting has finished
            for result in sort_results:
                result.get()
            #Need to add the _1 for a bogus chunk ID
            merge_results.append(pool.apply_async(merge_sorted_multimap_pe,(sorted_files,sample+"_"+str(library)+"_1")))
            #total_unique += merge_sorted_multimap(sorted_files,sample+"_"+library,num_procs)
        pool.close()
        pool.join()
        subprocess.check_call(shlex.split("rm "+" ".join(sorted_files)))
        for result in merge_results:
            total_unique+=result.get()
            
        print_checkpoint("There are " + str(total_unique) + " uniquely mapping reads, " + str(float(total_unique) / total_reads*100) + " percent remaining")
        
        
        pool = multiprocessing.Pool(num_procs)
        #I avoid sorting by chromosome or strand because it's not strictly necessary and this speeds up 
        #the sorting. I had to add some special logic in the processing though
        processed_library_files = glob.glob(path_to_output+sample+"_"+str(current_library)+"_*_no_multimap_*")
        for filen in processed_library_files:
            cmd = shlex.split("sort" + sort_mem + " -t '\t' -k 4n -k 8n -o "+filen+" "+filen)
            pool.apply_async(subprocess.check_call,(cmd,))
        pool.close()
        pool.join()
        
        clonal_results = []
        pool = multiprocessing.Pool(num_procs)

        for library in set(libraries):
            clonal_results.append(pool.apply_async(collapse_clonal_reads,(reference_fasta,path_to_samtools,num_procs,sample,library), {"sort_mem":sort_mem,"path_to_files":path_to_output,"remove_clonal":remove_clonal,"path_to_picard":path_to_picard}))
        pool.close()
        pool.join()
        
        for result in clonal_results:
            total_clonal += result.get()
            
    if remove_clonal == True:
        print_checkpoint("There are " + str(total_clonal/2) + " non-clonal reads, " + str(float(total_clonal/2) / total_reads*100) + " percent remaining") 
        library_files = [path_to_output+sample+"_processed_reads_"+str(library)+"_no_clonal.bam" for library in set(libraries)]
        if len(library_files) > 1:
            merge_bam_files(library_files,path_to_output+sample+"_processed_reads_no_clonal.bam",path_to_samtools)
            subprocess.check_call(shlex.split("rm "+" ".join(library_files)))
        else:
            subprocess.check_call(shlex.split("mv "+library_files[0]+" "+path_to_output+sample+"_processed_reads_no_clonal.bam"))
    else:
        library_files = [path_to_output+sample+"_processed_reads_"+str(library)+".bam" for library in set(libraries)]
        if len(library_files) > 1:
            merge_bam_files(library_files,path_to_output+sample+"_processed_reads.bam",path_to_samtools)
            subprocess.check_call(shlex.split("rm "+" ".join(library_files)))
        else:
            subprocess.check_call(shlex.split("mv "+library_files[0]+" "+path_to_output+sample+"_processed_reads.bam"))

    #Calling methylated sites
    print_checkpoint("Begin calling mCs")
    if remove_clonal == True:
        call_methylated_sites_pe(sample+"_processed_reads_no_clonal.bam",sample,reference_fasta,unmethylated_control,quality_version,sig_cutoff=sig_cutoff,num_procs=num_procs,min_cov=min_cov,binom_test=binom_test,bh=bh,sort_mem=sort_mem,path_to_files=path_to_output,path_to_samtools=path_to_samtools,min_base_quality=min_base_quality)
    else:
        call_methylated_sites_pe(sample+"_processed_reads.bam",sample,reference_fasta,unmethylated_control,quality_version,sig_cutoff=sig_cutoff,num_procs=num_procs,min_cov=min_cov,binom_test=binom_test,bh=bh,sort_mem=sort_mem,path_to_files=path_to_output,path_to_samtools=path_to_samtools,min_base_quality=min_base_quality)
    print_checkpoint("Done")


def run_mapping_pe(current_library,library_read1_files,library_read2_files,
                   sample,forward_reference,reverse_reference,reference_fasta,
                   path_to_samtools="",path_to_bowtie="",
                   bowtie_options=[],
                   num_procs=1,trim_reads=True,path_to_cutadapt="",
                   adapter_seq_R1 = "AGATCGGAAGAGCACACGTCTGAAC",
                   adapter_seq_R2 = "AGATCGGAAGAGCGTCGTGTAGGGA",
                   max_adapter_removal=None,overlap_length=None,zero_cap=None,
                   quality_base=None,error_rate=None,min_qual_score=10,
                   min_read_len=30,keep_temp_files=False,num_reads=-1, 
                   bowtie2=True, sort_mem="500M",path_to_output="",
                   save_space=True):
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
    
    path_to_bowtie is a string indicating the path to the folder in which bowtie resides. Bowtie
        is assumed to be in your path if this option isn't used
            
    bowtie_options is a list of strings indicating options you'd like passed to bowtie2 (or bowtie)
    
    num_procs is an integer indicating how many num_procs you'd like to run this function over
    
    trim_reads is a boolean indicating that you want to have reads trimmed by cutadapt

    path_to_cutadapt is the path to the cutadapt execuatable. Otherwise this is assumed to be in your
        path.

    adapter_seq_R1:
        Sequence of an adapter that was ligated to the 3' end of read 1. The adapter itself and anything that follows is
        trimmed.

    adapter_seq_R2:
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
        
    num_reads is an integer indicating how many reads you'd like to process at a time. Breaking up the processing this
        way can reduce the size of some of the intermediate files that are created. Use this if you're strapped for space.
    
    bowtie2 specifies whether to use the bowtie2 aligner instead of bowtie
    
    sort_mem is the parameter to pass to unix sort with -S/--buffer-size command
    """
    total_unique = 0
    chunk = 1
    file_name = sample+"_"+str(current_library)+"_"+str(chunk)
    file_path = path_to_output+file_name
    save_point_r1 = (False,0)
    save_point_r2 = (False,0)
    while save_point_r1 != False:#We expect R1 and R2 have the same number of lines
        print_checkpoint("Begin chunk "+str(chunk)+" for library "+str(current_library))
        #Split files
        print_checkpoint("Begin splitting reads for "+file_name)
        save_point_r1 = split_fastq_file(num_procs,library_read1_files,file_path+"_R1_split_",
                                         num_reads=num_reads,current_file=save_point_r1[0],current_pointer=save_point_r1[1])
        save_point_r2 = split_fastq_file(num_procs,library_read2_files,file_path+"_R2_split_",
                                         num_reads=num_reads,current_file=save_point_r2[0],current_pointer=save_point_r2[1])
        if trim_reads:
            #Trimming
            print_checkpoint("Begin trimming reads for "+file_name)  
            quality_trim_pe(
                inputf_R1=[file_path+"_R1_split_"+str(i) for i in xrange(0,num_procs)],
                outputf_R1=[file_path+"_R1_split_trimmed_"+str(i) for i in xrange(0,num_procs)],
                inputf_R2=[file_path+"_R2_split_"+str(i) for i in xrange(0,num_procs)],
                outputf_R2=[file_path+"_R2_split_trimmed_"+str(i) for i in xrange(0,num_procs)],
                adapter_seq_R1=adapter_seq_R1,
                adapter_seq_R2=adapter_seq_R2,
                error_rate=error_rate,
                quality_base = quality_base,
                min_qual_score=min_qual_score,
                min_read_len=min_read_len,
                format="fastq",
                num_procs=num_procs,
                max_adapter_removal=max_adapter_removal,
                overlap_length=overlap_length,
                zero_cap=zero_cap,
                path_to_cutadapt=path_to_cutadapt)        
            #Conversion
            print_checkpoint("Begin converting reads for "+file_name)
            pool = multiprocessing.Pool(num_procs)#R1
            for inputf,output in zip([file_path+"_R1_split_trimmed_"+str(i) for i in xrange(0,num_procs)],
                                     [file_path+"_R1_split_trimmed_converted_"+str(i) for i in xrange(0,num_procs)]):
                pool.apply_async(convert_reads_pe,(inputf,output))
            for inputf,output in zip([file_path+"_R2_split_trimmed_"+str(i) for i in xrange(0,num_procs)],
                                 [file_path+"_R2_split_trimmed_converted_"+str(i) for i in xrange(0,num_procs)]):
                pool.apply_async(convert_reads_pe,(inputf,output,True))
            pool.close()
            pool.join()
            subprocess.check_call(
                shlex.split("rm "+" ".join([file_path+"_R1_split_trimmed_"+str(i) for i in xrange(0,num_procs)]))
            )
            subprocess.check_call(
                shlex.split("rm "+" ".join([file_path+"_R2_split_trimmed_"+str(i) for i in xrange(0,num_procs)]))
            )        
            #Run bowtie
            print_checkpoint("Begin Running Bowtie for "+file_name)
            total_unique += run_bowtie_pe([file_path+"_R1_split_trimmed_converted_"+str(i) for i in xrange(0,num_procs)],
                                          [file_path+"_R2_split_trimmed_converted_"+str(i) for i in xrange(0,num_procs)],
                                          forward_reference,reverse_reference,file_path,options=bowtie_options,
                                          path_to_bowtie=path_to_bowtie,num_procs=num_procs,
                                          keep_temp_files=keep_temp_files, bowtie2=bowtie2, sort_mem=sort_mem,
                                          save_space=save_space)
        else:
            print_checkpoint("No trimming applied on reads")  
            #Conversion
            print_checkpoint("Begin converting reads for "+file_name)
            pool = multiprocessing.Pool(num_procs)#R1
            for inputf,output in zip([file_path+"_R1_split_"+str(i) for i in xrange(0,num_procs)],
                                     [file_path+"_R1_split_converted_"+str(i) for i in xrange(0,num_procs)]):
                pool.apply_async(convert_reads_pe,(inputf,output))
            for inputf,output in zip([file_path+"_R2_split_"+str(i) for i in xrange(0,num_procs)],
                                 [file_path+"_R2_split_converted_"+str(i) for i in xrange(0,num_procs)]):
                pool.apply_async(convert_reads_pe,(inputf,output,True))
            pool.close()
            pool.join()
            #Run bowtie
            print_checkpoint("Begin Running Bowtie for "+file_name)
            total_unique += run_bowtie_pe([file_path+"_R1_split_converted_"+str(i) for i in xrange(0,num_procs)],
                                          [file_path+"_R2_split_converted_"+str(i) for i in xrange(0,num_procs)],
                                          forward_reference,reverse_reference,file_path,options=bowtie_options,
                                          path_to_bowtie=path_to_bowtie,num_procs=num_procs,
                                          keep_temp_files=keep_temp_files, bowtie2=bowtie2, sort_mem=sort_mem,
                                          save_space=save_space)
        chunk+=1

    subprocess.check_call(shlex.split("rm "+" ".join([file_path+"_R1_split_"+str(i) for i in xrange(0,num_procs)])))
    subprocess.check_call(shlex.split("rm "+" ".join([file_path+"_R2_split_"+str(i) for i in xrange(0,num_procs)])))
    
    return total_unique

def run_bowtie_pe(library_read1_files,library_read2_files,forward_reference,reverse_reference,prefix,options="",
                  path_to_bowtie="",num_procs=1,keep_temp_files=False, bowtie2=True, sort_mem="500M",save_space=True):
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
    
    options is a list of strings indicating options you'd like passed to bowtie2 (or bowtie)
    
    path_to_bowtie is a string indicating the path to the folder in which bowtie resides. Bowtie
        is assumed to be in your path if this option isn't used
    
    num_procs is an integer indicating the number of processors you'd like used for removing multi
        mapping reads and for bowtie mapping
    
    keep_temp_files is a boolean indicating that you'd like to keep the intermediate files generated
        by this function. This can be useful for debugging, but in general should be left False.
    
    bowtie2 specifies whether to use the bowtie2 aligner instead of bowtie
    
    sort_mem is the parameter to pass to unix sort with -S/--buffer-size command
    """
    if sort_mem:
        if sort_mem.find("-S") == -1:
            sort_mem = " -S " + sort_mem
    else:
        sort_mem = ""
    if " ".join(options).find(" -p ") == -1:
        options.append("-p "+str(num_procs))
    if bowtie2:
        args = [path_to_bowtie+"bowtie2"]
        args.extend(options)
        args.append("--norc")
        args.append("-x "+forward_reference)
        args.append("-1 "+",".join(library_read1_files))
        args.append("-2 "+",".join(library_read2_files))
        args.append("-S "+prefix+"_forward_strand_hits.sam")
    else:
        args = [path_to_bowtie+"bowtie"]
        args.extend(options)
        args.append("--norc")
        args.append(forward_reference)
        args.append("-1 "+",".join(library_read1_files))
        args.append("-2 "+",".join(library_read2_files))
        args.append(prefix+"_forward_strand_hits.sam")
    subprocess.check_call(shlex.split(" ".join(args)))
    print_checkpoint("Processing forward strand hits")
    find_multi_mappers_pe(prefix+"_forward_strand_hits.sam",prefix,num_procs=num_procs,keep_temp_files=keep_temp_files)
    if bowtie2:
        args = [path_to_bowtie+"bowtie2"]
        args.extend(options)
        args.append("--nofw")
        args.append("-x "+reverse_reference)
        args.append("-1 "+",".join(library_read1_files))
        args.append("-2 "+",".join(library_read2_files))
        args.append("-S "+prefix+"_reverse_strand_hits.sam")
    else:
        args = [path_to_bowtie+"bowtie"]
        args.extend(options)
        args.append("--nofw")
        args.append(reverse_reference)
        args.append("-1 "+",".join(library_read1_files))
        args.append("-2 "+",".join(library_read2_files))
        args.append(prefix+"_reverse_strand_hits.sam")
    subprocess.check_call(shlex.split(" ".join(args)))
    
    print_checkpoint("Processing reverse strand hits")
    sam_header = find_multi_mappers_pe(prefix+"_reverse_strand_hits.sam",prefix,num_procs=num_procs,append=True,keep_temp_files=keep_temp_files)
    if keep_temp_files==False:
        subprocess.check_call(shlex.split("rm "+" ".join(library_read1_files+library_read2_files)))
    if save_space == True:
        pool = multiprocessing.Pool(num_procs)
        for file_num in xrange(0,num_procs):
            pool.apply_async(subprocess.check_call,(shlex.split("env LC_COLLATE=C sort" + sort_mem + " -t '\t' -k 1 -o "+prefix+"_sorted_"+str(file_num)+" "+prefix+"_sorted_"+str(file_num)),))
        pool.close()
        pool.join()
        print_checkpoint("Finding multimappers")

        total_unique = merge_sorted_multimap_pe([prefix+"_sorted_"+str(file_num) for file_num in xrange(0,num_procs)],prefix)
        subprocess.check_call(shlex.split("rm "+" ".join([prefix+"_sorted_"+str(file_num) for file_num in xrange(0,num_procs)])))
        return total_unique
    else:
        return 0


def find_multi_mappers_pe(inputf,output,num_procs=1,keep_temp_files=False,append=False):
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
        if line[0] == "@":
            continue
        
        fields = line.split("\t")

        #To deal with the way chromosomes were named in some of our older references
        fields[2] = fields[2].replace("_f","")
        fields[2] = fields[2].replace("_r","")

        if int(fields[1]) & 2 != 0:
            header = fields[0].split("!")
            #BIG ASSUMPTION!! NO TABS IN FASTQ HEADER LINES EXCEPT THE ONES I ADD!
            if (int(fields[1]) & 16) == 16:
                strand = "-"
            else:
                strand = "+"
            if (int(fields[1]) & 128) == 128:
                is_R2 = True
            else:
                is_R2 = False
            seq = decode_converted_positions(fields[9],header[-1],strand,is_R2)
            file_handles[cycle.next()].write(" ".join(header[:-1])+"\t"+"\t".join(fields[1:9])+"\t"+seq+"\t"+"\t".join(fields[10:]))
            #file_handles[cycle.next()].write("\t".join(fields[0:9])+"\t"+seq+"\t"+"\t".join(fields[10:]))
    f.close()
    if keep_temp_files == False:
        subprocess.check_call(shlex.split("rm "+inputf))
        pass
    for file_num in xrange(0,num_procs):
        file_handles[file_num].close()

def merge_sorted_multimap_pe(files, output):
    """
    This function takes the files from find_multi_mappers and outputs the uniquely mapping reads
    
    files is a list of filenames containing the output of find_multi_mappers
    
    output is a prefix you'd like prepended to the file containing the uniquely mapping reads
        This file will be named as <output>+"_no_multimap_"+<index_num>
    """
    lines = {}
    fields = {}
    output_handles = {}
    file_handles = {}
    
    total_unique = 0
    count= 0
    cycle = itertools.cycle(range(0,len(files)))
    
    for index,filen in enumerate(files):
        output_handles[index] = open(output+"_no_multimap_"+str(index),'a')
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
                #Need to modify this in order to deal with PE data
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
            index = cycle.next()
            #Output
            output_handles[index].write(current_line_1)
            output_handles[index].write(current_line_2)
            output_handles[index].flush()
            total_unique += 1
            
    for index,filen in enumerate(files):
        output_handles[index].close()
        file_handles[filen].close()

    ##Yupeng debug
#    exit()
    
    return total_unique

def convert_reads_pe(inputf,output,is_R2=False):
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
    encoding = encode_converted_positions(seq,is_R2=is_R2)
    
    if is_R2 == False:
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
            encoding = encode_converted_positions(seq,is_R2=is_R2)
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
            encoding = encode_converted_positions(seq,is_R2=is_R2)
    f.close()
    g.close()

def quality_trim(inputf, output = None, quality_base = None, min_qual_score = None, min_read_len = None, 
                 adapter_seq = "AGATCGGAAGAGCACACGTCTG", num_procs = 1, format = None, error_rate = None, 
                 max_adapter_removal = None, overlap_length = None, zero_cap = False, path_to_cutadapt = ""):
    """
    Information from cutadapt documentation:
    format:
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
        if isinstance(inputf, basestring):
            inputf = [inputf]
        else:
            sys.exit("input must be a list of strings")
    if not isinstance(output, list):
        if isinstance(output, basestring):
            output = [output]
        else:
            sys.exit("output must be a list of strings")
            
    if len(output) != len(inputf):
        sys.exit("Must provide an equal number of input and output files")
    base_cmd = path_to_cutadapt
    options = ""
    if zero_cap:
        zero = "-z "
    else:
        zero = ""  
    
    if format:
        options += " -f " + format
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
    pool = multiprocessing.Pool(num_procs)
    #adapter trimming
    for current_input,current_output in zip(inputf,output):
        if output:
            options += " -o " + current_output + " "
        pool.apply_async(subprocess.check_call,(base_cmd + options + current_input,),{"shell":True})
    pool.close()
    pool.join()

def quality_trim_pe(inputf_R1, outputf_R1,inputf_R2, outputf_R2,quality_base = None, min_qual_score = 10,
                    min_read_len = 30,adapter_seq_R1 = "AGATCGGAAGAGCACACGTCTGAAC",
                    adapter_seq_R2 = "AGATCGGAAGAGCGTCGTGTAGGGA",num_procs = 1, format = None,
                    error_rate = None, max_adapter_removal = None,overlap_length = None, zero_cap = False,
                    path_to_cutadapt = ""):
    """
    Information from cutadapt documentation:
    format:
        Input file format; can be either 'fasta', 'fastq' or 'sra-fastq'. Ignored when reading csfasta/qual files
        (default: auto-detect from file name extension).

    inputf_R1,inputf_R2:
        list of filenames for read 1 and read 2 respectively

    outputf_R1,outputf_R2:
        Write the modified sequences to these files instead of standard output and send the summary report to
        standard output. The format is FASTQ if qualities are available, FASTA otherwise. outputf_R1 and outputf_R2
        specify the output filenames of read 1 and read 2 respectively.

    adapter_seq_R1:
        Sequence of an adapter that was ligated to the 3' end of read 1. The adapter itself and anything that follows is
        trimmed.

    adapter_seq_R2:
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
                 
    if not isinstance(inputf_R1, list):
        if isinstance(inputf_R1, basestring):
            inputf = [inputf_R1]
        else:
            sys.exit("inputf_R1 must be a list of strings")
    if not isinstance(inputf_R2, list):
        if isinstance(inputf_R2, basestring):
            inputf = [inputf_R2]
        else:
            sys.exit("inputf_R2 must be a list of strings")

    if not isinstance(outputf_R1, list):
        if isinstance(outputf_R1, basestring):
            output = [outputf_R1]
        else:
            sys.exit("outputf_R1 must be a list of strings")
    if not isinstance(outputf_R2, list):
        if isinstance(outputf_R2, basestring):
            output = [outputf_R2]
        else:
            sys.exit("outputf_R2 must be a list of strings")            
            
    if len(outputf_R1) != len(inputf_R2) or len(outputf_R1) != len(outputf_R1) or len(outputf_R1) != len(outputf_R2):
        sys.exit("Must provide an equal number of input and output files")

    base_cmd = path_to_cutadapt
    options = ""
    if zero_cap:
        zero = "-z "
    else:
        zero = ""  
    
    if format:
        options += " -f " + format
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
    options += " -a " + adapter_seq_R1
    options += " -A " + adapter_seq_R2
    options += " " + zero
    pool = multiprocessing.Pool(num_procs)
    #adapter trimming
    for current_input_R1,current_output_R1,current_input_R2,current_output_R2 in zip(inputf_R1,outputf_R1,inputf_R2,outputf_R2):
        options += " -o " + current_output_R1 + " " + " -p " + current_output_R2 + " "
        pool.apply_async(subprocess.check_call,(base_cmd + options + current_input_R1 + " " + current_input_R2,),{"shell":True})
    pool.close()
    pool.join()


def flip_R2_strand(input_file,output_file,path_to_samtools=""):
    """
    This function flips the strand of all read2s (R2s) of mapped paired-end
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

    
def call_methylated_sites_pe(inputf, sample, reference, control,quality_version,sig_cutoff=.01,num_procs = 1,
                             min_cov=1,binom_test=True,min_mc=0,path_to_samtools="",sort_mem="500M",bh=False,
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
    
    quality_version is either an integer indicating the base offset for the quality scores or a float indicating
        which version of casava was used to generate the fastq files.
    
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

    #Flip the strand of read 2 (R2) and create a new bam file
    ##input: sample+"_processed_reads_no_clonal.bam"
    ##output: sample+"_processed_reads_no_clonal_flipped.bam"
    print_checkpoint("Begin flipping the strand of R2 reads")
    flip_R2_strand(input_file = inputf,
                   output_file = inputf+".R2flipped.bam",
                   path_to_samtools=path_to_samtools)

        
    #Call methylated sites
    call_methylated_sites(inputf+".R2flipped.bam", sample, reference, control,quality_version,sig_cutoff,num_procs,
                          min_cov,binom_test,min_mc,path_to_samtools,sort_mem,bh,
                          path_to_files,min_base_quality)
    
    #Remove intermediate bam file
    try:
        subprocess.check_call(shlex.split("rm -f "+inputf+".R2flipped.bam"+
                                          " "+inputf+".R2flipped.bam.bai"))
    except:
        pass

    
def call_methylated_sites(inputf, sample, reference, control,quality_version,sig_cutoff=.01,num_procs = 1,
                          min_cov=1,binom_test=True,min_mc=0,path_to_samtools="",sort_mem="500M",bh=False,path_to_files="",min_base_quality=1):

    """
    inputf is the path to a bam file that contains mapped bisulfite sequencing reads
    
    sample is the name you'd like for the allc files. The files will be named like so:
        allc_<sample>_<chrom>.tsv
    
    reference is the path to a samtools indexed fasta file
    
    control is the name of the chromosome/region that you want to use to estimate the non-conversion rate of your 
        sample, or the non-conversion rate you'd like to use. Consequently, control is either a string, or a decimal
        If control is a string then it should be in the following format: "chrom:start-end". 
        If you'd like to specify an entire chromosome simply use "chrom:"
    
    quality_version is either an integer indicating the base offset for the quality scores or a float indicating
        which version of casava was used to generate the fastq files.
    
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
    #Figure out all the correct quality options based on the offset or CASAVA version given
    try:
        quality_version = float(quality_version)
        if quality_version >= 1.3 and quality_version<1.8:
            quality_base = 64
            mpileup_quality = "-6"
        elif quality_version >= 1.8:
            quality_base = 33
            mpileup_quality = ""
        elif quality_version < 1.3:
            quality_base = 64
            mpileup_quality = "-6"
    except:
        if int(quality_version) == 64:
            quality_base = 64
            mpileup_quality = "-6"
        elif int(quality_version) == 33:
            quality_base = 33
            mpileup_quality = ""
        else:
            sys.exit("Unrecognized quality_version. Either enter the CASAVA version used to generate the fastq files, or enter the base quality offset (either 33 or 64).")

    try:
        num_procs = int(num_procs)
    except:
        sys.exit("num_procs must be an integer")
    if len(path_to_files)!=0:
        path_to_files+="/"

    mc_class_counts = {}
    for first in ["A","T","C","G","N"]:
        for second in ["A","T","C","G","N"]:
            mc_class_counts["C"+first+second]=0

    #figure out non-conversion rate if it isn't given
    try:
        non_conversion = float(control)
    except:
        if control.find(":") == -1:
            sys.exit("control MUST have a colon in it. Either chr: or chr:start-end")
        print_checkpoint("Begin calculating non-conversion rate for "+control)
        try:
            open(path_to_files+sample+"_mpileup_output_"+control+".tsv",'r')
        except IOError:
            try:
                #make sure bam file is indexed
                open(path_to_files+inputf+".bai",'r')
            except:
                print_checkpoint("Input not indexed. Indexing...")
                subprocess.check_call(shlex.split(path_to_samtools+"samtools index "+path_to_files+inputf))
            with open(path_to_files+sample+"_mpileup_output_"+control+".tsv",'w') as f:
                subprocess.check_call(shlex.split(path_to_samtools+"samtools mpileup -Q "+str(min_base_quality)+" -B "+mpileup_quality+" -r "+control+" -f "+reference+" "+path_to_files+inputf),stdout=f)
        else:
            print(sample+"_mpileup_output_"+control+".tsv already exists, using it for calculations.")
        f = open(path_to_files+sample+"_mpileup_output_"+control+".tsv" ,'r')

        unconverted_c = 0
        converted_c = 0
        total_phred = 0
        total_bases = 0
        total_control_sites = 0
        total_unconverted_c =0
        total_converted_c = 0
        pvalue_lookup={}
        control_sites = []
        for line in f:
            line = line.rstrip()
            fields = line.split("\t")
            if fields[3] != "0":
                total_phred += sum([ord(i)-quality_base for i in fields[5]]) 
                total_bases += len(fields[5])
                
                if fields[2] == "C":
                    unconverted_c = fields[4].count(".")
                    converted_c = fields[4].count("T")
                elif fields[2] == "G":
                    unconverted_c = fields[4].count(",")
                    converted_c = fields[4].count("a")
                else:
                    continue
                total_unconverted_c += unconverted_c
                total_converted_c += converted_c
                if binom_test==False and min_mc==0:
                    if unconverted_c+converted_c > 0:
                        total_control_sites += 1
                        control_sites.append((unconverted_c,unconverted_c+converted_c))
        f.close()
        #compute pvalues for control genome. Have to read through it completely twice to do this
        min_pvalue = 1
        avg_qual = total_phred / total_bases
        seq_error = 10 ** (avg_qual / -10)
        if seq_error > 1.0 or seq_error < 0.0:
            sys.exit("One of your quality values corresponds to a sequence error rate of "+str(seq_error)+". These error rates have to be between 0 and 1. Are you sure you chose the correct CASAVA version?")
        non_conversion = total_unconverted_c / float(total_converted_c + total_unconverted_c)
        control_pvalues =[]
        f = open(path_to_files+sample+"_mpileup_output_"+control+".tsv" ,'r')
        for line in f:

            fields = line.split("\t")
            if fields[2] != "C" and fields[2] != "G":
                continue
            if fields[2] == "C":
                unconverted_c = fields[4].count(".")
                converted_c = fields[4].count("T")
            elif fields[2] == "G":
                unconverted_c = fields[4].count(",")
                converted_c = fields[4].count("a")
            control_pvalues.append(sci.binom.sf(unconverted_c-1,unconverted_c+converted_c,non_conversion+seq_error))
            if binom_test==False and min_mc==0:
                coverage = unconverted_c+converted_c
                if (fields[2] == "C" or fields[2] == "G") and coverage > 0:
                    if coverage not in pvalue_lookup:
                        pvalue_lookup[coverage] = get_pvalue_dist(coverage,control_sites) 
                    if pvalue_lookup[coverage][unconverted_c] < min_pvalue:
                        min_pvalue=pvalue_lookup[coverage][unconverted_c]
        if binom_test==False and min_mc==0:
            print "The minimum p-value in the control genome is "+str(min_pvalue)
        print("\tThe non-conversion rate is "+str(non_conversion*100)+"%")
        print("\tThe estimated sequencing error rate is: "+str(seq_error))
        non_conversion+=seq_error
        f.close()
    #Avoids matching similarly named files
    allc_files = [filename for filename in glob.glob("allc_"+sample+"_*.tsv") if len(re.findall("^allc_"+sample+"_[0-9a-zA-Z]{1,3}.tsv$",filename))>0]
    if len(allc_files) != 0:
        print_checkpoint("allc files exist. Using these files for calling methylated cytosines.")
        print_checkpoint("Begin binomial tests")
        if num_procs > 1:
            results = []
            pool=multiprocessing.Pool(num_procs)
            for chunk in allc_files:
                results.append(pool.apply_async(allc_run_binom_tests,(chunk,non_conversion,min_cov), {"sort_mem":sort_mem}))
            pool.close()
            pool.join()
            for result in results:
                result_mc_class_counts = result.get()
                for mc_class in result_mc_class_counts:
                    mc_class_counts[mc_class]+=result_mc_class_counts[mc_class]
            if bh == False and isinstance(control,str):
                best_pvalues = calculate_control_fdr([path_to_files+filename+"_binom_results.tsv" for filename in allc_files],control_pvalues,sig_cutoff,mc_class_counts)
            else:
                if isinstance(control,float):
                    print "Must use benjamini-hochberg correction if no control genome is provided"
                best_pvalues = benjamini_hochberg_correction_call_methylated_sites([path_to_files+filename+"_binom_results.tsv" for filename in allc_files],mc_class_counts,sig_cutoff)
            pool=multiprocessing.Pool(num_procs)
            #Normally I can't do this because all the chromosomes are mixed together, but in this case it's fine.       
            for allc_file in [path_to_files+filename+"_binom_results.tsv" for filename in allc_files]:
                pool.apply_async(filter_files_by_pvalue,([allc_file],path_to_files+sample,best_pvalues,1),{"remove_file":True, "sort_mem":sort_mem})
            pool.close()
            pool.join()
        
        else:
            for allc_file in allc_files:
                mc_class_counts = allc_run_binom_tests(allc_file,non_conversion,min_cov,sort_mem=sort_mem)
            if bh == False and isinstance(control,str):
                best_pvalues = calculate_control_fdr([path_to_files+filename+"_binom_results.tsv" for filename in allc_files],control_pvalues,sig_cutoff,mc_class_counts)
            else:
                if isinstance(control,float):
                    print "Must use benjamini-hochberg correction if no control genome is provided"
                best_pvalues = benjamini_hochberg_correction_call_methylated_sites([path_to_files+filename+"_binom_results.tsv" for filename in allc_files],mc_class_counts,sig_cutoff)
            filter_files_by_pvalue([path_to_files+filename+"_binom_results.tsv" for filename in allc_files],path_to_files+sample,best_pvalues,num_procs,remove_file=True,sort_mem=sort_mem)
    else:
        try:
            f = open(path_to_files+sample+"_mpileup_output.tsv",'r')
            f.close()
            print(sample+"_mpileup_output.tsv exists. Using this file for calling methylated cytosines.")
        except IOError:
            print_checkpoint(sample+"_mpileup_output.tsv does not exist. Beginning mpileup.")
            with open(path_to_files+sample+"_mpileup_output.tsv",'w') as f:
                subprocess.check_call(shlex.split(path_to_samtools+"samtools mpileup -Q "+str(min_base_quality)+" -B "+mpileup_quality+" -f "+reference+" "+path_to_files+inputf),stdout=f)
        print_checkpoint("Begin calling methylated cytosines")
        if binom_test == True:
            if num_procs > 1:
                results = []
                try:
                    subprocess.check_call(shlex.split("rm -f "+path_to_files+sample+"_mpileup_output_chunk_[0-9]"))
                except:
                    pass
                try:
                    subprocess.check_call(shlex.split("rm -f "+path_to_files+sample+"_mpileup_output_chunk_[0-9][0-9]"))
                except:
                    pass
                try:
                    subprocess.check_call(shlex.split("rm -f "+path_to_files+sample+"_mpileup_output_chunk_[0-9][0-9][0-9]"))
                except:
                    pass
                split_mpileup_file(num_procs,path_to_files+sample+"_mpileup_output.tsv",path_to_files+sample+"_mpileup_output_chunk_")
                print_checkpoint("Begin binomial tests")
                pool=multiprocessing.Pool(num_procs)
                chunks = glob.glob(path_to_files+sample+"_mpileup_output_chunk_*")
                for chunk in chunks:
                    results.append(pool.apply_async(run_binom_tests,(chunk,non_conversion,min_cov),{"sort_mem":sort_mem}))
                pool.close()
                pool.join()
                pvalues = []
                for result in results:
                    result_mc_class_counts = result.get()
                    for mc_class in result_mc_class_counts:
                        mc_class_counts[mc_class]+=result_mc_class_counts[mc_class]
                #files = [path_to_files+i+"_binom_results.tsv" for i in chunks]
                files = [i+"_binom_results.tsv" for i in chunks]
                subprocess.check_call(shlex.split("rm "+" ".join(chunks)))
            else:
                print_checkpoint("Begin binomial tests")
                mc_class_counts = run_binom_tests(path_to_files+sample+"_mpileup_output.tsv",non_conversion,min_cov=min_cov,sort_mem=sort_mem)
                files = [path_to_files+sample+"_mpileup_output.tsv_binom_results.tsv"]
            print_checkpoint("Begin adjusting p-values")
            if bh == False and isinstance(control,str):
                best_pvalues = calculate_control_fdr(files,control_pvalues,sig_cutoff,mc_class_counts)
            else:
                if isinstance(control,float):
                    print "Must use benjamini-hochberg correction if no control genome is provided"
                best_pvalues = benjamini_hochberg_correction_call_methylated_sites(files,mc_class_counts,sig_cutoff)
            filter_files_by_pvalue(files,path_to_files+sample,best_pvalues,num_procs,remove_file=True,sort_mem=sort_mem)

        elif min_mc !=0:
            run_mc_filter(path_to_files+sample+"_mpileup_output.tsv",min_cov=min_cov,output=sample,min_mc=min_mc)
        else:
            run_pvalue_lookup_filter(path_to_files+sample+"_mpileup_output.tsv",pvalue_lookup,min_pvalue,control_sites,min_cov=min_cov,output=sample)
    print_checkpoint("Done")
    
def filter_files_by_pvalue(files,output,best_pvalues,num_procs,remove_file=True, sort_mem="500M"):
    """
    sort_mem is the parameter to pass to unix sort with -S/--buffer-size command
    """
    if sort_mem:
        if sort_mem.find("-S") == -1:
            sort_mem = " -S " + sort_mem
    else:
        sort_mem = ""
    output_files = {}
    for filen in files:
        f = open(filen,'r')
        for line in f:
            line = line.rstrip()
            fields = line.split("\t")
            if fields[0] not in output_files:
                output_files[fields[0]] = open(output[:output.rfind("/")+1]+"allc_"+output[output.rfind("/")+1:]+"_"+fields[0]+".tsv",'w')
                output_files[fields[0]].write("\t".join(["chr","pos","strand","mc_class","mc_count","total","methylated"])+"\n")
            if fields[6] != "2.0" and float(fields[6]) <= best_pvalues[fields[3]]:
                output_files[fields[0]].write("\t".join(fields[:6])+"\t1\n")
            else:
                output_files[fields[0]].write("\t".join(fields[:6])+"\t0\n")
        f.close()
        if remove_file == True:
            subprocess.check_call(shlex.split("rm "+filen))    
    print_checkpoint("Begin sorting files by position")
    if num_procs > 1:
        pool=multiprocessing.Pool(num_procs)
        for chrom in output_files:
            output_files[chrom].close()
            pool.apply_async(subprocess.check_call, (shlex.split("sort" + sort_mem + " -k 2n,2n -o "+output_files[chrom].name+" "+output_files[chrom].name),))
        pool.close()
        pool.join()
    else:
        for chrom in output_files:
            output_files[chrom].close()
            subprocess.check_call(shlex.split("sort" + sort_mem + " -k 2n,2n -o "+output_files[chrom].name+" "+output_files[chrom].name))
def calculate_control_fdr(files,control_pvalues,fdr_target,mc_class_counts):
    pvalues = {}
    file_handles = {}
    lines = {}
    fields = {}
    for filen in files:
        file_handles[filen] = open(filen,'r')
        lines[filen] = file_handles[filen].readline().rstrip()
        fields[filen] = lines[filen].split("\t")
        try:
            pvalues[filen] = float(fields[filen][6])
        except:
            pvalues[filen]= 2.0

    
    control_index = 0
    control_count = 0
    pvalue_count = {}
    best_fdr = {}
    best_pvalue = {}
    control_pvalues = sorted(control_pvalues)
    control_quantiles = mquantiles(control_pvalues,prob = [i / 1000000.0 for i in range(0,1000000)],alphap=1,betap=1)
    #I need this because I carry out way more actual tests than control tests
    #the number of control pvalues I see at a certain pvalue is based on the 
    #size of the control region, so the number I'd expect in the real dataset
    #needs to be scaled up.
    scale_factor = {}
    for first in ["A","T","C","G","N"]:
        for second in ["A","T","C","G","N"]:
            scale_factor["C"+first+second]= float(mc_class_counts["C"+first+second]) / len(control_pvalues)
    while [i for i in pvalues if pvalues[i]!=2.0]:
        min_key = min(pvalues,key=pvalues.get)
        pvalue = pvalues[min_key]
        mc_class = fields[min_key][3]
        if mc_class not in pvalue_count:
            pvalue_count[mc_class] = 0
            best_fdr[mc_class] = 0
            best_pvalue[mc_class] = 0
        pvalue_count[mc_class] += len([i for filen in pvalues if pvalues[filen]==pvalue and fields[filen][3]==mc_class])
        while control_index < len(control_pvalues) and control_pvalues[control_index]<= pvalue:
            control_count += 1
            control_index += 1
        fdr = len(control_pvalues) * scale_factor[mc_class] * (bisect.bisect_left(control_quantiles,pvalue) / 1000000.0) / float(pvalue_count[mc_class])
        #fdr = (control_count) * scale_factor[mc_class] / float(pvalue_count[mc_class])
        #if (fdr <= fdr_target and fdr >= best_fdr[mc_class]) or best_fdr[mc_class] == 0 or (fdr > fdr_target and fdr <= best_fdr[mc_class]):
        if (fdr <= fdr_target and fdr >= best_fdr[mc_class]):

            best_fdr[mc_class] = fdr
            best_pvalue[mc_class] = pvalue
        for filen in pvalues:
            if pvalues[filen] == pvalue and fields[filen][3]==mc_class:
                lines[filen] = file_handles[filen].readline().rstrip()
                fields[filen] = lines[filen].split("\t")
                try:
                    pvalues[filen] = float(fields[filen][6])
                except:
                    pvalues[filen]= 2.0            
    for mc_class in best_pvalue:
        if mc_class_counts[mc_class] > 0:
            print "The closest p-value cutoff for "+mc_class+" at your desired FDR is "+str(best_pvalue[mc_class])+" which corresponds to an FDR of "+str(best_fdr[mc_class])
    for filen in file_handles:
        file_handles[filen].close()
    return best_pvalue
def get_pvalue_dist(coverage,control_sites):
    """
    This function will calculate the p-values for a given coverage based on an input control genome
    
    coverage is an integer specifying a particular coverage
    
    control_sites is a list of sites from the control genomes. Each site consists of two integers
        the number of methylated reads at the at site and the total number of reads. Consequently,
        this argument is a list of lists like so: [[1,3],[2,5],...]
    """
    dist={}
    total_control_sites=len(control_sites)
    combinations = {}
    for mc_count in range(0,coverage+1):
        dist[mc_count]=0
        if coverage - mc_count in combinations:
            combinations[mc_count] = combinations[coverage-mc_count]
        else:
            combinations[mc_count] = choose(coverage,mc_count)
    for site in control_sites:
        mc_level = site[0] / float(site[1])
        for mc_count in range(0,coverage+1):   
            dist[mc_count]+=(combinations[mc_count])*(mc_level**mc_count)*((1-mc_level)**(coverage-mc_count))
    dist[coverage]=dist[coverage]/total_control_sites 
    for mc_count in range(0,coverage)[::-1]:
        dist[mc_count]=dist[mc_count]/total_control_sites
        dist[mc_count]+=dist[mc_count+1]
    return dist

def choose(n, k):
    """
    Taken from:
    http://stackoverflow.com/questions/3025162/statistics-combinations-in-python
    
    A fast way to calculate binomial coefficients by Andrew Dalke (contrib).
    """
    if 0 <= k <= n:
        ntok = 1
        ktok = 1
        for t in xrange(1, min(k, n - k) + 1):
            ntok *= n
            ktok *= t
            n -= 1
        return ntok // ktok
    else:
        return 0

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
def run_mc_filter(filen,min_cov=1,output="temp",min_mc=0):
    """
    This function is used by call_methylated_sites to decide which sites are methylated.
    file is a string containing the path to an mpileup file
    min_cov is the minimum number of reads a site must have to be tested
    output is a string indicating the name of the output file like so:
        allc_<output>_<chr>.tsv
    min_mc is the minimum number of mCs that must be observed to 
    """
    output_files={}
    reverse_complement = {"A":"T","C":"G","G":"C","T":"A","N":"N"}
    f = open(filen,'r')
    line1 = f.readline().rstrip()
    line2 = f.readline().rstrip()
    line3 = f.readline().rstrip()
    line4 = f.readline().rstrip()
    line5 = f.readline().rstrip()
    while line5:
        fields1 = line1.split("\t")
        fields2 = line2.split("\t")
        fields3 = line3.split("\t")
        fields4 = line4.split("\t")
        fields5 = line5.split("\t")
        
        fields1[2] = fields1[2].upper()
        fields2[2] = fields2[2].upper()
        fields3[2] = fields3[2].upper()
        fields4[2] = fields4[2].upper()
        fields5[2] = fields5[2].upper()
        
        if fields3[2] == "C":
            #make sure bases are contiguous
            if fields4[0] == fields3[0] and int(fields4[1]) == int(fields3[1]) + 1:
                if fields4[2] not in ["A","C","G","T","N"]:
                    mc_class = "CN"
                else:
                    mc_class = "C"+fields4[2]
            else:
                mc_class = "CN"
            if fields5[0] == fields3[0] and int(fields5[1]) == int(fields3[1]) + 2:
                if fields5[2] not in ["A","C","G","T","N"]:
                    mc_class += "N"
                else:
                    mc_class += fields5[2]
            else:
                mc_class += "N"
            chrom = fields3[0].replace("chr","")
            unconverted_c = fields3[4].count(".")
            converted_c = fields3[4].count("T")
            total = unconverted_c+ converted_c
            
            if chrom not in output_files:
                output_files[chrom] = open("allc_"+output+"_"+chrom+".tsv",'w')
                output_files[chrom].write("\t".join(["chr","pos","strand","mc_class","mc_count","total","methylated"])+"\n")

            if total>= min_cov and unconverted_c >=min_mc:
                output_files[chrom].write("\t".join(map(str,[chrom,fields3[1],"+",mc_class,unconverted_c,total,"1"]))+"\n")
            elif total!=0:
                output_files[chrom].write("\t".join(map(str,[chrom,fields3[1],"+",mc_class,unconverted_c,total,"0"]))+"\n")

                
        elif fields3[2] == "G":
            #make sure bases are contiguous
            if fields2[0] == fields3[0] and int(fields2[1]) == int(fields3[1]) - 1:
                try:
                    mc_class = "C"+reverse_complement[fields2[2]]
                except:
                    mc_class = "CN"
            else:
                mc_class = "CN"
            if fields1[0] == fields3[0] and int(fields1[1]) == int(fields3[1]) - 2:
                try:
                    mc_class += reverse_complement[fields1[2]]
                except:
                    mc_class += "N"
            else:
                mc_class += "N"
            chrom = fields3[0].replace("chr","")
            unconverted_c = fields3[4].count(",")
            converted_c = fields3[4].count("a") 
            total = unconverted_c + converted_c
        
            if chrom not in output_files:
                output_files[chrom] = open("allc_"+output+"_"+chrom+".tsv",'w')
                output_files[chrom].write("\t".join(["chr","pos","strand","mc_class","mc_count","total","methylated"])+"\n")

            if total>= min_cov and unconverted_c >=min_mc:
                output_files[chrom].write("\t".join(map(str,[chrom,fields3[1],"-",mc_class,unconverted_c,total,"1"]))+"\n")
            elif total!=0:
                output_files[chrom].write("\t".join(map(str,[chrom,fields3[1],"-",mc_class,unconverted_c,total,"0"]))+"\n")

        line1 = line2
        line2 = line3
        line3 = line4
        line4 = line5
        line5 = f.readline().rstrip()
def run_pvalue_lookup_filter(filen,pvalue_lookup,min_pvalue,control_sites,output="temp",min_cov=1):
    """
    This function is used by call_methylated_sites to parallelize the binomial test.
    file is a string containing the path to an mpileup file
    non_conversion is a float indicating the estimated non-conversion rate
    min_cov is the minimum number of reads a site must have to be tested
    """
    reverse_complement = {"A":"T","C":"G","G":"C","T":"A","N":"N"}
    f = open(filen,'r')
    output_files = {}
    line1 = f.readline().rstrip()
    line2 = f.readline().rstrip()
    line3 = f.readline().rstrip()
    line4 = f.readline().rstrip()
    line5 = f.readline().rstrip()
    
    fields1 = line1.split("\t")
    fields2 = line2.split("\t")
    fields3 = line3.split("\t")
    fields4 = line4.split("\t")
    fields5 = line5.split("\t")
    
    fields1[2] = fields1[2].upper()
    fields2[2] = fields2[2].upper()
    fields3[2] = fields3[2].upper()
    fields4[2] = fields4[2].upper()
    while line5:
        fields5[2] = fields5[2].upper()
        
        if fields3[2] == "C":
            #make sure bases are contiguous
            if fields4[0] == fields3[0] and int(fields4[1]) == int(fields3[1]) + 1:
                if fields4[2] not in ["A","C","G","T","N"]:
                    mc_class = "CN"
                else:
                    mc_class = "C"+fields4[2]
            else:
                mc_class = "CN"
            if fields5[0] == fields3[0] and int(fields5[1]) == int(fields3[1]) + 2:
                if fields5[2] not in ["A","C","G","T","N"]:
                    mc_class += "N"
                else:
                    mc_class += fields5[2]
            else:
                mc_class += "N"
            chrom = fields3[0].replace("chr","")
            unconverted_c = fields3[4].count(".")
            converted_c = fields3[4].count("T")
            total = unconverted_c+ converted_c
            if total >= min_cov:
                if total not in pvalue_lookup:
                    pvalue_lookup[total] = get_pvalue_dist(total,control_sites)
                p_value = pvalue_lookup[total][unconverted_c]
                if chrom not in output_files:
                    output_files[chrom] = open("allc_"+output+"_"+chrom+".tsv",'w')
                mc_indicator = str(int(p_value<=min_pvalue))
                output_files[chrom].write("\t".join(map(str,[chrom,fields3[1],"+",mc_class,unconverted_c,total,mc_indicator]))+"\n")
            elif total != 0:
                if chrom not in output_files:
                    output_files[chrom] = open("allc_"+output+"_"+chrom+".tsv",'w')
                mc_indicator = "0"
                output_files[chrom].write("\t".join(map(str,[chrom,fields3[1],"+",mc_class,unconverted_c,total,mc_indicator]))+"\n")

        elif fields3[2] == "G":
            #make sure bases are contiguous
            if fields2[0] == fields3[0] and int(fields2[1]) == int(fields3[1]) - 1:
                try:
                    mc_class = "C"+reverse_complement[fields2[2]]
                except:
                    mc_class = "CN"
            else:
                mc_class = "CN"
            if fields1[0] == fields3[0] and int(fields1[1]) == int(fields3[1]) - 2:
                try:
                    mc_class += reverse_complement[fields1[2]]
                except:
                    mc_class += "N"
            else:
                mc_class += "N"
            chrom = fields3[0].replace("chr","")
            unconverted_c = fields3[4].count(",")
            converted_c = fields3[4].count("a") 
            total = unconverted_c + converted_c
            if total >= min_cov:
                if total not in pvalue_lookup:
                    pvalue_lookup[total] = get_pvalue_dist(total,control_sites)
                p_value = pvalue_lookup[total][unconverted_c]
                if chrom not in output_files:
                    output_files[chrom] = open("allc_"+output+"_"+chrom+".tsv",'w')
                mc_indicator = str(int(p_value<=min_pvalue))
                output_files[chrom].write("\t".join(map(str,[chrom,fields3[1],"-",mc_class,unconverted_c,total,mc_indicator]))+"\n")
            elif total != 0:
                if chrom not in output_files:
                    output_files[chrom] = open("allc_"+output+"_"+chrom+".tsv",'w')
                mc_indicator = "0"
                output_files[chrom].write("\t".join(map(str,[chrom,fields3[1],"-",mc_class,unconverted_c,total,mc_indicator]))+"\n")
        line1 = line2
        line2 = line3
        line3 = line4
        line4 = line5
        line5 = f.readline().rstrip()
        
        fields1 = fields2
        fields2 = fields3
        fields3 = fields4
        fields4 = fields5
        fields5 = line5.split("\t")
        
    for chrom in output_files:
        output_files[chrom].close()
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
def run_binom_tests(filen,non_conversion,min_cov=1,sort_mem="500M"):
    """
    This function is used by call_methylated_sites to parallelize the binomial test.
    file is a string containing the path to an mpileup file
    
    non_conversion is a float indicating the estimated non-conversion rate and sequencing error
        
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
    line1 = f.readline().rstrip('\n')
    line2 = f.readline().rstrip('\n')
    line3 = f.readline().rstrip('\n')
    line4 = f.readline().rstrip('\n')
    line5 = f.readline().rstrip('\n')
    g = open(filen+"_binom_results.tsv",'w')
    
    #Deal with edge case of the first two positions
    fields1 = line1.split("\t")
    fields2 = line2.split("\t")
    fields3 = line3.split("\t")
    fields4 = line4.split("\t")
    fields5 = line5.split("\t")
    #For the first position
    if fields1[2] == "C" and fields1[3]!="0":
        #make sure bases are contiguous
        if fields2[0] == fields1[0] and int(fields2[1]) == int(fields1[1]) + 1:
            if fields2[2] not in ["A","C","G","T","N"]:
                mc_class = "CN"
            else:
                mc_class = "C"+fields2[2]
        else:
            mc_class = "CN"
        if fields3[0] == fields1[0] and int(fields3[1]) == int(fields1[1]) + 2:
            if fields3[2] not in ["A","C","G","T","N"]:
                mc_class += "N"
            else:
                mc_class += fields3[2]
        else:
            mc_class += "N"
        chrom = fields1[0].replace("chr","")    
        unconverted_c = fields1[4].count(".")
        converted_c = fields1[4].count("T")
        total = unconverted_c+ converted_c
        if total >= min_cov and unconverted_c != 0:
            try:
                p_value = obs_pvalues[(unconverted_c,total)]
            except:
                p_value = sci.binom.sf(unconverted_c-1,total,non_conversion)
                obs_pvalues[(unconverted_c,total)] = p_value
            g.write("\t".join(map(str,[chrom,fields1[1],"+",mc_class,unconverted_c,total,p_value]))+"\n")
            mc_class_counts[mc_class]+=1
        elif total != 0:
            #a dummy value that will always sort to the bottom of the BH correction and be interpreted as
            #a unmethylated site
            p_value = 2.0
            g.write("\t".join(map(str,[chrom,fields1[1],"+",mc_class,unconverted_c,total,p_value]))+"\n")
    elif fields1[2] == "G" and fields1[3]!="0":
        #position is right on the edge so has to be CNN
        mc_class = "CNN"
        chrom = fields1[0].replace("chr","")
        unconverted_c = fields1[4].count(",")
        converted_c = fields1[4].count("a") 
        total = unconverted_c + converted_c
        if total >= min_cov and unconverted_c != 0:
            try:
                p_value = obs_pvalues[(unconverted_c,total)]
            except:
                p_value = sci.binom.sf(unconverted_c-1,total,non_conversion)
                obs_pvalues[(unconverted_c,total)] = p_value
            g.write("\t".join(map(str,[chrom,fields1[1],"-",mc_class,unconverted_c,total,p_value]))+"\n")
            mc_class_counts[mc_class]+=1
        elif total != 0:
            #a dummy value that will always sort to the bottom of the BH correction and be interpreted as
            #a unmethylated site
            p_value = 2.0
            g.write("\t".join(map(str,[chrom,fields1[1],"-",mc_class,unconverted_c,total,p_value]))+"\n")        
    #FOR THE SECOND POSITION
    if fields2[2] == "C" and fields2[3]!="0":
        #make sure bases are contiguous
        if fields3[0] == fields2[0] and int(fields3[1]) == int(fields2[1]) + 1:
            if fields3[2] not in ["A","C","G","T","N"]:
                mc_class = "CN"
            else:
                mc_class = "C"+fields3[2]
        else:
            mc_class = "CN"
        if fields4[0] == fields2[0] and int(fields4[1]) == int(fields2[1]) + 2:
            if fields4[2] not in ["A","C","G","T","N"]:
                mc_class += "N"
            else:
                mc_class += fields4[2]
        else:
            mc_class += "N"
        chrom = fields2[0].replace("chr","")    
        unconverted_c = fields2[4].count(".")
        converted_c = fields2[4].count("T")
        total = unconverted_c+ converted_c
        if total >= min_cov and unconverted_c != 0:
            try:
                p_value = obs_pvalues[(unconverted_c,total)]
            except:
                p_value = sci.binom.sf(unconverted_c-1,total,non_conversion)
                obs_pvalues[(unconverted_c,total)] = p_value
            g.write("\t".join(map(str,[chrom,fields2[1],"+",mc_class,unconverted_c,total,p_value]))+"\n")
            mc_class_counts[mc_class]+=1
        elif total != 0:
            #a dummy value that will always sort to the bottom of the BH correction and be interpreted as
            #a unmethylated site
            p_value = 2.0
            g.write("\t".join(map(str,[chrom,fields2[1],"+",mc_class,unconverted_c,total,p_value]))+"\n")
    elif fields2[2] == "G" and fields2[3]!="0":
        #make sure bases are contiguous
        if fields1[0] == fields2[0] and int(fields1[1]) == int(fields2[1]) - 1:
            try:
                mc_class = "C"+reverse_complement[fields1[2]]
            except:
                mc_class = "CN"
        else:
            mc_class = "CN"
        #This position is at the edge so it'll always be C[ACTG]N
        mc_class += "N"
        chrom = fields2[0].replace("chr","")
        unconverted_c = fields2[4].count(",")
        converted_c = fields2[4].count("a") 
        total = unconverted_c + converted_c
        if total >= min_cov and unconverted_c != 0:
            try:
                p_value = obs_pvalues[(unconverted_c,total)]
            except:
                p_value = sci.binom.sf(unconverted_c-1,total,non_conversion)
                obs_pvalues[(unconverted_c,total)] = p_value
            g.write("\t".join(map(str,[chrom,fields2[1],"-",mc_class,unconverted_c,total,p_value]))+"\n")
            mc_class_counts[mc_class]+=1
        elif total != 0:
            #a dummy value that will always sort to the bottom of the BH correction and be interpreted as
            #a unmethylated site
            p_value = 2.0
            g.write("\t".join(map(str,[chrom,fields2[1],"-",mc_class,unconverted_c,total,p_value]))+"\n")
    
    while line5:
        fields1 = line1.split("\t")
        fields2 = line2.split("\t")
        fields3 = line3.split("\t")
        fields4 = line4.split("\t")
        fields5 = line5.split("\t")
        
        fields1[2] = fields1[2].upper()
        fields2[2] = fields2[2].upper()
        fields3[2] = fields3[2].upper()
        fields4[2] = fields4[2].upper()
        fields5[2] = fields5[2].upper()
        
        if fields3[2] == "C" and fields3[3]!="0":
            #make sure bases are contiguous
            if fields4[0] == fields3[0] and int(fields4[1]) == int(fields3[1]) + 1:
                if fields4[2] not in ["A","C","G","T","N"]:
                    mc_class = "CN"
                else:
                    mc_class = "C"+fields4[2]
            else:
                mc_class = "CN"
            if fields5[0] == fields3[0] and int(fields5[1]) == int(fields3[1]) + 2:
                if fields5[2] not in ["A","C","G","T","N"]:
                    mc_class += "N"
                else:
                    mc_class += fields5[2]
            else:
                mc_class += "N"
            chrom = fields3[0].replace("chr","")    
            unconverted_c = fields3[4].count(".")
            converted_c = fields3[4].count("T")
            total = unconverted_c+ converted_c
            if total >= min_cov and unconverted_c != 0:
                try:
                    p_value = obs_pvalues[(unconverted_c,total)]
                except:
                    p_value = sci.binom.sf(unconverted_c-1,total,non_conversion)
                    obs_pvalues[(unconverted_c,total)] = p_value
                g.write("\t".join(map(str,[chrom,fields3[1],"+",mc_class,unconverted_c,total,p_value]))+"\n")
                mc_class_counts[mc_class]+=1
            elif total != 0:
                #a dummy value that will always sort to the bottom of the BH correction and be interpreted as
                #a unmethylated site
                p_value = 2.0
                g.write("\t".join(map(str,[chrom,fields3[1],"+",mc_class,unconverted_c,total,p_value]))+"\n")
        elif fields3[2] == "G" and fields3[3]!="0":
            #make sure bases are contiguous
            if fields2[0] == fields3[0] and int(fields2[1]) == int(fields3[1]) - 1:
                try:
                    mc_class = "C"+reverse_complement[fields2[2]]
                except:
                    mc_class = "CN"
            else:
                mc_class = "CN"
            if fields1[0] == fields3[0] and int(fields1[1]) == int(fields3[1]) - 2:
                try:
                    mc_class += reverse_complement[fields1[2]]
                except:
                    mc_class += "N"
            else:
                mc_class += "N"
            chrom = fields3[0].replace("chr","")
            unconverted_c = fields3[4].count(",")
            converted_c = fields3[4].count("a") 
            total = unconverted_c + converted_c
            if total >= min_cov and unconverted_c != 0:
                try:
                    p_value = obs_pvalues[(unconverted_c,total)]
                except:
                    p_value = sci.binom.sf(unconverted_c-1,total,non_conversion)
                    obs_pvalues[(unconverted_c,total)] = p_value
                g.write("\t".join(map(str,[chrom,fields3[1],"-",mc_class,unconverted_c,total,p_value]))+"\n")
                mc_class_counts[mc_class]+=1
            elif total != 0:
                #a dummy value that will always sort to the bottom of the BH correction and be interpreted as
                #a unmethylated site
                p_value = 2.0
                g.write("\t".join(map(str,[chrom,fields3[1],"-",mc_class,unconverted_c,total,p_value]))+"\n")
        line1 = line2
        line2 = line3
        line3 = line4
        line4 = line5
        line5 = f.readline().rstrip('\n')
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
     parser_pipeline.add_argument('--quality_version', type=str, default="1.8", help="either an integer indicating the base offset for \
        the quality scores or a float indicating which version of casava was used to generate the fastq files.")
     parser_pipeline.add_argument('--path_to_samtools', type=str, default="", help='Path to samtools installation (default is current dir)')     
     parser_pipeline.add_argument('--path_to_bowtie', type=str, default="", help='Path to bowtie installation (default is current dir)')
     parser_pipeline.add_argument('--bowtie_options', type=str, nargs='+', help="list of strings indicating options you'd like passed to bowtie \
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
     parser_pipeline.add_argument('--num_reads', type=int, default=-1, help="integer indicating how many reads you'd like to process at a time. \
        Breaking up the processing this way can reduce the size of some of the intermediate files that are created. Use this if you're strapped for space.")
     parser_pipeline.add_argument('--save_space', type=bool, default=True, help="indicates whether or not you'd like to perform read collapsing \
        right after mapping or once all the libraries have been mapped. If you wait until after everything has been mapped, the collapsing can be \
        parallelized. Otherwise the collapsing will have to be done serially. The trade-off is that you must keep all the mapped files around, rather \
        than deleting them as they are processed, which can take up a considerable amount of space. It's safest to set this to True.")
     parser_pipeline.add_argument('--bowtie2', type=bool, default=False, help="Specifies whether to use the bowtie2 aligner instead of bowtie")
     parser_pipeline.add_argument('--sort_mem', type=str, help="Parameter to pass to unix sort with -S/--buffer-size command")
     parser_pipeline.add_argument('--path_to_output', type=str, default="", help="Path to a directory where you would like the output to be stored. \
        The default is the same directory as the input fastqs.")
     parser_pipeline.add_argument('--in_mem', type=bool, default=False, help="Does not write files as part of run_mapping, keeps all input/output in memory")
     parser_pipeline.add_argument('--path_to_picard', type=str, default=False, help="The path to MarkDuplicates jar from picard. Default is false indicating that you don't want to use this jar for duplication removal")
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
         if not args.bowtie_options:
             args.bowtie_options = ["-S","-k 1","-m 1","--chunkmbs 3072","--best","--strata","-o 4","-e 80","-l 20","-n 0"]
        
         run_methylation_pipeline(args.files,args.libraries,args.sample,args.forward_ref,args.reverse_ref,args.ref_fasta,
                                  args.unmethylated_control,args.quality_version,args.path_to_samtools,args.path_to_bowtie,
                                  args.bowtie_options,args.num_procs,args.trim_reads,
                                  args.path_to_cutadapt,args.adapter_seq,args.max_adapter_removal,
                                  args.overlap_length,args.zero_cap,args.error_rate,args.min_qual_score,args.min_read_len,
                                  args.sig_cutoff,args.min_cov,args.binom_test,args.keep_temp_files,args.num_reads,
                                  args.save_space,
                                  args.bowtie2,args.sort_mem,args.path_to_output,args.in_mem)
                                  
     elif args.command == "call_methylated_sites":
         call_methylated_sites(args.inputf, args.sample, args.reference, args.control, args.casava_version, args.sig_cutoff,
                              args.num_procs, args.min_cov, args.binom_test, args.min_mc, args.path_to_samtools, args.sort_mem,
                              args.bh, args.path_to_files)
    
if __name__ == '__main__':
    parse_args()
