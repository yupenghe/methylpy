import sys
import multiprocessing
import subprocess
import scipy.stats as sci
from scipy.stats.mstats import mquantiles
from methylpy.utilities import print_checkpoint, print_error
from methylpy.utilities import split_fastq_file
from methylpy.utilities import split_fastq_file_pbat
import pdb
import shlex
import itertools
import re
import glob
import cStringIO as cStr
import bisect
import gzip

def run_methylation_pipeline(read_files, libraries, sample,
                             forward_reference, reverse_reference, reference_fasta,
                             unmethylated_control="chrL:",
                             path_to_output="", sig_cutoff=0.01,
                             num_procs=1, sort_mem="500M",
                             num_upstr_bases=1, num_downstr_bases=2,
                             generate_allc_file=True,split_allc_files=False,
                             generate_mpileup_file=True, compress_output=True,
                             binom_test=True, bh=True, min_cov=2,
                             trim_reads=True, path_to_cutadapt="",
                             pbat=False,
                             bowtie2=False, path_to_aligner="", aligner_options=[],
                             remove_clonal=True,keep_clonal_stats=False,
                             path_to_picard="",java_options="-Xmx20g",
                             path_to_samtools="",
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

    #Default bowtie option
    if len(aligner_options) == 0:
        if not bowtie2:
            aligner_options = ["-S", "-k 1", "-m 1", "--chunkmbs 3072",
                               "--best", "--strata", "-o 4", "-e 80",
                               "-l 20", "-n 0"]

    # CASAVA >= 1.8
    aligner_options.append("--phred33-quals")
    quality_base = 33

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
                                            aligner_options=aligner_options,
                                            pbat=pbat,
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
                                            bowtie2=bowtie2, sort_mem=sort_mem)
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

    print_checkpoint("Begin calling mCs")
    if remove_clonal == True:
        output_bam_file = path_to_output+sample+"_processed_reads_no_clonal.bam"
    else:
        output_bam_file = path_to_output+sample+"_processed_reads.bam"

    if generate_allc_file:
        call_methylated_sites(output_bam_file,
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

def run_mapping(current_library, library_files, sample,
                forward_reference, reverse_reference, reference_fasta,
                path_to_output="",
                path_to_samtools="", path_to_aligner="",
                aligner_options=[],pbat=False,
                num_procs=1, trim_reads=True, path_to_cutadapt="",
                adapter_seq="AGATCGGAAGAGCACACGTCTG",
                max_adapter_removal=None, overlap_length=None, zero_cap=None,
                quality_base=None, error_rate=None,
                min_qual_score=10, min_read_len=30,
                keep_temp_files=False, bowtie2=False,
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

    if len(aligner_options) == 0:
        if not bowtie2:
            aligner_options=["-S", "-k 1", "-m 1", "--chunkmbs 3072",
                             "--best", "--strata", "-o 4", "-e 80", "-l 20", "-n 0"]

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
        quality_trim([file_path+"_split_"+str(i) for i in xrange(0, num_procs)],
                     output=[file_path+"_split_trimmed_"+str(i)
                             for i in xrange(0, num_procs)],
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
                                                          for i in xrange(0,num_procs)])))

        print_checkpoint("Begin converting reads for "+file_name)
        if num_procs > 1:
            pool = multiprocessing.Pool(num_procs)
            for inputf, output in zip([file_path+"_split_trimmed_"+str(i) for i in xrange(0, num_procs)],
                                     [file_path+"_split_trimmed_converted_"+str(i)
                                      for i in xrange(0, num_procs)]):
                pool.apply_async(convert_reads,(inputf,output))
            pool.close()
            pool.join()
        else:
            for inputf, output in zip([file_path+"_split_trimmed_"+str(i) for i in xrange(0, num_procs)],
                                     [file_path+"_split_trimmed_converted_"+str(i)
                                      for i in xrange(0, num_procs)]):
                                          convert_reads(inputf,output)
        subprocess.check_call(shlex.split("rm "+
                                          " ".join([file_path+"_split_trimmed_"+str(i)
                                                    for i in xrange(0,num_procs)])))
        input_fastq = [file_path+"_split_trimmed_converted_"+str(i) for i in xrange(0, num_procs)]
    else:
        print_checkpoint("No trimming on reads")
        print_checkpoint("Begin converting reads for "+file_name)
        if num_procs > 1:
            pool = multiprocessing.Pool(num_procs)
            for inputf, output in zip([file_path+"_split_"+str(i) for i in xrange(0, num_procs)],
                                     [file_path+"_split_converted_"+str(i) for i in xrange(0, num_procs)]):
                pool.apply_async(convert_reads, (inputf, output))
            pool.close()
            pool.join()
        else:
            for inputf, output in zip([file_path+"_split_"+str(i) for i in xrange(0, num_procs)],
                                     [file_path+"_split_converted_"+str(i) for i in xrange(0, num_procs)]):
                convert_reads(inputf, output)
        subprocess.check_call(shlex.split("rm "+" ".join([file_path+"_split_"+str(i)
                                                          for i in xrange(0, num_procs)])))
        input_fastq = [file_path+"_split_converted_"+str(i) for i in xrange(0, num_procs)]

    print_checkpoint("Begin Running Bowtie for "+current_library)
    total_unique = run_bowtie(current_library,
                              input_fastq,
                              sample,
                              forward_reference, reverse_reference, reference_fasta,
                              path_to_output=path_to_output,
                              aligner_options=aligner_options,
                              path_to_aligner=path_to_aligner, num_procs=num_procs,
                              keep_temp_files=keep_temp_files,
                              bowtie2=bowtie2, sort_mem=sort_mem)

    subprocess.check_call(shlex.split("rm " + " ".join(input_fastq)))

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

def encode_converted_positions(seq):
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

def decode_converted_positions(seq,indexes,strand):
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

def run_bowtie(current_library,library_read_files,
               sample,
               forward_reference,reverse_reference,reference_fasta,
               path_to_output="",
               path_to_samtools="",
               path_to_aligner="",aligner_options="",
               num_procs=1,keep_temp_files=False, bowtie2=False, sort_mem="500M"):
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
    if len(aligner_options) == 0:
        if not bowtie2:
            aligner_options=["-S","-k 1","-m 1","--chunkmbs 3072",
                             "--best","--strata","-o 4","-e 80","-l 20","-n 0"]
    options = aligner_options

    if len(path_to_output) !=0:
        path_to_output+="/"

    prefix = path_to_output+sample+"_"+str(current_library)

    if sort_mem:
        if sort_mem.find("-S") == -1:
            sort_mem = " -S " + sort_mem
    else:
        sort_mem = ""
    if " ".join(options).find(" -p ") == -1:
        options.append("-p "+str(num_procs))
    
    if bowtie2:
        args = [path_to_aligner+"bowtie2"]
        args.extend(options)
        args.append("--norc")
        args.append("-x "+forward_reference)
        args.append("-U "+",".join(library_read_files))
        args.append("-S "+prefix+"_forward_strand_hits.sam")
    else:
        args = [path_to_aligner+"bowtie"]
        args.extend(options)
        args.append("--norc")
        args.append(forward_reference)
        args.append(",".join(library_read_files))
        args.append(prefix+"_forward_strand_hits.sam")
    subprocess.check_call(shlex.split(" ".join(args)))
    print_checkpoint("Processing forward strand hits")
    find_multi_mappers(prefix+"_forward_strand_hits.sam",prefix,num_procs=num_procs,keep_temp_files=keep_temp_files)

    if bowtie2:
        args = [path_to_aligner+"bowtie2"]
        args.extend(options)
        args.append("--nofw")
        args.append("-x "+reverse_reference)
        args.append("-U "+",".join(library_read_files))
        args.append("-S "+prefix+"_reverse_strand_hits.sam")
    else:
        args = [path_to_aligner+"bowtie"]
        args.extend(options)
        args.append("--nofw")
        args.append(reverse_reference)
        args.append(",".join(library_read_files))
        args.append(prefix+"_reverse_strand_hits.sam")
    subprocess.check_call(shlex.split(" ".join(args)))    
    print_checkpoint("Processing reverse strand hits")
    sam_header = find_multi_mappers(prefix+"_reverse_strand_hits.sam",prefix,num_procs=num_procs,append=True,keep_temp_files=keep_temp_files)
    ## Clear temporary files
    if num_procs > 1:
        pool = multiprocessing.Pool(num_procs)
        for file_num in xrange(0,num_procs):
            pool.apply_async(subprocess.check_call,(shlex.split("env LC_COLLATE=C sort" + sort_mem + " -t '\t' -k 1 -o "+prefix+"_sorted_"+str(file_num)+" "+prefix+"_sorted_"+str(file_num)),))
        pool.close()
        pool.join()
    else:
        for file_num in xrange(0,num_procs):
            subprocess.check_call(shlex.split("env LC_COLLATE=C sort" + sort_mem + " -t '\t' -k 1 -o "+prefix+"_sorted_"+str(file_num)+" "+prefix+"_sorted_"+str(file_num)))

    print_checkpoint("Finding multimappers")

    total_unique = merge_sorted_multimap(current_library,
                                         [prefix+"_sorted_"+str(file_num) for file_num in xrange(0,num_procs)],
                                         prefix,
                                         reference_fasta,
                                         path_to_samtools="")
    subprocess.check_call(shlex.split("rm "+" ".join([prefix+"_sorted_"+str(file_num) for file_num in xrange(0,num_procs)])))
    return total_unique
 
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
    subprocess.check_call(shlex.split(path_to_samtools+"samtools sort "+output_bam_file+
                                      " -o "+output_bam_file))
    
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
            total_clonal = fields[5]
        else:
            total_clonal = fields[4]
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
        header = header.next()[1:].strip().replace("chr","")
        if header != query_chrom:
            continue
        # join all sequence lines to one.
        seq = "".join(s.strip() for s in faiter.next())
        return seq

def call_methylated_sites(inputf, sample, reference_fasta, control,sig_cutoff=.01,num_procs = 1,
                          num_upstr_bases=0,num_downstr_bases=2,
                          generate_mpileup_file=True,
                          compress_output=True,
                          split_allc_files=False,
                          buffer_line_number = 100000,
                          min_cov=1,binom_test=True,min_mc=0,path_to_samtools="",
                          sort_mem="500M",bh=True,path_to_files="",min_base_quality=1):

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
        open(path_to_files+inputf+".bai",'r')
    except:
        print_checkpoint("Input not indexed. Indexing...")
        subprocess.check_call(shlex.split(path_to_samtools+"samtools index "+path_to_files+inputf))

    ## Input
    if not generate_mpileup_file:
        cmd = path_to_samtools+"samtools mpileup -Q "+str(min_base_quality)+" -B -f "+reference_fasta+" "+path_to_files+inputf
        pipes = subprocess.Popen(shlex.split(cmd),
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE,
                                 universal_newlines=True)
        fhandle = pipes.stdout
    else:
        with open(path_to_files+sample+"_mpileup_output.tsv",'w') as f:
            subprocess.check_call(
                shlex.split(
                    path_to_samtools+"samtools mpileup -Q "+
                    str(min_base_quality)+" -B -f "+reference_fasta+" "+path_to_files+inputf),
                stdout=f)
        fhandle = open(path_to_files+sample+"_mpileup_output.tsv" ,'r')

    ## Output
    if compress_output:
        output_filehandler = gzip.open(path_to_files+"allc_"+sample+".tsv.gz",'w')
    else:
        output_filehandler = open(path_to_files+"allc_"+sample+".tsv",'w')
        
    complement = {"A":"T","C":"G","G":"C","T":"A","N":"N"}
    cur_chrom = ""
    cur_chrom_nochr = ""
    line_counts = 0
    out = ""
    for line in fhandle:
        fields = line.split("\t")
        if fields[0] != cur_chrom:
            cur_chrom = fields[0]
            cur_chrom_nochr = cur_chrom.replace("chr","")
            seq = fasta_iter(reference_fasta,cur_chrom_nochr)
            if seq != None:
                seq = seq.upper()

        if seq == None:
            continue
    
        if fields[2] == "C":
            pos = int(fields[1])-1
            context = seq[(pos-num_upstr_bases):(pos+num_downstr_bases+1)]
            unconverted_c = fields[4].count(".")
            converted_c = fields[4].count("T")
            cov = unconverted_c+converted_c
            if cov > 0:
                line_counts += 1
                #output_filehandler.write("\t".join([cur_chrom_nochr,str(pos+1),"+",context,
                #                                    str(unconverted_c),str(cov),"1"])+"\n")
                out += "\t".join([cur_chrom_nochr,str(pos+1),"+",context,
                                  str(unconverted_c),str(cov),"1"])+"\n"
        elif fields[2] == "G":
            pos = int(fields[1])-1
            context = "".join([complement[base]
                               for base in reversed(
                                       seq[(pos-num_downstr_bases):(pos+num_upstr_bases+1)]
                               )]
            )
            unconverted_c = fields[4].count(",")
            converted_c = fields[4].count("a")
            cov = unconverted_c+converted_c
            if cov > 0:
                line_counts += 1
                #output_filehandler.write("\t".join([cur_chrom_nochr,str(pos+1),"-",context,
                #                                    str(unconverted_c),str(cov),"1"])+"\n")
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

    if generate_mpileup_file:
        subprocess.check_call(shlex.split("rm -f "+path_to_files+sample+"_mpileup_output.tsv"))

    
    if not split_allc_files:
        return 0
    
    # Split allc files
    if compress_output:
        fhandle = gzip.open(path_to_files+"allc_"+sample+".tsv.gz",'r')
    else:
        fhandle = open(path_to_files+"allc_"+sample+".tsv",'r')
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
                output_handle = gzip.open(path_to_files+"allc_"+sample+"_"+cur_chrom+".tsv.gz",'w')
            else:
                output_handle = open(path_to_files+"allc_"+sample+"_"+cur_chrom+".tsv",'w')
        else:
            # data
            line_counts += 1
            out += line

    if line_counts > 0:
        output_handle.write(out)
        line_counts = 0
        out = ""
    output_handle.close()

    if compress_output:
        subprocess.check_call(shlex.split("rm -f "+path_to_files+"allc_"+sample+".tsv.gz",'w'))
    else:
        subprocess.check_call(shlex.split("rm -f "+path_to_files+"allc_"+sample+".tsv",'w'))

    return 0

def bam_quality_mch_filter(inputf,
                           outputf,
                           reference_fasta,
                           quality_cutoff = 30,
                           min_ch = 3,
                           max_mch_ratio = 0.7,
                           buffer_line_number = 100000,
                           path_to_samtools = ""):
    """
    
    """

    min_ch = int(min_ch)
    max_mch_ratio = float(max_mch_ratio)
    
    # quality filter
    cmd = path_to_samtools+"samtools view -q"+str(quality_cutoff)+" "+inputf
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
            cur_chrom_nochr = cur_chrom.replace("chr","")
            seq = fasta_iter(reference_fasta,cur_chrom_nochr)
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
        else: # - strand
            for ind,base in enumerate(fields[9]):
                pos = ind + read_pos
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

        # apply filter
        tot_ch = float(mch+uch)
        if tot_ch >= min_ch and float(mch)/float(tot_ch) > max_mch_ratio:
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
