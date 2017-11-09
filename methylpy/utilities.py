import sys
import os
import pdb
import math
import shlex
import numpy as np
import subprocess
import time
import operator
import itertools
import multiprocessing
import bz2
import gzip
import collections
from pkg_resources import parse_version


def get_executable_version(exec_name):
    try:
        out = subprocess.check_output(shlex.split(exec_name+" --version"))
    except:
        print_error("Failed to run %s\n" %(exec_name))
    first_line = out.split("\n")[0]
    fields = first_line.split(" ")
    return(fields[-1])

def check_call_mc_dependencies(path_to_samtools="",
                               trim_reads=True,
                               path_to_cutadapt="",
                               bowtie2=True,
                               path_to_aligner="",
                               remove_clonal=True,
                               path_to_picard=""):
    
    # check samtools version
    if len(path_to_samtools) != 0:
        path_to_samtools += "/"
    # check picard
    samtools_version = get_executable_version(path_to_samtools+"samtools")
    if parse_version(samtools_version) < parse_version("1.2"):
        print_error("samtools version %s found.\nmethylpy need at least samtools 1.3\nExit!\n"
                    %(samtools_version) )

    # check bowtie/bowtie2
    if len(path_to_aligner) != 0:
        path_to_aligner += "/"
    if bowtie2:
        aligner_version = get_executable_version(path_to_aligner+"bowtie2")
    else:
        aligner_version = get_executable_version(path_to_aligner+"bowtie")

    # check cutadapt
    if trim_reads:
        if len(path_to_cutadapt) != 0:
            path_to_cutadapt += "/"
        cutadapt_version = get_executable_version(path_to_cutadapt+"cutadapt")
        if parse_version(cutadapt_version) < parse_version("1.8"):
            print_error("cutadapt version %s found.\nmethylpy need at least cutadapt 1.9\nExit!\n"
                        %(cutadapt_version) )

    # check picard
    if remove_clonal:
        if len(path_to_picard) != 0:
            path_to_cutadapt += "/"
        # check java
        try:
            exec_name = "java"
            out = subprocess.check_output(shlex.split(exec_name+" -version"),
                                          stderr=subprocess.PIPE)
        except:
            print_error("java not found.\n"
                        +"methylpy need java to run picard to remove PCR duplicates\n"
                        +"Exit!\n")
        # check picard
        if not os.path.isfile(path_to_picard+"/picard.jar"):
            print_error("picard is not found at\n\""
                        +path_to_picard+"\"\n"
                        +"Exit!\n")

    return(True)

def convert_allc_to_bigwig(input_allc_file,
                           output_file,
                           reference_fasta,
                           mc_type="CGN",
                           bin_size = 100,
                           path_to_wigToBigWig="",
                           path_to_samtools="",
                           min_sites = 0,
                           min_cov = 0,
                           remove_chr_prefix=False
                           ):
    if not isinstance(mc_type, list):
        if isinstance(mc_type, str):
            mc_type = [mc_type]
        else:
            exit("mc_type must be a list of string(s)")

    if len(path_to_wigToBigWig):
        path_to_wigToBigWig += "/"

    mc_class = expand_nucleotide_code(mc_type)

    # prepare wig file
    cur_chrom = ""
    bin_end = bin_size
    bin_mc, bin_h, bin_site = 0, 0, 0
    g = open(output_file+".wig",'w')
    with open_allc_file(input_allc_file) as f:
        for line in f:
            fields = line.split("\t")
            if fields[3] not in mc_class:
                continue
            pos = int(fields[1])
            if cur_chrom != fields[0] or pos >= bin_end:
                if bin_h > 0 and bin_site >= min_sites and bin_h >= min_cov:
                    mc_level = str(float(bin_mc)/float(bin_h))
                    g.write("\t".join([cur_chrom,
                                       str(bin_end-bin_size),
                                       str(bin_end),
                                       mc_level])+"\n")
                # reset
                cur_chrom = fields[0]
                bin_end = int(pos/100+1) * bin_size
                bin_mc, bin_h, bin_site = 0, 0, 0

            # update
            bin_mc += int(fields[4])
            bin_h += int(fields[5])
            bin_site += 1

    if bin_h > 0 and bin_site >= min_sites and bin_h >= min_cov:
        mc_level = str(float(bin_mc)/float(bin_h))
        g.write("\t".join([cur_chrom,
                           str(bin_end-bin_size),
                           str(bin_end),
                           mc_level])+"\n")
    g.close()

    # chromosome size
    try:
        f = open(reference_fasta+".fai",'r')
    except:
        print("Reference fasta not indexed. Indexing.")
        try:
            subprocess.check_call(shlex.split(path_to_samtools
                                              +"samtools faidx "
                                              +reference_fasta))
            f = open(reference_fasta+".fai",'r')
        except:
            sys.exit("Reference fasta wasn't indexed, and couldn't be indexed. "
                     +"Please try indexing it manually and running methylpy again.")
    g = open(output_file+".chrom_size",'w')
    for line in f:
        fields = line.split("\t")
        if remove_chr_prefix and fields[0].startswith("chr"):
            fields[0] = fields[0][3:]
        g.write(fields[0]+"\t"+fields[1]+"\n")
    g.close()

    # generate bigwig file
    subprocess.check_call(shlex.split(path_to_wigToBigWig + "wigToBigWig -clip "
                                      +"%s.wig " %(output_file)
                                      +"%s.chrom_size " %(output_file)
                                      +output_file),stderr=subprocess.PIPE)
    subprocess.check_call(shlex.split("rm "+output_file+".wig "+output_file+".chrom_size"))


def filter_allc_file(allc_file,
                     output_file,
                     mc_type="CGN",
                     chroms = None,
                     compress_output=True,
                     min_cov=0,
                     buffer_line_number=100000):

    mc_class = expand_nucleotide_code(mc_type)
    # input & output
    f = open_allc_file(allc_file)
    if compress_output:
        output_fhandler = gzip.open(output_file,'wt')
    else:
        output_fhandler = open(output_file,'w')
    # begin
    line_counts = 0
    out = ""
    for line in f:
        fields = line.split("\t")
        if (chroms is None or fields[0] in chroms) \
           and fields[3] in mc_class \
           and int(fields[5]) >= min_cov:
            line_counts += 1
            out += line
        if line_counts > buffer_line_number:
            output_fhandler.write(out)
            line_counts = 0
            out = ""
    if line_counts > 0:
        output_fhandler.write(out)
        out = ""
    f.close()
    output_fhandler.close()

def merge_allc_files_minibatch(allc_files,
                               output_file,
                               mini_batch=20,
                               compress_output=True):
    
    # User input checks
    if not isinstance(allc_files, list):
        exit("allc_files must be a list of string(s)")

    try:
        merge_allc_files(allc_files=allc_files,
                         output_file=output_file,
                         compress_output=compress_output)
        index_allc_file(output_file)
        return 0
    except:
        print("Failed to merge all allc files once. Do minibatch merging")
        
    # init
    remaining_allc_files = list(allc_files[mini_batch:])
    output_tmp_file = output_file + ".tmp"
    merge_allc_files(allc_files=allc_files[:mini_batch],
                     output_file=output_file,
                     compress_output=compress_output)
    # batch merge
    while len(remaining_allc_files) > 0:
        processing_allc_files = [output_file]
        while len(remaining_allc_files) > 0 \
              and len(processing_allc_files) < mini_batch:
            processing_allc_files.append(remaining_allc_files.pop())
        merge_allc_files(allc_files=processing_allc_files,
                         output_file=output_tmp_file,
                         compress_output=compress_output)
        subprocess.check_call(["mv",output_tmp_file,output_file])
    # index output allc file
    index_allc_file(output_file)
    return 0

def merge_allc_files(allc_files,
                     output_file,
                     compress_output=True,
                     buffer_line_number=100000):
    #User input checks
    if not isinstance(allc_files, list):
        exit("allc_files must be a list of string(s)")

    # scan allc file to set up a table for fast look-up of lines belong
    # to different chromosomes
    fhandles = []
    chrom_pointer = []
    chroms = set([])
    for index,allc_file in enumerate(allc_files):
        fhandles.append(open_allc_file(allc_file))
        chrom_pointer.append(read_allc_index(allc_file))
        for chrom in chrom_pointer[index].keys():
            chroms.add(chrom)
    # output
    if compress_output:
        g = gzip.open(output_file,'wt')
    else:
        g = open(output_file,'w')
    # merge allc files
    chroms = list(map(str,chroms))
    for chrom in chroms:
        cur_pos = np.array([np.nan for index in range(len(allc_files))])
        cur_fields = []
        num_remaining_allc = 0
        # init
        for index,allc_file in enumerate(allc_files):
            fhandles[index].seek(chrom_pointer[index].get(chrom,0))
            line = fhandles[index].readline()
            fields = line.split("\t")
            cur_fields.append(None)
            if fields[0] == chrom:
                cur_pos[index] = int(fields[1])
                cur_fields[index] = fields
                num_remaining_allc += 1
        # merge
        out = ""
        line_counts = 0
        while num_remaining_allc > 0:
            mc, h = 0, 0
            c_info = None
            for index in np.where(cur_pos == np.nanmin(cur_pos))[0]:
                mc += int(cur_fields[index][4])
                h += int(cur_fields[index][5])
                if c_info is None:
                    c_info = "\t".join(cur_fields[index][:4])
                # update
                line = fhandles[index].readline()
                fields = line.split("\t")
                if fields[0] == chrom:
                    cur_pos[index] = int(fields[1])
                    cur_fields[index] = fields
                else:
                    cur_pos[index] = np.nan
                    num_remaining_allc -= 1
            # output
            out += c_info+"\t"+str(mc)+"\t"+str(h)+"\t1\n"
            line_counts += 1
            if line_counts > buffer_line_number:
                g.write(out)
                line_counts = 0
                out = ""
    if line_counts > 0:
        g.write(out)
        out = ""
    g.close()
    for index in range(len(allc_files)):
        fhandles[index].close()
    return 0

def get_index_file_name(allc_file):
    if allc_file[-4:] == ".tsv":
        index_file = allc_file[:-4]+".idx"
    elif allc_file[-7:] == ".tsv.gz":
        index_file = allc_file[:-7]+".idx"
    else:
        index_file = allc_file+".idx"
    return index_file

def index_allc_file_batch(allc_files,num_procs=1,no_reindex=False):
    if num_procs == 1:
        for allc_file in allc_files:
            index_allc_file(allc_file,no_reindex)
    else:
        pool = multiprocessing.Pool(min(num_procs,len(allc_files),100))
        for allc_file in allc_files:
            pool.apply_async(index_allc_file,(allc_file,no_reindex))
        pool.close()
        pool.join()
    return 0

def index_allc_file(allc_file,no_reindex=False):
    index_file = get_index_file_name(allc_file)
    # do not reindex if the index file is available
    if no_reindex and os.path.exists(index_file):
        return 0
    g = open(index_file,'w')
    f = open_allc_file(allc_file)
    cur_chrom = ""
    # check header
    line = f.readline()
    try:
        fields = line.split("\t")
        int(fields[1])
        int(fields[4])
        int(fields[5])
        # no header, continue to start from the beginning of allc file
        f.seek(0)
        cur_pointer = 0
    except:
        # find header, skip it
        cur_pointer = f.tell()
    # find chrom pointer
    while True:
        line = f.readline()
        if not line: break
        fields = line.split("\t")
        if fields[0] != cur_chrom:
            g.write(fields[0]+"\t"+str(cur_pointer)+"\n")
            cur_chrom = fields[0]
        cur_pointer = f.tell()
    f.close()
    g.close()
    return 0

def read_allc_index(allc_file):
    index_file = get_index_file_name(allc_file)
    try:
        f = open(index_file,'r')
    except:
        index_allc_file(allc_file)
        f = open(index_file,'r')
    chrom_pointer = {}
    for line in f:
        fields = line.rstrip().split("\t")
        chrom_pointer[fields[0]] = int(fields[1])
    f.close()
    return(chrom_pointer)

def remove_allc_index(allc_file):
    index_file = get_index_file_name(allc_file)
    subprocess.check_call(["rm",index_file])

def expand_nucleotide_code(mc_type):
    iub_dict = {"N":["A","C","G","T"],
                "H":["A","C","T"],
                "D":["A","G","T"],
                "B":["C","G","T"],
                "A":["A","C","G"],
                "R":["A","G"],
                "Y":["C","T"],
                "K":["G","T"],
                "M":["A","C"],
                "S":["G","C"],
                "W":["A","T"],
                "C":["C"],
                "G":["G"],
                "T":["T"],
                "A":["A"]}
    
    mc_class = list(mc_type) # copy
    if "C" in mc_type:
        mc_class.extend(["CGN", "CHG", "CHH","CNN"])
    elif "CG" in mc_type:
        mc_class.extend(["CGN"])

    for motif in mc_type:
        mc_class.extend(["".join(i) for i in
                         itertools.product(*[iub_dict[nuc] for nuc in motif])])
    return(set(mc_class))

def split_fastq_file(num_chunks, input_files, output_prefix):
    """
    This function mimics the unix split utility.
    """
    if not isinstance(input_files, list):
        if isinstance(input_files, str):
            input_files = [input_files]
        else:
            sys.exit("input_files must be a list of strings")
    file_handles = {}
    for index in range(0,num_chunks):
        file_handles[index]=open(output_prefix+str(index),'w')
    cycle = itertools.cycle(list(range(0,num_chunks)))
    total_reads=0
    for inputf in input_files:
        if inputf[-3:] == ".gz":
            f = gzip.open(inputf,'rt')
        elif inputf[-4:] == ".bz2":
            f = bz2.BZ2File(inputf,'r')
        else:
            f = open(inputf,'r')

        while True:
            current_file = next(cycle)
            # processing read id
            # remove any string after the first space character
            line = f.readline()
            if not line:
                break
            line = line.rstrip()
            file_handles[current_file].write(line.split(" ")[0]+"\n")
            total_reads += 1
            # seq
            line = f.readline()
            file_handles[current_file].write(line)
            # seq
            line = f.readline()
            file_handles[current_file].write("+\n")
            # qual
            line = f.readline()
            file_handles[current_file].write(line)
        f.close()

    for index in range(0,num_chunks):
        file_handles[index].close()

    return(total_reads)

def split_fastq_file_pbat(num_chunks, input_files, output_prefix):
    """
    This function mimics the unix split utility.
    """

    def reverse_complement(dna):
        complement = {"A":"T","C":"G","G":"C","T":"A","N":"N"}
        return("".join([complement[base] for base in reversed(dna)]))

    if not isinstance(input_files, list):
        if isinstance(input_files, str):
            input_files = [input_files]
        else:
            sys.exit("input_files must be a list of strings")
    file_handles = {}
    for index in range(0,num_chunks):
        file_handles[index]=open(output_prefix+str(index),'w')
    cycle = itertools.cycle(list(range(0,num_chunks)))
    total_reads=0
    for inputf in input_files:
        if inputf[-3:] == ".gz":
            f = gzip.open(inputf,'rt')
        elif inputf[-4:] == ".bz2":
            f = bz2.BZ2File(inputf,'r')
        else:
            f = open(inputf,'r')
                
        while True:
            current_file = next(cycle)
            # processing read id
            # remove any string after the first space character
            line = f.readline()
            if not line:
                break
            # read id
            line = line.rstrip()
            file_handles[current_file].write(line.split(" ")[0]+"\n")
            total_reads += 1
            # seq
            line = f.readline()
            line = line.rstrip()
            file_handles[current_file].write(reverse_complement(line)+"\n")
            # 
            line = f.readline()
            file_handles[current_file].write(line)
            # qual
            line = f.readline()
            line = line.rstrip()
            file_handles[current_file].write(line[::-1]+"\n")

        f.close()

    for index in range(0,num_chunks):
        file_handles[index].close()

    return(total_reads)

def split_mpileup_file(num_chunks,inputf,output_prefix):
    """
    This function mimics the unix split utility.
    preserve_order guarantees that the chunk files are split in the same order as the input file
    """
    total_lines = int(subprocess.check_output(["wc", "-l", inputf]).rstrip().split()[0])
    chunk_size = math.ceil(float(total_lines) / num_chunks)
    if chunk_size == 0:
        sys.exit("No lines in "+inputf)
    chunk_num = 0
    try:
        f = gzip.open(inputf,'rt')
        f.readline()
    except:
        try:
            f = bz2.BZ2File(inputf,'r')
            f.readline()
        except:
            f = open(inputf,'r')
    finally:
        f.seek(0)
    g = open(output_prefix+str(chunk_num),'w')
    count = 0
    prev_lines = collections.deque(maxlen=2)
    line = f.readline()
    while line:
        if count >= chunk_size:
            g.close()
            chunk_num+=1
            #This code looks a bit weird, but it's to handle a nasty edge case.
            #There's a problem if there is a C in one of the last two positions of a chunk
            #This information is needed at the end of the current file (to figure out the context of previous cytosines)
            #and the beginning of the next (to create new positions in the allc file). Consequently I write
            #the last two positions of a file again to the beginning of the next file.
            g = open(output_prefix+str(chunk_num),'w')
            g.write("".join(prev_lines))
            count = 0
        g.write(line)
        count+=1
        prev_lines.append(line)
        line = f.readline()
    g.close()
    f.close()
  
def parallel_count_lines(filename,
                         sample_chrom_pointer, chrom,
                         min_cov, mc_class):
    """
    Parallel helper to count lines in files. Used in split_files_by_position.
    """
    f = open_allc_file(filename)
    f.seek(sample_chrom_pointer[chrom])
    num_lines = 0

    for line in f:
        fields = line.split("\t")
        try:
            int(fields[1])
            int(fields[5])
        except:
            print("WARNING: One of the lines in "+filename+" is not formatted correctly:\n"+line+"Skipping it.")
            continue
        if fields[0] != chrom: break
        if fields[3] in mc_class and int(fields[5]) >= min_cov:
            num_lines += 1
    
    return((num_lines, filename))
    
def split_files_by_position(files,samples,
                            chunks,mc_class,
                            chrom_pointer,chrom,
                            num_procs=1, min_cov = 0,pool=False,max_dist=0,
                            weight_by_dist=False):
    """
    This function will split a group of files into chunks number of subfiles with the guarantee
    that two subfiles with the same suffix will contain the same range of positions. These files are
    assumed to be from the same chromosome. ASSUMES THAT ALLC FILES ARE SORTED.
    files - a list of file names
    chunks - the number of subfiles you wish to create
    mc_class - a list of 3 character strings indicating the nucleotide contexts you wish to consider.
    nrange - only consider positions between these coordinates
    max_dist - mc and h values from positions within max_dist of one another will be combined
        this helps leverage the correlation of methylation across a distance to make up for low
        coverage experiments.
    weight_by_dist - weight nearby position influence on mc & h values by their distance (otherwise sum them with equal weight)
    """
    #Note: only the first file will be truly evenly split, the other files will be split such that each
    #file contains the same range (but not necessarily the same number) of positions (e.g., from chr1:1000-2000)
    #num_lines = subprocess.check_output("wc -l "+files[0],shell=True).split(" ")[0]
    if not isinstance(files, list):
        if isinstance(files, str):
            files = [files]
        else:
            sys.exit("files must be a list")
    #use multiprocessing to count lines in each file and return (num_lines, filename) tuple
    if num_procs > 1:
        line_results = []
        #If a pool has not been passed to this function 
        if pool == False:
            pool_new = multiprocessing.Pool(num_procs)
            for input_file,sample in zip(files,samples):
                line_results.append(pool_new.apply_async(parallel_count_lines,
                                                     (input_file,chrom_pointer[sample],
                                                      chrom,min_cov,mc_class)))
            pool_new.close()
            pool_new.join()
        #If a pool has been passed to this function make sure it doesn't get closed by this funciton
        else:
            for input_file,sample in zip(files,samples):
                line_results.append(pool.apply_async(parallel_count_lines,
                                                     (input_file,chrom_pointer[sample],
                                                      chrom,min_cov,mc_class)))
        lines = [r.get() for r in line_results] #get() call blocks until all are finished so no need for wait()
    else:
        lines = []
        for input_file,sample in zip(files,samples):
            lines.append(parallel_count_lines(input_file,chrom_pointer[sample],
                                              chrom,min_cov,mc_class))

    #I have to guarantee that the smallest file is the one the chunk splitting is based on
    #if I don't do this, the smallest file might not have enough lines to be split into enough
    #chunks and then bad things happen.
    min_lines = 100000000000000
    min_file = None
    for index in range(len(lines)):
        if lines[index][0] > 0 and min_lines > lines[index][0]:
            min_lines,min_file = lines[index]
    chunk_size = math.ceil(float(min_lines) / chunks)
    if chunk_size == 0 or min_file is None:
        return 0
        #sys.exit("No lines match your range and/or mc_class criteria")
    min_sample = 0
    for input_file,sample in zip(files,samples):
        if input_file == min_file:
            min_sample = sample
            break
    chunk_num = 0
    count = 0
    #This list stores the position cutoffs for each chunk
    cutoffs = {}
    with open_allc_file(min_file) as f:
        f.seek(chrom_pointer[min_sample][chrom])
        for line in f:
            line = line.rstrip("\n")
            fields = line.split("\t")
            try:
                int(fields[1])
                int(fields[5])
            except:
                print("WARNING: One of the lines in "+min_file+" is not formatted correctly:\n"+line+"\nSkipping it.")
                continue
            if fields[0] != chrom:
                break
            elif int(fields[5]) >= min_cov and fields[3] in mc_class:
                count += 1
                if count > chunk_size:
                    cutoffs[chunk_num]=int(fields[1])
                    chunk_num += 1
                    count = 0
    #Want to make sure that for the last chunk all remaining lines get output
    #There's an edge case where the initial sample has positions which aren't as far
    #along the chromosome as other samples
    cutoffs[chunk_num] = 50000000000
    #enable multiprocessing
    if num_procs > 1:
        #If a pool has not been passed to this function 
        if pool == False:
            pool_new = multiprocessing.Pool(num_procs)
            for input_file,sample in zip(files,samples):
                pool_new.apply_async(parallel_split_files_by_position,(input_file,cutoffs,
                                                                   chrom_pointer[sample][chrom],
                                                                   chrom,
                                                                   mc_class),
                                 {"min_cov":min_cov, "max_dist":max_dist, "weight_by_dist":weight_by_dist})
            pool_new.close()
            pool_new.join()
        #If a pool has been passed to this function make sure it doesn't get closed by this funciton
        else:
            results = []
            for input_file,sample in zip(files,samples):
                results.append(pool.apply_async(parallel_split_files_by_position,
                                                (input_file,cutoffs,
                                                 chrom_pointer[sample][chrom],
                                                 chrom,
                                                 mc_class),
                                                {"min_cov":min_cov,
                                                 "max_dist":max_dist,
                                                 "weight_by_dist":weight_by_dist}))
            for result in results:
                result.wait()
    else:
        for input_file,sample in zip(files,samples):
            parallel_split_files_by_position(input_file,cutoffs,
                                             chrom_pointer[sample][chrom],
                                             chrom,
                                             mc_class,
                                             min_cov=min_cov,max_dist=max_dist,
                                             weight_by_dist=weight_by_dist)

def parallel_split_files_by_position(filen,cutoffs,
                                     sample_chrom_pointer,chrom,
                                     mc_class,
                                     min_cov=0,max_dist=0,
                                     weight_by_dist=False):
    chunk_num = 0
    f = open_allc_file(filen)
    f.seek(sample_chrom_pointer)
    #Just in case there's a header line in the file
    #I assume that if the position field can't be cast as an int
    #it must be a header
    fields_deque = collections.deque()
    added_values_deque = collections.deque()
    g = open(filen+"_"+chrom+"_"+str(chunk_num)+".tsv",'w')
    for line in f:
        line = line.rstrip("\n")
        fields = line.split("\t")
        try:
            int(fields[1])
            int(fields[5])
        except:
            print("WARNING: One of the lines in "+files[0]+" is not formatted correctly:\n"+line+"\nSkipping it.")
            continue
        if fields[0] != chrom: break
        if int(fields[1]) >= cutoffs[chunk_num]:
            g.close()
            chunk_num += 1
            g = open(filen+"_"+chrom+"_"+str(chunk_num)+".tsv",'w')
        if fields[3] in mc_class and int(fields[5]) >= min_cov:
            fields_deque.append(fields)
            mc = int(fields[4])
            h=int(fields[5])
            added_values_deque.append([mc,h])
            while len(fields_deque)> 0 and (int(fields_deque[-1][1]) - int(fields_deque[0][1])) >= max_dist:
                fields = fields_deque.popleft()
                added_values = added_values_deque.popleft()
                if added_values[0] != int(fields[4]):
                    g.write("\t".join(fields[0:4])+"\t"+str(added_values[0])+"\t"+str(added_values[1])+"\t1\n")
                else:
                    g.write("\t".join(fields[0:4])+"\t"+str(added_values[0])+"\t"+str(added_values[1])+"\t"+fields[6]+"\n")
            if len(fields_deque) > 1:
                for index in range(0,len(added_values_deque)-1):
                    if weight_by_dist:
                        weighting_factor = float(max_dist - abs(int(fields_deque[index][1])-
                                                                int(fields_deque[-1][1])))/max_dist
                    else:
                        weighting_factor = 1.0
                    added_values_deque[index][0]+= int(round(mc * weighting_factor))
                    added_values_deque[index][1]+=int(round(h * weighting_factor))
                    added_values_deque[-1][0]+=int(round(int(fields_deque[index][4])* weighting_factor))
                    added_values_deque[-1][1]+=int(round(int(fields_deque[index][5])* weighting_factor))

    if len(added_values_deque) > 0:
        fields = fields_deque.popleft()
        added_values = added_values_deque.popleft()
        if added_values[0] != int(fields[4]):
            g.write("\t".join(fields[0:4])+"\t"+str(added_values[0])+"\t"+str(added_values[1])+"\t1\n")
        else:
            g.write("\t".join(fields[0:4])+"\t"+str(added_values[0])+"\t"+str(added_values[1])+"\t"+fields[6]+"\n") 
    g.close()
    f.close()

def open_allc_file(allc_file):
    if allc_file[-3:] == ".gz":
        f = gzip.open(allc_file,'rt')
    elif allc_file[-4:] == ".bz2":
        f = bz2.BZ2File(allc_file,'rt')
    else:
        f = open(allc_file,'r')
    return(f)

def print_checkpoint(message):
    """
    Print message and current time
    """
    print(message)
    tabs = message.count("\t")
    print(("\t" * tabs) + time.asctime(time.localtime(time.time())) + "\n")
    sys.stdout.flush()


def print_error(error_message=""):
    sys.stderr.write("Error:\n" + error_message)
    sys.exit(1)

if __name__ == '__main__':
    pass
