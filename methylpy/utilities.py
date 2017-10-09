'''
Created on Dec 16, 2011

@author: Matt
'''

import sys
import pdb
import math
import shlex

try:
    import numpy as np
except:
    sys.exit("methylpy.utilities requires the numpy module")

try:
    import subprocess
except:
    sys.exit("methylpy.utilities requires the subprocess module")    

try:
    import time
except:
    sys.exit("methylpy.utilities requires the time module")
try:
    import operator
except:
    sys.exit("methylpy.utilities requires the operator module")
try:
    import itertools
except:
    sys.exit("methylpy.utilities requires the itertools module")
try:
    import multiprocessing
except:
    sys.exit("methylpy.utilities requires the multiprocessing module")
try:
    import bz2
except Exception,e:
    exc_type, exc_obj, exc_tb = sys.exc_info()
    print(exc_type, exc_tb.tb_lineno)
    print e
    sys.exit("methylpy.DMRfind requires the bz2 module")

try:
    import gzip
except Exception,e:
    exc_type, exc_obj, exc_tb = sys.exc_info()
    print(exc_type, exc_tb.tb_lineno)
    print e
    sys.exit("methylpy.utilities requires the gzip module")
try:
    import collections
except Exception,e:
    exc_type, exc_obj, exc_tb = sys.exc_info()
    print(exc_type, exc_tb.tb_lineno)
    print e
    sys.exit("methylpy.utilities requires the collections module") 
def browser_file(collapsed, uncollapsed):
    """
    NOT FINISHED
    Create genome browser files for each sample from collapsed and uncollapsed files.
    format: chr, start, stop, methylation level
    """
    #Make dictionary of all samples and where they are methylated
    methyl_coords = {}
    f = open(collapsed, 'r')
    f.readline()    
    for line in f:
        line = line.strip()
        line = line.split('\t')
        if len(line) < 5: #if there's no samples, continue
            continue
        methyl_coords[(line[0], line[1], line[2])] = {} #[chr, start, stop]
        methylated_samples = line[4].split(",") #list of methylated samples
        if methylated_samples!=['']:
            for sample in methylated_samples:
                methyl_coords[(line[0], line[1], line[2])][sample] = None                       
    f.close()
    #Find methylation level for each location in a sample's list
    f2 = open(uncollapsed, 'r')
    
    fields = {} #mapping of header values to field number
    header = f2.readline()
    header = header.strip().split('\t')
    for value in header:
        fields[value] = header.index(value)
    
    for line in f2:
        line = line.strip()
        line = line.split("\t")
        #create list of blocks that contained this position
        coord_keys = [x for x in methyl_coords.keys() if x[0]==line[0] and line[1]>=x[1] and line[1]<=x[2]]
        for key in coord_keys:
            samples = methyl_coords[key].keys() #all the methylated samples that were in this region
            for sample in samples:
                if not methyl_coords[key][sample]:
                    methyl_coords[key][sample] = [0, 0]
                pdb.set_trace()
                methyl_coords[key][sample][0] += int(line[fields["mc_" + sample]]) 
                methyl_coords[key][sample][1] += int(line[fields["h_" + sample]])
    f2.close()
    
    #Write output to files
    file_handles = {}
    for coord in methyl_coords:
        for sample in methyl_coords[coord]:
            if sample not in file_handles:
                file_handles[sample] = open("_".join([sample, "browserfile.bed"]), 'w')
                file_handles[sample].write("chr\tstart\tstop\tmethyl_level\n")
            try:
                methyl_level = float(methyl_coords[coord][sample][0]) / methyl_coords[coord][sample][1]
            except:
                methyl_level = "NA"
            line = "\t".join([str(coord[0]), str(coord[1]), str(coord[2]), str(methyl_level), "\n"])
            file_handles[sample].write(line)
    
    for filen in file_handles:
        file_handles[filen].close()

def expand_nucleotide_code(mc_type):
    iub_dict = {"N":["A","C","G","T"],
                "H":["A","C","T"],
                "C":["C"],
                "G":["G"],
                "T":["T"],
                "A":["A"]}
    
    mc_class = list(mc_type) # copy
    #
    if "C" in mc_type:
        mc_class.extend(["CG", "CHG", "CHH","CNN"])
	if "CG" in mc_type:
	    mc_class.extend(["CGN"])
    #       
    for motif in mc_type:
	mc_class.extend(["".join(i) for i in
			itertools.product(*[iub_dict[nuc] for nuc in motif])])
    return(set(mc_class))

def split_allc_window(num_chunks, inputf, output_prefix, max_dist, jump_size, window_size):
    """
    This function mimics the unix split utility. This splits the allc file into num_chunks files.
        max_dist = block must meet distance from either side of block from center
        jump_size = amount to jump centers when doing new overlaps
        window_size = the amount of positions to consider in each window to calc mc/h value
        
    Assumption that the input file only has one chromosome information inside
    """
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
        f.seek(0)
    
    total_lines = int(subprocess.check_output(["wc", "-l", inputf]).rstrip().split()[0])
    chunk_size = math.ceil(float(total_lines) / num_chunks)
    if chunk_size == 0:
        sys.exit("No lines in allc files")
    line_queue = [] #maintains list of previous lines to draw from for windows
    chunk_num = 0
    g = open(output_prefix+str(chunk_num),'w')
    
    #start with whole window loaded into queue
    for i in xrange(0, window_size):
        line = f.readline()
        line = line.strip().split('\t')
        line_queue.append(line)
    
    print_count = 0
    while True:   
        center = line_queue[window_size/2] #add all mc and h if it meets distance requirements
        if center == "NA": #end loop if center no longer exists
            break
        mc = 0
        h = 0
        for line in filter(lambda x: x!="NA", line_queue): #filter out NA's to calc mc & h
            if math.fabs(int(line[1]) - int(center[1])) <= max_dist: #if it meets distance req
                mc += int(line[4])
                h += int(line[5])
        
        if print_count >= chunk_size: #print out center information for this window
            g.close()
            chunk_num+=1
            g = open(output_prefix+str(chunk_num),'w')
            print_count = 0
        g.write("\t".join([center[0], center[1], center[2], center[3], str(mc), str(h), center[6]]) + "\n")
        print_count+=1
        
        #shift window by amount specified by jump_size
        for i in xrange(0, jump_size):
            line = f.readline()
            if not line: #very last window will have missing lines, fill in with "NA"
                line_queue.pop(0)
                line_queue.append("NA")
            else:
                line = line.strip().split('\t')
                line_queue.pop(0) #remove oldest line from queue
                line_queue.append(line) #add line to queue
    
    g.close()
    f.close()

def split_allc_file(num_chunks, inputf, output_prefix, window=False):
    """
    This function mimics the unix split utility.
    """
    total_lines = int(subprocess.check_output(["wc", "-l", inputf]).rstrip().split()[0])
    chunk_size = math.ceil(float(total_lines) / num_chunks)
    if chunk_size == 0:
        sys.exit("No lines in allc_file "+inputf)
    chunk_num = 0
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
        f.seek(0)
    g = open(output_prefix+str(chunk_num),'w')
    count = 0
    line = f.readline()
    while line:
        if count >= chunk_size:
            g.close()
            chunk_num+=1
            g = open(output_prefix+str(chunk_num),'w')
            count = 0
        g.write(line)
        count+=1
        line = f.readline()
    g.close()
    f.close()

def split_fastq_file(num_chunks, input_files, output_prefix):
    """
    This function mimics the unix split utility.
    """
    if not isinstance(input_files, list):
        if isinstance(input_files, basestring):
            input_files = [input_files]
        else:
            sys.exit("input_files must be a list of strings")
    file_handles = {}
    for index in xrange(0,num_chunks):
        file_handles[index]=open(output_prefix+str(index),'w')
    cycle = itertools.cycle(range(0,num_chunks))
    total_reads=0
    for inputf in input_files:
        if inputf[-3:] == ".gz":
            f = gzip.open(inputf,'r')
        elif inputf[-4:] == ".bz2":
            f = bz2.BZ2File(inputf,'r')
        else:
            f = open(inputf,'r')
                
        while True:
            current_file = cycle.next()
            # processing read id
            # remove any string after the first space character
            line = f.readline()
            if not line:
                break
            line = line.rstrip()
            file_handles[current_file].write(line.split(" ")[0]+"\n")
            total_reads += 1
            for index in xrange(0,3):
                line = f.readline()
                file_handles[current_file].write(line)
        f.close()

    for index in xrange(0,num_chunks):
        file_handles[index].close()

    return total_reads

def split_fastq_file_pbat(num_chunks, input_files, output_prefix):
    """
    This function mimics the unix split utility.
    """

    def reverse_complement(dna):
        complement = {"A":"T","C":"G","G":"C","T":"A","N":"N"}
        return "".join([complement[base] for base in dna[::-1]])

    if not isinstance(input_files, list):
        if isinstance(input_files, basestring):
            input_files = [input_files]
        else:
            sys.exit("input_files must be a list of strings")
    file_handles = {}
    for index in xrange(0,num_chunks):
        file_handles[index]=open(output_prefix+str(index),'w')
    cycle = itertools.cycle(range(0,num_chunks))
    total_reads=0
    for inputf in input_files:
        if inputf[-3:] == ".gz":
            f = gzip.open(inputf,'r')
        elif inputf[-4:] == ".bz2":
            f = bz2.BZ2File(inputf,'r')
        else:
            f = open(inputf,'r')
                
        while True:
            current_file = cycle.next()
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

    for index in xrange(0,num_chunks):
        file_handles[index].close()

    return total_reads

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
        f = gzip.open(inputf,'r')
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
    try:
        f = gzip.open(filename,'r')
        f.readline()
    except:
        try:
            f = bz2.BZ2File(filename,'r')
            f.readline()
        except:
            f = open(filename,'r')
    finally:
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
    
    return (num_lines, filename)
    
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
        if isinstance(files, basestring):
            files = [files]
        else:
            sys.exit("files must be a list")
    #I have to guarantee that the smallest file is the one the chunk splitting is based on
    #if I don't do this, the smallest file might not have enough lines to be split into enough
    #chunks and then bad things happen.
    min_lines = 100000000000000

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
    min_entry = min(lines, key=lambda x: x[0])
    min_lines = min(min_entry[0], min_lines)
    min_file = min_entry[1]
    min_sample = 0
    for input_file,sample in zip(files,samples):
        if input_file == min_file:
            min_sample = sample
            break
    chunk_size = math.ceil(float(min_lines) / chunks)
    if chunk_size == 0:
        sys.exit("No lines match your range and/or mc_class criteria")    
    chunk_num = 0
    count = 0
    #This list stores the position cutoffs for each chunk
    cutoffs = {}
    with open(min_file,'r') as f:
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
    try:
        f = gzip.open(filen,'r')
        f.readline()
    except:
        try:
            f = bz2.BZ2File(filen,'r')
            f.readline()
        except:
            f = open(filen,'r')
    finally:
        f.seek(sample_chrom_pointer)
    #Just in case there's a header line in the file
    #I assume that if the position field can't be cast as an int
    #it must be a header
    fields_deque = collections.deque()
    added_values_deque = collections.deque()
    g = open(filen+"_"+chrom+"_"+str(chunk_num),'w')
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
            g = open(filen+"_"+chrom+"_"+str(chunk_num),'w')
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
                for index in xrange(0,len(added_values_deque)-1):
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

def print_checkpoint(message):
    """
    Print message and current time
    """
    print message
    tabs = message.count("\t")
    print ("\t" * tabs) + time.asctime(time.localtime(time.time())) + "\n"
    sys.stdout.flush()


def print_error(error_message=""):
    sys.stderr.write("Error:\n" + error_message)
    sys.exit(1)

if __name__ == '__main__':
    pass
