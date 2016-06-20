"""
Make randomly generated allc files
"""

from sys import exit
from random import randint, uniform, shuffle
from numpy.random import binomial

def make_allc(output, context, chr, strand, dmr_size, between_size, dmr_diff, num_dmr, error, coverage, nonDMR_methyl, 
              num_between_pos, num_pos=2, num_files=2, balanced_mcs = True, replicates = False):
    """
    Write 2 allc files with num_dmr DMRs that have a methylation difference of dmr_diff and DMR sizes of dmr_size. 
    The h value is determined by coverage. All C's have the same chr, context, and strand.
    
    output - filename used for generated allc file
    context - methylation context for all C's
    chr - string of chr for all blocks
    strand - strand used for all blocks
    dmr_size - position span of DMR blocks
    between_size - distance between DMR blocks
    dmr_diff - percentage difference between two sample DMRs
    num_dmr - number of DMRs in allc
    error - chance of introducing errors in methylation
    coverage - number of sites that each position is covered by
    num_pos - number of positions within a block, minimum is 2
    num_files - number of allc files to output, default 2
    balanced_mcs - indicates whether you want an equal number of max & min blocks among the files or random number
    replicates - indicates whether you want replicate files made
    nonDMR_methyl - true value used for all nonDMR methylation
    num_between_pos - number of positions between each DMR block
    """
    if dmr_diff >= 1 or dmr_diff < 0:
        exit("dmr_diff not within bounds")
    if num_pos < 2 or num_pos > dmr_size:
        exit("Number of positions must be between 2 and dmr_size")
    if num_files < 2:
        exit("num_files must be greater than 2")
    if num_between_pos > between_size:
        exit("The number of positions between DMRs must be less than or equal to the between_size")
    if nonDMR_methyl < 0 or nonDMR_methyl > 1:
        exit("nonDMR_methyl must be between 0 and 1")
    
    if dmr_diff + 2*error > 1:
        raise Exception("(dmr_diff + 2*error) must be <= 1.")
    
    if dmr_diff + error >= 1:
        max_vals = (dmr_diff, 1 - error)
    else:
        max_vals = (dmr_diff + error, 1 - error)
    
    pos = 1
    chr = str(chr)
    extra_pos = num_pos - 2 #subtract 2 for start and end pos, these are anchored at edges of block
      
    output_files = [] #list of files
    for x in xrange(num_files):
        output_files.append(open(output+str(x), 'w'))
    
    if replicates:
        output_replicate = [] #list of files
        for x in xrange(num_files):
            output_replicate.append(open(output+str(x)+"_replicate", 'w'))
    
    file_category = [] #list of categories for files that indicates whether it gets maximum methylation or minimum
    if balanced_mcs:
        for x in xrange(num_files):
            file_category.append(x%2)
        shuffle(file_category)
    else:
        for x in xrange(num_files):
            file_category.append(randint(0,1))
        if 0 not in file_category or 1 not in file_category: #minimum 1 of each type
            file_category[0] = abs(file_category[0] - 1)
            
    g = open(output + "_allc_sim_summary", 'w')
    header = "chr\tstart\tstop\tstrand\tcontext\tmax_methyl_value\tmin_methyl_value"
    for i in xrange(num_files):
        header += "\t" + "file" + str(i)
    g.write(header+"\n")
    
    in_block = True
    max_methyl = uniform(max_vals[0], max_vals[1]) #get true value of max and min for this block
    min_methyl = uniform(max_methyl-dmr_diff, max_methyl-dmr_diff)
    if min_methyl < 0: min_methyl = 0
    
    # Loop goes through twice for each DMR, once for start pos, once for end pos
    for i in xrange(num_dmr*2):
          
        for index in xrange(len(output_files)):
            f = output_files[index]
            if file_category[index]:
                max_err = uniform(max_methyl-error, max_methyl+error) #introduce error for this pos
                mc, h = (binomial(coverage, max_err), coverage)
                output = [chr, str(pos), strand, context, str(mc), str(h), "1"]
                f.write("\t".join(output) + "\n")
                if replicates:
                    max_err = uniform(max_methyl-error, max_methyl+error) #introduce error for this pos
                    mc, h = (binomial(coverage, max_err), coverage)
                    output = [chr, str(pos), strand, context, str(mc), str(h), "1"]
                    output_replicate[index].write("\t".join(output) + "\n")                
            else:
                min_err = uniform(min_methyl-error, min_methyl+error) #introduce error for this pos
                mc, h = (binomial(coverage, min_err), coverage)
                output = [chr, str(pos), strand, context, str(mc), str(h), "1"]
                f.write("\t".join(output) + "\n")
                if replicates:
                    min_err = uniform(min_methyl-error, min_methyl+error) #introduce error for this pos
                    mc, h = (binomial(coverage, min_err), coverage)
                    output = [chr, str(pos), strand, context, str(mc), str(h), "1"]
                    output_replicate[index].write("\t".join(output) + "\n")                    
        
        if in_block:
            #write extra positions inside this block with random # pos between them
            save_pos = pos
            for j in xrange(extra_pos):
                pos += (dmr_size - 1) / (extra_pos+1)
                for index in xrange(len(output_files)):
                    f = output_files[index]
                    if file_category[index]:
                        max_err = uniform(max_methyl-error, max_methyl+error) #introduce error for this pos
                        mc, h = (binomial(coverage, max_err), coverage)
                        output = [chr, str(pos), strand, context, str(mc), str(h), "1"]
                        f.write("\t".join(output) + "\n")
                        if replicates:
                            max_err = uniform(max_methyl-error, max_methyl+error) #introduce error for this pos
                            mc, h = (binomial(coverage, max_err), coverage)
                            output = [chr, str(pos), strand, context, str(mc), str(h), "1"]
                            output_replicate[index].write("\t".join(output) + "\n")                            
                    else:
                        min_err = uniform(min_methyl-error, min_methyl+error) #introduce error for this pos
                        mc, h = (binomial(coverage, min_err), coverage)
                        output = [chr, str(pos), strand, context, str(mc), str(h), "1"]
                        f.write("\t".join(output) + "\n")
                        if replicates:
                            min_err = uniform(min_methyl-error, min_methyl+error) #introduce error for this pos
                            mc, h = (binomial(coverage, min_err), coverage)
                            output = [chr, str(pos), strand, context, str(mc), str(h), "1"]
                            output_replicate[index].write("\t".join(output) + "\n")                            
            pos = save_pos + dmr_size - 1
        else:
            #write summary of previous block
            summary_out = [chr, str(pos-dmr_size+1), str(pos), strand, context, str(max_methyl), str(min_methyl)]
            for cat in file_category:
                summary_out.append(str(cat)) 
            g.write("\t".join(summary_out) + "\n")
            
            #write between block positions
            save_pos = pos
            for k in xrange(num_between_pos):
                pos += between_size / num_between_pos
                for index in xrange(len(output_files)):
                    f = output_files[index]
                    err = uniform(nonDMR_methyl-error, nonDMR_methyl+error) #introduce error for this pos
                    mc, h = (binomial(coverage, err), coverage)
                    output = [chr, str(pos), strand, context, str(mc), str(h), "1"]
                    f.write("\t".join(output) + "\n")
                    if replicates:
                        err = uniform(nonDMR_methyl-error, nonDMR_methyl+error) #introduce error for this pos
                        mc, h = (binomial(coverage, err), coverage)
                        output = [chr, str(pos), strand, context, str(mc), str(h), "1"]
                        output_replicate[index].write("\t".join(output) + "\n")  
            pos = save_pos + between_size + 1
            
            max_methyl = uniform(max_vals[0], max_vals[1]) #true value of max and min for next block
            min_methyl = uniform(max_methyl-dmr_diff, max_methyl-dmr_diff)
            if min_methyl < 0: min_methyl = 0            
            #shuffle categories for next block
            if balanced_mcs:
                shuffle(file_category)
            else:
                file_category[:] = []
                for x in xrange(num_files):
                    file_category.append(randint(0,1))
                if 0 not in file_category or 1 not in file_category: #minimum 1 of each type
                    file_category[0] = abs(file_category[0] - 1)
        in_block = not in_block
    
    for f in output_files:
        f.close()
    if replicates:
        for f2 in output_replicate:
            f2.close()
    g.close()
