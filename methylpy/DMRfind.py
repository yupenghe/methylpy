'''
Created on Dec 16, 2011

@author: Matt
'''
try:
    from sys import exc_info, exit
except Exception,e:
    exc_type, exc_obj, exc_tb = exc_info()
    print(exc_type, exc_tb.tb_lineno)
    print e
    exit("methylpy.DMRfind requires exec_info and exit from the sys module")
try:
    from pdb import set_trace
except Exception,e:
    exc_type, exc_obj, exc_tb = exc_info()
    print(exc_type, exc_tb.tb_lineno)
    print e
    exit("methylpy.DMRfind requires set_trace from the pdb module")
try:
    from subprocess import check_call,check_output
except Exception,e:
    exc_type, exc_obj, exc_tb = exc_info()
    print(exc_type, exc_tb.tb_lineno)
    print e
    exit("methylpy.DMRfind requires check_call and check_output from the subprocess module")

try:
    from shlex import split
except Exception,e:
    exc_type, exc_obj, exc_tb = exc_info()
    print(exc_type, exc_tb.tb_lineno)
    print e
    exit("methylpy.DMRfind requires split from the shlex module")
try:
    from glob import glob
except Exception,e:
    exc_type, exc_obj, exc_tb = exc_info()
    print(exc_type, exc_tb.tb_lineno)
    print e
    exit("methylpy.DMRfind requires the glob from the glob module")
try:
    from math import ceil, log10
except Exception,e:
    exc_type, exc_obj, exc_tb = exc_info()
    print(exc_type, exc_tb.tb_lineno)
    print e
    exit("methylpy.DMRfind requires ceil and log10 from the math module")

try:
    from multiprocessing import Pool
except Exception,e:
    exc_type, exc_obj, exc_tb = exc_info()
    print(exc_type, exc_tb.tb_lineno)
    print e
    exit("methylpy.DMRfind requires pool from the multiprocessing module")

try:
    from methylpy.utilities import print_checkpoint,expand_nucleotide_code,split_files_by_position
except Exception,e:
    exc_type, exc_obj, exc_tb = exc_info()
    print(exc_type, exc_tb.tb_lineno)
    print e
    exit("methylpy.DMRfind requires print_checkpoint, expand_nucleotide_code, and split_files_by_position from the methylpy.utilities module")

try:
    from scipy.stats import scoreatpercentile
except Exception,e:
    exc_type, exc_obj, exc_tb = exc_info()
    print(exc_type, exc_tb.tb_lineno)
    print e
    exit("methylpy.DMRfind requires scipy.stats module")

import subprocess
import shlex
import gzip

def DMRfind(allc_files, samples,
            mc_type, chroms,
            output_prefix,
            num_procs=1, 
            min_cov=0,keep_temp_files=False,mc_max_dist=0,
            dmr_max_dist=250,resid_cutoff=.01,sig_cutoff=.01,
            num_sims=3000,num_sig_tests=100,seed=-1,
            min_num_dms=0,collapse_samples=False,
            sample_category=False, min_cluster=0,
            max_iterations=1000,convergence_diff=1,
            buffer_line_number=100000):
    """
    This function will take a set of allc files, and look for differentially methylated regions. Note that in the output file, 
    -1 is used to represent missing data.
    
    mc_type is a list of the mc nucleotide contexts for which you want to look for DMRs. These classifications 
        may use the wildcards "H" (indicating anything but a G) and "N" (indicating any nucleotide). For example,
        CAH would look for DMRs in methylated cytosines followed by AT, AC, or AA.
    region_dict is a dictionary that has keys of the chromosomes that you want to search that points to a list 
        of [start,end] for that chromosome. The list elements should be positive integers.
    samples is a list of sample names to help find the allc files. See path_to_allc for more details.
    path_to_allc is a string indicating the path to the tab separated files containing methylation information
        for all C nucleotides in the genome. These files will be of the format: chromosome, position, strand, nucleotide context,
        reads with a methylated C at this position, total reads covering this position, and a binary value 0/1 indicating whether or not this
        site was statistically significantly methylated (e.g., as determined by a binomial test based on factors like sequencing
        error and non-conversion rate). A file for each chromosome is expected and their names will be in this format:
        path_to_allc/allc_<sample name>_<chromosome number>.tsv. Note that chromosome number should be an integer or letter (i.e.,
        exclude the string "chr")
    num_procs is the number of processors you wish to use to parallelize this function
    save_result is a string indicating the prefix for various result files produced by this function
    min_cov is the minimum number of reads that must cover a site for it to be considered in DMR finding
    keep_temp_files indicates that you'd like to keep all the intermediate files this function generates along
        the way. This can be useful for debugging.
    num_sims indicates the number of permutation tests you'd like to run to estimate the p-values of the differential methylation tests
    num_sig_tests is an integer indicating how many permutations can return a result more significant than the original statistic
        before permutation testing is abandoned for that particular site.
    mc_max_dist is an integer indicating the maximum distance two sites can be from one another for their methylation counts to be combined
        this option helps with low coverage experiments where you may want to leverage the correlation of methylation between sites to get 
        more statistical power.
    *Options that are passed through to collapse_dmr_windows*
    dmr_max_dist is the maximum distance two significant sites can be to be included in the same block
    resid_cutoff - If this option is specified, not only will a result have to be significant to be included in a window, but it will also,
        have to show deviations in the contingency table in the same direction as the rest of the window. For example, if sample A is generally
        showing a degree of methylation higher than expected and sample B lower than expected, any new site will have to have these same properties.
        Furthermore, these deviations have to be at least as extreme as resid_cutoff (in the positive or negative direction). This value is determined
        by looking at the distribution of residuals at sites with non-significant p-values. The value specified here indicates the percentile (0.00-1.00) 
        in the distribution of non-significant p-values to look for to determine a residual cutoff.
    sig_cutoff - the cutoff for determining whether a row is significant or not
    seed - a seed to provide to the random number generator for permutation testing. Only change this if you are debugging and want to make sure
        the permutation output is consistent
    min_num_dms - the minimum number of differentially methylated sites that a differentially methylated region needs to contain to be reported
    collapse_samples - a list of samples for collapsing blocks. If the methylation status (hypermethylated/hypomethylated) of one of the samples switch
        for a DMS, a new block will be created. NOTE: it is best to run DMRfind with all the samples you will potentially want to consider, and then
        to run collapse_dmr_windows for the individual cases you are interested in. In other words, run DMRfind with the sample argument set to:
        ["A","B","C"] and then run collapse_dmr_windows 3 times with collapse_samples set to ["A","B"],["A","C"], and finally ["B","C"]
    sample_category - A list of categories that each respective sample belongs to; the categories must begin at 0 and increase by
        1 for each category added. ex: samples [A,B,C] categories [0,1,2] or categories [0, 1, 0] 
    min_cluster - The minimum number of each sample category that must be present in every block that is output.
    max_iterations is the maximum number of iterations performed by the algorithm described in the paper cited above
    convergence_diff determines when the algorithm will terminate. When the current m0 estimate and the last m0 estimate differ by no
        more than convergence_diff.
    """
    
    #User input checks
    if not isinstance(allc_files, list):
        exit("mc_type must be a list of string(s)")
    if not isinstance(mc_type, list):
        if isinstance(mc_type, basestring):
            mc_type = [mc_type]
        else:
            exit("mc_type must be a list of string(s)")
    if not isinstance(samples, list):
        exit("samples must be a list of string(s)")
    try:
        num_procs = int(num_procs)
    except:
        exit("In DMRfind, num_procs must be an integer")

    try:
        mc_max_dist = int(mc_max_dist)
    except:
        exit("In DMRfind, mc_max_dist must be an integer")
    try:
        dmr_max_dist = int(dmr_max_dist)
    except:
        exit("In DMRfind, dmr_max_dist must be an integer")

    try:
        min_cov = int(min_cov)
    except:
        exit("In DMRfind, min_cov must be an integer")
    if isinstance(chroms,list) == False:
        exit("chroms must be a list of string(s)")
    chroms = map(str,chroms)
        
    if collapse_samples != False:
        if not isinstance(collapse_samples, list):
            exit("collapse_samples must be a list of strings")
        if sample_category and not isinstance(sample_category, list):
            exit("sample_category must be a list of strings")
        for sample in collapse_samples:
            if sample not in samples:
                exit("There is a sample in collapse_samples that is not in samples."+
                     "collapse_samples MUST be a subset of samples.")   
    elif sample_category != False:
        exit("In order to use sample_category, you must specify a corresponding"+
             " list of samples in collapse_samples!")

    if mc_max_dist <0:
        exit("In DMRfind, mc_max_dist must be greater than or equal to 0")
    if dmr_max_dist <0:
        exit("In DMRfind, dmr_max_dist must be greater than 0")
    
    #This code creates all variations of the shorthand C contexts (e.g., CHG->CHG,CAG,CCG,CTG)
    mc_class = expand_nucleotide_code(mc_type)
    # scan allc file to set up a table for fast look-up of lines belong
    # to different chromosomes
    chrom_pointer = {}
    for allc_file,sample in zip(allc_files,samples):
        cp_dict = {}
        with open(allc_file,'r') as f:
            cur_chrom = ""
            cur_pointer = 0
            while True:
                line = f.readline()
                if not line: break
                fields = line.split("\t")
                if fields[0] != cur_chrom:
                    cp_dict[fields[0]] = cur_pointer
                    cur_chrom = fields[0]
                cur_pointer = f.tell()
        chrom_pointer[sample] = cp_dict

    if num_procs > 1:
        pool = Pool(num_procs)
    else:
        pool = False

    try:
        for chr_key in chroms:
            chrom = str(chr_key).replace("chr","")
            results = []
            print_checkpoint("Splitting allc files for chromosome "+str(chrom))
            split_files_by_position(allc_files,samples,
                                    num_procs,mc_class,
                                    chrom_pointer = chrom_pointer,
                                    chrom=chrom,
                                    num_procs=num_procs,min_cov=min_cov,
                                    pool=pool,max_dist=mc_max_dist)
            print_checkpoint("Running rms tests for chromosome "+str(chrom))
            if num_procs > 1:
                for chunk in xrange(0,num_procs):
                    filenames = []
                    for allc_file in allc_files:
                        filenames.extend(glob(allc_file+"_"+chrom+"_"+str(chunk)))
                    if len(filenames) == 0:
                        print "Nothing to run for chunk "+str(chunk)
                        continue
                    pool.apply_async(run_rms_tests,
                                     (filenames,
                                      output_prefix+"_rms_results_for_"+str(chrom)+
                                      "_chunk_"+str(chunk)+".tsv",samples),
                                     {"min_cov":min_cov,
                                      "num_sims":num_sims,
                                      "num_sig_tests":num_sig_tests,
                                      "seed":seed,"keep_temp_files":keep_temp_files})
            else:
                filenames = []
                for allc_file in allc_files:
                    filenames.extend(glob(allc_file+"_"+chrom+"_0"))
                if len(filenames) == 0:
                    print "Nothing to run for chunk "+str(chunk)
                    continue
                run_rms_tests(filenames,output_prefix+"_rms_results_for_"+str(chrom)+"_chunk_0.tsv",
                              samples,min_cov=min_cov,num_sims=num_sims,num_sig_tests=num_sig_tests,
                              seed=seed,keep_temp_files=keep_temp_files)
        if pool != False:
            pool.close()
            pool.join()

    except Exception,e:
        exc_type, exc_obj, exc_tb = exc_info()
        print(exc_type, exc_tb.tb_lineno)
        print e
        try:
            pool.terminate()
            pool.join()
        except:
            pass
        exit("Running RMS tests failed.")
        
    print_checkpoint("Merging sorted "+output_prefix+"_rms_results.tsv files.")
    header = "\t".join(["\t".join(["chr","pos","strand","mc_class","pvalue"]),
                        "\t".join(["mc_"+sample for sample in samples]),
                        "\t".join(["h_"+sample for sample in samples]),
                        "\t".join(["frac_"+sample for sample in samples]),
                        "\t".join(["mc_residual_"+sample for sample in samples]),
                        "\t".join(["uc_residual_"+sample for sample in samples]),
                        "num_simulations_sig\tnum_simulations_run"+"\n"])
    #I put this up here because it's actually a lot harder to prepend a header file than you might think
    g = open(output_prefix+"_rms_results.tsv",'w')
    g.write(header)
    for chr_key in sorted(chroms):
        chrom = str(chr_key).replace("chr","")
        for chunk in xrange(0,num_procs):
            try:
                with open(output_prefix+"_rms_results_for_"+chrom+"_chunk_"+str(chunk)+".tsv",'r') as f:
                    for line in f:
                        g.write(line)
            except:
                pass
    g.close()
    if keep_temp_files == False:
        basecmd = ['rm'] 
        file_paths = glob(
            output_prefix+"_rms_results_for_*_chunk_[0-9].tsv") + glob(
                output_prefix+"_rms_results_for_*_chunk_[0-9][0-9].tsv") + glob(
                    output_prefix+"_rms_results_for_*_chunk_[0-9][0-9][0-9].tsv")
        if file_paths:
            try:
                check_call(basecmd + file_paths)
            except:
                pass
    print_checkpoint("Begin FDR Correction")
    pvalue_cutoff=histogram_correction_DMRfind(output_prefix+"_rms_results.tsv",
                                               num_sims,num_sig_tests,
                                               target_fdr =sig_cutoff,
                                               max_iterations=max_iterations,
                                               convergence_diff=convergence_diff)

    print_checkpoint("Calculating Residual Cutoff")
    resid_cutoff = get_resid_cutoff(resid_cutoff, pvalue_cutoff,
                                    len(samples), output_prefix+"_rms_results.tsv")

    print_checkpoint("Begin Defining Windows")
    collapse_dmr_windows(output_prefix+"_rms_results.tsv",
                         output_prefix+"_rms_results_collapsed.tsv",
                         column=4,sig_cutoff=pvalue_cutoff,
                         max_dist=dmr_max_dist,
                         resid_cutoff=resid_cutoff,
                         min_num_dms=min_num_dms,
                         collapse_samples=collapse_samples,
                         sample_category=sample_category,
                         min_cluster=min_cluster)
    get_methylation_levels_DMRfind(output_prefix+"_rms_results_collapsed.tsv",
                                   output_prefix+"_rms_results_collapsed_with_levels.tsv",
                                   allc_files,
                                   samples,
                                   mc_type=mc_type,
                                   num_procs=num_procs,
                                   buffer_line_number=buffer_line_number)
    subprocess.check_call(shlex.split("mv "+output_prefix+"_rms_results_collapsed_with_levels.tsv "+
                                      output_prefix+"_rms_results_collapsed.tsv"))
    print_checkpoint("Done")

def filter_collapsed(filen,output,min_level_diff=0,min_DMS=0,samples = [],hyper_samples=[],hypo_samples=[],strict=False):
    """
    This function filters collapsed DMR files based on the parameters passed and writes to output.
    
    filen is the path to the collapsed file that was output by DMRfind
    output is the path that you want to write the results to
    min_level_diff is the minimum methylation level difference between max and min sites in a block
    min_DMS is the minimum number of sites that need to be in a block for it to be output
    samples is a list of samples for which you want DMRs. These samples can be hyper or hypo methylated.
    hyper_samples is a list of samples that are hypermethylated to be included in the output.
    hypo_samples is a list of samples that are hypomethylated to be included in the output.
    strict is a boolean that if set to True means that all entries in samples/hyper_samples/hypo_samples
        must be present. If set to False, only one sample is necessary as long as the other samples are 
        not in the opposite category (e.g., a sample in hypo_samples is in the hypermethylated list)
    """
    f = open(filen, 'r')
    g = open(output, 'w')
    
    for sample in hypo_samples:
        if sample in hyper_samples:
            exit(sample+" found in both hyper and hypo sample list. These lists must not overlap.")
    for sample in samples:
        if sample in hyper_samples or sample in hypo_samples:
            exit(sample+" found in the hyper or hypo list. The samples list must not overlap with either the hyper or hypo list.")
    
    g.write(f.readline()) #write header in new file
    for line in f:
        write = True
        fields = line.rstrip().split('\t')
        
        levels = [float(x) for x in fields[6:] if x!="NA"]
        if int(fields[3]) < min_DMS or max(levels) - min(levels) < min_level_diff:
            continue
        else:
            hyper_count = 0
            hypo_count = 0
            sample_count = 0
            
            hyper = fields[4].split(',')
            hypo = fields[5].split(',')                
            
            for sample in samples:
                if sample in hyper:
                    hyper_count += 1
                if sample in hypo:
                    hypo_count += 1
            
            for sample in hyper_samples:
                if sample in hypo:
                    write = False
                    break
                if sample in hyper:
                    hyper_count += 1
            
            for sample in hypo_samples:
                if sample in hyper:
                    write = False
                    break
                if sample in hypo:
                    hypo_count += 1
            
            if strict == True and (hyper_count + hypo_count) != (len(hyper_samples) + len(hypo_samples) + len(samples)):
                write = False
            if strict == False and (hypo_count == 0 or hyper_count == 0):
                write = False
                    
        if write:
            g.write(line)
        
def run_rms_tests(files,output,samples,min_cov = 0,num_sims="10000",num_sig_tests="100",seed=-1,keep_temp_files=False):
    """
    This function allows for the parallelization of the C++ code that runs the root mean square tests
    """
    c_exe_path = __file__[:__file__.rfind("/")]+"/run_rms_tests.out"
    check_call(split(" ".join([c_exe_path,",".join(files),output,",".join(samples),str(min_cov),str(num_sims),str(num_sig_tests),str(seed)])))
    try:
        with open(output, 'r'):
            check_call(split("sort -k 2n,2n -o "+output+" "+output))
    except IOError:
        print "No results in "+output
    if keep_temp_files == False:
        check_call(split("rm "+" ".join(files)))

def histogram_correction_DMRfind(rms_results,num_sims,num_sig_tests,target_fdr =.01,max_iterations=1000,convergence_diff=.01):
    """
    This function uses the method outlined here:
    http://www.stat.iastate.edu/preprint/articles/2009-04.pdf
    to calculate the p-value cutoff for a particular FDR.
    
    rms_results is the path to a file generated by the DMRfind function. In other words, the results for the root mean square
        tests at each site
    num_sims is the maximum number of permutations the root mean square tests were run with
    num_sig_tests is an integer indicating how many permutations can return a result more significant than the original statistic
        before permutation testing is abandonded for that particular site.
    target_fdr is the FDR for which you wish to calculate a p-value cutoff. Note, this is a target because the exact FDR you desire
        may not be possible. This function will give you the p-value cutoff for the largest FDR that is less than or equal to your target FDR.
    max_iterations is the maximum number of iterations performed by the algorithm described in the paper cited above
    convergence_diff determines when the algorithm will terminate. When the current m0 estimate and the last m0 estimate differ by no
        more than convergence_diff.
    """
    #The minus 1 is for the header
    total_tests=int(check_output(split("wc -l "+rms_results)).split(" ")[0])-1
    #A dict of lists
    #the list is in the format [expected_count, observed_count,total_significant] the key is the pvalue
    table = {(1,num_sims):[(1.0/num_sims)*total_tests,0]}
    precision = int(ceil(log10(num_sims)))
    last_pvalue = 0
    sorted_pvalues = [(1,num_sims)]
    for numerator in xrange(2,num_sig_tests+1):
        pvalue = round(float(numerator) / num_sims,precision)
        expected_value=(pvalue-last_pvalue)*total_tests
        table[(numerator,num_sims)]=[expected_value,0]
        if pvalue not in sorted_pvalues:
            sorted_pvalues.append((numerator,num_sims))
        last_pvalue = pvalue
    for denominator in reversed(xrange(num_sig_tests,num_sims)):
        pvalue = round(float(num_sig_tests) / denominator,precision)
        expected_value=(pvalue-last_pvalue)*total_tests
        table[(num_sig_tests,denominator)]=[expected_value,0]
        if pvalue not in sorted_pvalues:
            sorted_pvalues.append((num_sig_tests,denominator))
        last_pvalue = pvalue
    f = open(rms_results,'r')
    f.readline()
    for line in f:
        line = line.rstrip()
        fields = line.split("\t")
        table[(int(fields[-2]),int(fields[-1]))][1]+=1
    f.close()
    
    last_m0_estim = 0
    m0_estim = total_tests
    m0_estim_diff = 10000
    iter = 0
    while iter < max_iterations and m0_estim_diff > convergence_diff:
        print "m0 estimate for iteration "+str(iter)+": "+str(m0_estim)
        difference = 0
        for pvalue in sorted_pvalues:
            if table[pvalue][1] < table[pvalue][0]:
                break
            difference += table[pvalue][1] - table[pvalue][0]
        last_m0_estim = m0_estim
        m0_estim = total_tests - difference
        m0_estim_diff = abs(float(last_m0_estim)-m0_estim)
        iter+=1
        #update expected pvalue
        last_pvalue = 0
        for pvalue in sorted_pvalues:
            pvalue_frac = round(float(pvalue[0]) / pvalue[1],precision)
            expected_value=(pvalue_frac-last_pvalue)*m0_estim
            table[pvalue][0]=expected_value
            last_pvalue = pvalue_frac
    if iter == 1 or (iter == max_iterations and m0_estim_diff > convergence_diff):
        print "Histogram FDR correction did not converge. Switching to Benjamini-Hochberg."
        return benjamini_hochberg_correction_DMRfind(rms_results,target_fdr)
    print "m0 estimate for iteration "+str(iter)+": "+str(m0_estim)
    #get the number of tests that would be declared significant at each p-value
    #and figure out the FDR
    
    running_count = table[sorted_pvalues[0]][1]
    diff = 1

    frac = round(float(sorted_pvalues[0][0]) / sorted_pvalues[0][1],precision)
    pvalue_cutoff = sorted_pvalues[0]
    best_fdr = frac * m0_estim / running_count
    for pvalue in sorted_pvalues[1:]:
        running_count+= table[pvalue][1]
        frac = round(float(pvalue[0]) / pvalue[1],precision)
        fdr = frac * m0_estim / running_count
        #Require that the best FDR is less than or equal to the target FDR
        if abs(target_fdr - fdr) < diff and fdr <= target_fdr:
            best_fdr = fdr
            pvalue_cutoff = pvalue
            diff = abs(target_fdr - fdr)
    print "Difference between best and last m0 estimate: "+str(m0_estim_diff)
    print "The closest p-value cutoff for your desired FDR is "+str(float(pvalue_cutoff[0])/pvalue_cutoff[1])+" which corresponds to an FDR of "+str(best_fdr)
    return float(pvalue_cutoff[0])/pvalue_cutoff[1]

def get_resid_cutoff(resid_cutoff, pvalue_cutoff, num_samples, rms_file):
    """
    Get all residual values from rms_file that are above the p-value cutoff and return the
    resid_cutoff percentile number.
    
    resid_cutoff - fraction that you want the cutoff to be from the distribution (.01 -> 1% extreme is cutoff)
    pvalue_cutoff - includes only entries that are above pvalue cutoff to find the resid cutoff
    num_samples - number of samples
    rms_file
    """
    index = 5 + num_samples * 3 #get the starting index for the residual fields
    resid_cutoff = 100 - resid_cutoff*100
    residuals = []
    f = open(rms_file, 'r')
    f.readline()
    for line in f:
        #All of the negative results should be equivalent so only
        #grab 10M of them.
        if len(residuals) > 10000000:
            break
        line = line.split("\t")
        if float(line[4]) > pvalue_cutoff:
            for resid in line[index:index+num_samples*2]:
                residuals.append(float(resid))
    f.close()
    if len(residuals) == 0:
        print "There are no null residuals to calculate resid_cutoff. Using 2 as the cutoff."
        return 2
    return scoreatpercentile(residuals, resid_cutoff)
    

def benjamini_hochberg_correction_DMRfind(filen,sig_cutoff):
    """
    This function is similar to the one defined here:
    http://stats.stackexchange.com/questions/870/multiple-hypothesis-testing-correction-with-benjamini-hochberg-p-values-or-q-va
    But takes advantage of the fact that the elements provided to it are in a sorted file.
    This way, it doesn't have to load much into memory.
    This link:
    http://brainder.org/2011/09/05/fdr-corrected-fdr-adjusted-p-values/
    was also helpful as the monotonicity correction from stats.stackexchange is not correct.
    
    file is a string indicating the path to an rms_results file (generated by DMRfind).
    sig_cutoff is the FDR cutoff you'd like to use to indicate if a site is significant.
    """
    total_tests=int(check_output(split("wc -l "+filen)).split(" ")[0])
    check_call(split("sort -k 5g,5g -o "+filen+" "+filen))
    f = open(filen,'r')
    line = f.readline()
    test_num = 1
    prev_bh_value = 0
    best_fdr = 0
    best_pvalue = 0
    for line in f:
        fields = line.rstrip().split("\t")
        bh_value = float(fields[4]) * total_tests / (test_num + 1)
        # Sometimes this correction can give values greater than 1,
        # so we set those values at 1
        bh_value = min(bh_value, 1)
        prev_bh_value = bh_value
        if bh_value <= sig_cutoff and bh_value >= best_fdr:
            best_fdr = bh_value
            best_pvalue = fields[4]
            
        test_num += 1
    f.close()
    print "The closest p-value cutoff for your desired FDR is "+str(best_pvalue)+" which corresponds to an FDR of "+str(best_fdr)
    print_checkpoint("Sorting "+filen+" by position")
    check_call(split("sort -k 1n,1n -k 2n,2n "+filen+" -o "+filen))
    return float(best_pvalue)

def check_clusters(category_dict, min_cluster, block):
    """
    Check that the given block has the minimum amount of clustered category. 
    category_dict has a key of sample that points to a category number.
    min_cluster is the minimum amount of a category needed in a list to be valid.
    block is a list for a certain block sent from collapse_dmr_windows
    
    If any of the categories is valid, return true. Otherwise return false.
    """  
    #This is to ensure that if a sample category has fewer than min_cluster samples
    #in total (regardless of hyper/hypo status) that it won't get left out.
    total_category = [0] * len(set(category_dict.values()))
    for category_number in category_dict.values():
        total_category[category_number] += 1
    #create a list of counters for hypermethylated and hypomethylated sites

    meth_count = [0] * len(set(category_dict.values()))
    unmeth_count = [0] * len(set(category_dict.values()))
    
    #split lists into samples
    meth_list = block[4].split(",")
    unmeth_list = block[5].split(",")
    
    #check if there are no samples in a certain list
    if filter(None, meth_list):
        for sample in meth_list: 
            meth_count[category_dict[sample]] += 1 #increment the counter
    if filter(None, unmeth_list):
        for sample in unmeth_list:
            unmeth_count[category_dict[sample]] += 1
    
    #check if at least two sample groups have the minimum number of samples in the unmethylated
    #and methylated categories
    if [count for index,count in enumerate(meth_count) if count >= min_cluster or (total_category[index] < min_cluster and total_category[index] == count)] and [count for index,count in enumerate(unmeth_count) if count >= min_cluster or (total_category[index] < min_cluster and total_category[index] == count)]:
        return True
        
    return False

def collapse_dmr_windows(inputf, output, column=4, max_dist=100, resid_cutoff=False, sig_cutoff=.01, min_num_dms=0,
                         collapse_samples=False, sample_category=False, min_cluster=0): 
    """
    input - name of tab separated input file
    output - name of file for collapsed output
    column - the field (STARTING NUMBERED FROM ZERO!) that contains the pvalue from the rms test.
    max_dist - the maximum number of nucleotides two sites with 1s in the "column" field can be and be combined into a block
    resid_cutoff - If this option is specified, not only will a result have to be significant to be included in a window, but it will also,
        have to show deviations in the contingency table in the same direction as the rest of the window. For example, if sample A is generally
        showing a degree of methylation higher than expected and sample B lower than expected, any new site will have to have these same properties.
        Furthermore, these deviations have to be at least as extreme as resid_cutoff (in the positive or negative direction). Two is a good value to start
        with here.
    sig_cutoff sets the threshold for considering sites to be combined. All sites with a p-value greater than this cutoff will be excluded.
    min_num_dms -
    collapse_samples - Base the sample filtration of hypermethylated and hypomethylated sites only if they are on this samples list
    sample_category - A list of categories that each respective sample belongs to; the categories must begin at 0 and increase by
        1 for each category added. ex: samples [A,B,C] categories [0,1,2] or categories [0, 1, 0] 
    min_cluster - The minimum number of each sample category that must be present in every block that is output.
    """
    #Argument checks  
    try:
        resid_cutoff = float(resid_cutoff)
    except:
        exit("resid_cutoff must be a float or something that can be cast as a float")
    try:
        sig_cutoff = float(sig_cutoff)
    except:
        exit("sig_cutoff must be a float or something that can be cast as a float")
    try:
        max_dist = int(max_dist)
    except:
        exit("max_dist must be an integer or something that can be cast as an integer")
    
    if collapse_samples and sample_category:
        category_dict = dict(zip(collapse_samples, sample_category))
    elif sample_category:
        exit("In order to use sample_category, you must specify a corresponding list of samples in collapse_samples!")
    
    #this is where sample specific fields begin. Since the number of fields before these fields is fixed
    #I created a variable in case I add a fixed field later.
    fields_offset = 5  
    #Collapse DMRs into windows
    f = open(inputf,'r')
    g = open(output,'w')
    g.write("chr\tstart\tend\tnumber_of_dms\thypermethylated_samples\thypomethylated_samples\n")
    line = f.readline()
    fields = line.split("\t")
    num_samples = len(fields[fields_offset:])/5
    samples = [mc_field[3:] for mc_field in fields[fields_offset:fields_offset+num_samples]]
    
    line = f.readline().rstrip()
    fields = line.split("\t")
    #Number of significant blocks
    block_count = 0
    #some fields may contain a less than sign because of the permutation procedure used to
    #generate the p-values. These sites are automatically assumed to be significant
    while line and (fields[column].find("<")== -1 and float(fields[column]) > sig_cutoff):
        line = f.readline()
        fields = line.split("\t")
    #Edge case, need to be sure we haven't just read through the whole file (i.e., there are no
    #significant results)
    if line:
        #chromosome, start, end, number of mCs, hypermethylated samples, hypomethylated samples
        block = [fields[0],fields[1],fields[1],1,[],[]]
        line = f.readline()
        line = line.rstrip()
    while line:
        fields = line.split("\t")
        try:
            #some fields may contain a less than sign because of the permutation procedure used to
            #generate the p-values. These sites are automatically assumed to be significant
            if fields[column].find("<")!= -1 or float(fields[column]) <= sig_cutoff:
                #Check that block is within a certain distance and that the methylation level order is the same
                if block[0] == fields[0] and int(block[2]) + max_dist >= int(fields[1]):
                    #check that either the hypermethylated residual is contributing significantly and positively to the chisquare or the hypomethylated residual is contributing
                    #significantly and negatively. This checks not only significance but the direction of the significance as well
                    hypermethylated_samples = [sample for index,sample in enumerate(samples) if (collapse_samples==False or sample in collapse_samples) and\
                         fields[fields_offset+index+num_samples*3] != "NA" and fields[fields_offset+index+num_samples*4] != "NA" \
                         and ((float(fields[fields_offset+index+num_samples*3]) >= resid_cutoff and float(fields[fields_offset+index+num_samples*4]) < resid_cutoff)\
                          or (-1 * float(fields[fields_offset+index+num_samples*4]) >= resid_cutoff and -1 * float(fields[fields_offset+index+num_samples*3]) < resid_cutoff))
                         ]
                    hypomethylated_samples = [sample for index,sample in enumerate(samples) if (collapse_samples==False or sample in collapse_samples) and\
                         fields[fields_offset+index+num_samples*3] != "NA" and fields[fields_offset+index+num_samples*4] != "NA" \
                         and ((-1 * float(fields[fields_offset+index+num_samples*3]) >= resid_cutoff and -1 * float(fields[fields_offset+index+num_samples*4]) < resid_cutoff)\
                          or (float(fields[fields_offset+index+num_samples*4]) >= resid_cutoff and float(fields[fields_offset+index+num_samples*3]) < resid_cutoff))
                         ]                    
                    #If one of the hypermethylated/hypomethylated samples is in the wrong list, the list comrehension will return a non-empty list
                    #which has a bool value of True. This will be not'ed to False and fail the if statement
                    if resid_cutoff != False and not([i for i in hypomethylated_samples if (collapse_samples == False or i in collapse_samples) and i in block[4]]) and not([i for i in hypermethylated_samples if (collapse_samples == False or i in collapse_samples) and i in block[5]]):
                        block[2] = fields[1]
                        block[3]+=1
                        
                        for sample in samples:
                            if sample in hypermethylated_samples and sample not in block[4]:
                                block[4].append(sample)                            
                            if sample in hypomethylated_samples and sample not in block[5]:
                                block[5].append(sample)
                    elif resid_cutoff == False:
                        block[2] = fields[1]
                        block[3]+=1
                    else:
                        #check to make sure all the samples aren't shoved into one category
                        #no sense in reporting that. Could happen if collapse_samples is being used
                        if (collapse_samples and len(block[4]) != len(collapse_samples) and len(block[5]) != len(collapse_samples)) or (not collapse_samples and len(block[4]) != num_samples and len(block[5]) != num_samples):
                            block[4] = ",".join(block[4])
                            block[5] = ",".join(block[5])
                            if block[3] >= min_num_dms:
                                if not sample_category or check_clusters(category_dict, min_cluster, block):
                                    g.write("\t".join(map(str,block))+"\n")
                                    block_count+=1                                        
                        hypermethylated_samples = [sample for index,sample in enumerate(samples) if (collapse_samples==False or sample in collapse_samples) and\
                             fields[fields_offset+index+num_samples*3] != "NA" and fields[fields_offset+index+num_samples*4] != "NA" \
                             and ((float(fields[fields_offset+index+num_samples*3]) >= resid_cutoff and float(fields[fields_offset+index+num_samples*4]) < resid_cutoff)\
                              or (-1 * float(fields[fields_offset+index+num_samples*4]) >= resid_cutoff and -1 * float(fields[fields_offset+index+num_samples*3]) < resid_cutoff))
                             ]
                        hypomethylated_samples = [sample for index,sample in enumerate(samples) if (collapse_samples==False or sample in collapse_samples) and\
                             fields[fields_offset+index+num_samples*3] != "NA" and fields[fields_offset+index+num_samples*4] != "NA" \
                             and ((-1 * float(fields[fields_offset+index+num_samples*3]) >= resid_cutoff and -1 * float(fields[fields_offset+index+num_samples*4]) < resid_cutoff)\
                              or (float(fields[fields_offset+index+num_samples*4]) >= resid_cutoff and float(fields[fields_offset+index+num_samples*3]) < resid_cutoff))
                             ] 
                        block = [fields[0],fields[1],fields[1],1,hypermethylated_samples,hypomethylated_samples]
                else:
                    #check to make sure all the samples aren't shoved into one category
                    #no sense in reporting that. Could happen if collapse_samples is being used
                    if (collapse_samples and len(block[4]) != len(collapse_samples) and len(block[5]) != len(collapse_samples)) or (not collapse_samples and len(block[4]) != num_samples and len(block[5]) != num_samples):
                        block[4] = ",".join(block[4])
                        block[5] = ",".join(block[5])
                        if block[3] >= min_num_dms:
                            if not sample_category or check_clusters(category_dict, min_cluster, block):                                                    
                                g.write("\t".join(map(str,block))+"\n")
                                block_count+=1                                    
                    hypermethylated_samples = [sample for index,sample in enumerate(samples) if (collapse_samples==False or sample in collapse_samples) and\
                         fields[fields_offset+index+num_samples*3] != "NA" and fields[fields_offset+index+num_samples*4] != "NA"\
                         and ((float(fields[fields_offset+index+num_samples*3]) >= resid_cutoff and float(fields[fields_offset+index+num_samples*4]) < resid_cutoff)\
                          or (-1 * float(fields[fields_offset+index+num_samples*4]) >= resid_cutoff and -1 * float(fields[fields_offset+index+num_samples*3]) < resid_cutoff))
                         ]
                    hypomethylated_samples = [sample for index,sample in enumerate(samples) if (collapse_samples==False or sample in collapse_samples) and\
                         fields[fields_offset+index+num_samples*3] != "NA" and fields[fields_offset+index+num_samples*4] != "NA" \
                         and ((-1 * float(fields[fields_offset+index+num_samples*3]) >= resid_cutoff and -1 * float(fields[fields_offset+index+num_samples*4]) < resid_cutoff)\
                          or (float(fields[fields_offset+index+num_samples*4]) >= resid_cutoff and float(fields[fields_offset+index+num_samples*3]) < resid_cutoff))
                         ] 
                    block = [fields[0],fields[1],fields[1],1,hypermethylated_samples,hypomethylated_samples]               
                     
            line = f.readline()
            line = line.rstrip()
        except Exception, e:
            exc_type, exc_obj, exc_tb = exc_info()
            print(exc_type, exc_tb.tb_lineno)
            print e
            set_trace()
        
    f.close()
    #Write out the last block
    if (block) and ((collapse_samples and len(block[4]) != len(collapse_samples) and len(block[5]) != len(collapse_samples)) or (not collapse_samples and len(block[4]) != num_samples and len(block[5]) != num_samples)):
        block[4] = ",".join(block[4])
        block[5] = ",".join(block[5])
        if block[3] >= min_num_dms:
            if not sample_category or check_clusters(category_dict, min_cluster, block):                                                    
                g.write("\t".join(map(str,block))+"\n")
                block_count+=1
    g.close()
    return block_count


def get_methylation_levels_DMRfind(input_tsv_file,
                                   output,
                                   input_allc_files,
                                   samples,
                                   mc_type=["CGN"],
                                   num_procs=1,
                                   buffer_line_number=100000):
    """
    This function assumes that allc files are of the format allc_<sample>_<chr>.tsv
    input is the path to a file containing collapsed DMR results
    output is the path to a file where the methylation values should be stored
    samples is a list of samples you'd like to compute the methylation level for
    path_to_allc is the path to the directory containing the allc files for these samples
    num_procs is an integer indicating the number of processors you'd like to use for calculating
        methylation level. This function can be parallelized up to the number of samples
    """
    # need to check
    # 1. whether the tsv file is sorted by coordinate
    # 1. check whether allc files are stored in a list
    # 2. check whether samples are stored in a list
    # 4. whether samples and allc files have the same length
    
    mc_class = expand_nucleotide_code(mc_type)
    
    with open(output,'w') as g:
        # header
        f = open(input_tsv_file,'r')
        line = f.readline()
        line = line.rstrip("\n")
        fields = line.split("\t")
        g.write("\t".join(fields)+"\t"+
                "\t".join(["methylation_level_"+sample for sample in samples])+"\n")
        # 
        methylation_levels = {}
        if num_procs == 1:
            for allc_file,sample in zip(input_allc_files,samples):
                get_methylation_level_DMRfind_worker(
                    input_tsv_file,
                    allc_file,
                    sample,
                    output,
                    mc_class,
                    buffer_line_number)
        else:
            pool = Pool(min(num_procs,len(samples)))
            results = {}
            #for allc_file,sample in zip(input_allc_files,samples):
            for ind in xrange(len(samples)):
                pool.apply_async(
                    get_methylation_level_DMRfind_worker,
                    (input_tsv_file,
                     input_allc_files[ind],
                     samples[ind],
                     output,
                     mc_class,
                     buffer_line_number)
                )
            pool.close()
            pool.join()
        temp_files = {}
        for sample in samples:
            temp_files[sample]=open(output.replace(".tsv","")+"_"+sample+"_temp_methylation_levels.tsv",'r')

        f.seek(0)
        line = f.readline()
        for line in f:
            g.write(line.rstrip("\n"))
            for sample in samples:
                g.write("\t"+temp_files[sample].readline().rstrip("\n"))
            g.write("\n")
        for sample in samples:
            temp_files[sample].close()
            subprocess.check_call(shlex.split("rm "+output.replace(".tsv","")+"_"+sample+"_temp_methylation_levels.tsv"))

def get_methylation_level_DMRfind_worker(inputf_tsv,
                                         inputf_allc,
                                         sample,
                                         output,
                                         mc_class,
                                         buffer_line_number=100000):
    # open allc file
    if inputf_allc[-3:] == ".gz":
        allc_file = gzip.open(inputf_allc,'r')
    elif inputf_allc[-4:] == ".bz2":
        allc_file = bz2.BZ2File(inputf_allc,'r')
    else:
        allc_file = open(inputf_allc,'r')
    # scan allc file to set up a table for fast look-up of lines belong
    # to different chromosomes
    chrom_pointer = {}
    cur_chrom = ""
    cur_pointer = 0
    while True:
        line = allc_file.readline()
        if not line: break
        fields = line.split("\t")
        if fields[0] != cur_chrom:
            chrom_pointer[fields[0]] = cur_pointer
            cur_chrom = fields[0]
        cur_pointer = allc_file.tell()

    # init
    prev_chrom = ""
    prev_end = ""

    # get methylation level for each DMR
    out = ""
    line_counts = 0
    with open(inputf_tsv,'r') as f,open(output.replace(".tsv","")+"_"+sample+
                                        "_temp_methylation_levels.tsv",'w') as g:
        f.readline() # skip header
        for line in f:
            line = line.rstrip("\n")
            fields = line.split("\t")
            dmr_chr=fields[0]
            dmr_start = int(fields[1])
            dmr_end = int(fields[2])
            
            # get to new chromosome and redirect pointer to related lines in allc file
            if prev_chrom != fields[0]:
                allc_file.seek(chrom_pointer[dmr_chr])
                allc_prevbyte = chrom_pointer[dmr_chr]
                # update allc line
                allc_line=allc_file.readline()
                allc_field=allc_line.split("\t")


            #If this new dmr overlaps with the previous, begin where the previous start was found
            elif prev_end and dmr_start < prev_end:
                allc_file.seek(allc_prevbyte)
                allc_line = allc_file.readline()
                allc_field = allc_line.split("\t")
                
            #read up to the beginning of the dmr
            byte = allc_prevbyte #in case the new dmr never enters the loop, keep byte the same as previously
            while allc_line and int(allc_field[1]) < dmr_start:
                byte = allc_file.tell()
                allc_line=allc_file.readline()
                allc_field=allc_line.split("\t")
            allc_prevbyte = byte #record the byte where dmr_start was found
            
            mc = 0
            h = 0
            while allc_line and int(allc_field[1]) >= dmr_start and int(allc_field[1]) <= dmr_end:
                if allc_field[3] in mc_class:
                    mc += int(allc_field[4])
                    h += int(allc_field[5])
                allc_line=allc_file.readline()
                allc_field=allc_line.split("\t")
            if h != 0:
                methylation_level = str(float(mc) / h)
            else:
                methylation_level = "NA"

            out += str(methylation_level)+"\n"
            line_counts += 1
            if line_counts > buffer_line_number:
                g.write(out)
                line_counts = 0
                out = ""
                
            prev_chrom = dmr_chr
            prev_end = dmr_end

        if line_counts > 0:
            g.write(out)
            line_counts = 0
