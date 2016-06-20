from simulation import make_allc
from DMRfind import DMRfind
import subprocess
import shlex
import sys
import matplotlib.pyplot as plt
import os
import pybedtools as pbt

def main(path, procs):
    """
    Below is how to run for DMRfind, for others use the same syntax, except call run_misc and then call
    graph() with a different output string
    """
    path += "/"
    errors = [0.05]
    coverages = [2,5,10,20,30]
    
    #DMRfind options:
    interval = True
    mc_max_dist = [0]
    DMRfind_options = [y for y in mc_max_dist] #every combination of options
    DMRfind_output = path + "graph_simulation_results_DMRfind"
    
    #BSmooth options:
    Bsmooth_output = path + "graph_simulation_results_Bsmooth"
    ns = [7]
    h = [1000]
    Bsmooth_options = [(x,y) for x in ns for y in h]

    #methylkit options:
    methylkit_output = path + "graph_simulation_results_methylkit"
    
    DMRfind_results = [] #list of (label for graph, results from run())
    Bsmooth_results = []
    methylkit_results = []
     
    for error in errors:
        print "Running for error " + str(error)
        names = simulate_allc_files(coverages, path, error) #simulate allc files
        samples = correct_name(path, names, coverages, replicates=True) #adjust file names for DMRfind
        correct_bed_files = summary_to_bed(names, coverages)
        
        #DMRfind:
        for dist in DMRfind_options:
            print "Running mc_max_dist: " + str(dist)
            DMRfind_results.append(("Error: " + str(error) + ", Dist: " + str(dist), 
                            run(correct_bed_files, coverages, path, names, samples, mc_max_dist=dist, 
                                interval=interval, procs=procs)))
        
        #Bsmooth:
        for x,y in Bsmooth_options:
            print "Running Bsmooth with ns:" + str(x) + ", h:" + str(y)
            Bsmooth_results.append(("Error: " + str(error) + ", ns:" + str(x) + ", h:" + str(y), 
                            run_misc(correct_bed_files, coverages, path, names, samples, procs, "Bsmooth",ns=x, h=y)))
        
        #methylkit            
        methylkit_results.append(("Error: " + str(error), 
                        run_misc(correct_bed_files, coverages, path, names, samples, procs, "methylkit")))
        clean(path, names) #Clean other files from other pipelines?      
            
    graph(coverages, DMRfind_results, DMRfind_output, interval)
    graph(coverages, Bsmooth_results, Bsmooth_output, False)
    graph(coverages, methylkit_results, methylkit_output, False)

def run(correct_bed_files, coverages, path, names, samples, mc_max_dist=0, interval=False, procs=1):
    """
    Run with a certain mc_max_dist
    """    
    cutoffs = run_DMRfind(path, samples, coverages, mc_max_dist, procs)
    
    if interval:
        bed_files = zip(correct_bed_files, collapsed_to_bed(path, coverages))
        performance = judge_interval(bed_files)
    else:
        performance = judge(path, coverages, names, cutoffs)
        
    print_per(coverages, performance)
    return performance

def run_misc(correct_bed_files, coverages, path, names, samples, procs, program, ns="", h=""):
    """
    Run a miscellaneous function and then return the performance results
        coverages - list of coverages
        path - same as path_to_allc
        error - error for simulated reads
        procs - number of processors to run
        func - some function used to call your pipeline
        args, kwargs - arguments to the function
    
    Equivalent of run()
    """
    if program == "Bsmooth":
        bed_files = run_Bsmooth_template(path, samples, coverages, procs, ns, h) #call your function here and return list of bedfile names
    elif program == "methylkit":
        bed_files = run_methylkit_template(path, samples, coverages, procs) #call your function here and return list of bedfile names
    bed_files = zip(correct_bed_files, bed_files) #zip up correct files with the files returned by pipeline
    performance = judge_interval(bed_files) #pass in tuple of correct/actual bedfiles to get performance results
    
    print_per(coverages, performance)
    return performance

def run_Bsmooth_template(path, samples, coverages, num_procs, ns, h, 
                         script_path="/home/mschultz/utilities/my_python_modules/methylpy-paper/bsmooth_dmrfind.py",
                         Rscript_path="/home/mschultz/utilities/R/R-2.15.3/bin/Rscript"):
    """
    RUN SCRIPT HERE; equivalent of run_DMRfind()
        path ~ path_to_allc
        samples ~ samples
    """
    script_path = "/home/matteo/methylpy-paper/bsmooth_dmrfind.py"
    mc_type = "CGA"
    
    bed_files = []
    for cov,sample_list in zip(coverages,samples):
        print "Running Bsmooth for " + str(cov) + "x coverage"
        group1 = sample_list[0] + "," + sample_list[0] + "rep"
        group2 = sample_list[1] + "," + sample_list[1] + "rep"
        output_prefix = path+"cov_"+str(cov)
        subprocess.check_call(["python2.7", script_path, mc_type, path, num_procs, "1", group1, group2, output_prefix, 
                               Rscript_path, str(ns), str(h)])
        bed_files.append(output_prefix+"_ns_"+str(ns)+"_h_"+str(h)+"_"+"bsmooth_dmrs.bed")
    return bed_files

def run_methylkit_template(path, samples, coverages, num_procs,
                           script_path="/home/mschultz/utilities/my_python_modules/methylpy-paper/methylkit_dmrfind.py",
                           Rscript_path="/home/mschultz/utilities/R/R-2.15.3/bin/Rscript"):
    """
    RUN SCRIPT HERE; equivalent of run_DMRfind()
        path ~ path_to_allc
        samples ~ samples
    """
    script_path = "/home/matteo/methylpy-paper/methylkit_dmrfind.py"
    mc_type = "CG"
    
    bed_files = []
    for cov,sample_list in zip(coverages,samples):
        print "Running methylkit for " + str(cov) + "x coverage"
        group1 = sample_list[0] + "," + sample_list[0] + "rep"
        group2 = sample_list[1] + "," + sample_list[1] + "rep"
        output_prefix = path+"cov_"+str(cov)
        subprocess.check_call(["python2.7", script_path, mc_type, path, num_procs, "1", group1, group2, output_prefix,
                               Rscript_path])
        collapse_methylkit(output_prefix+"_"+"methylkit_dmrs.bed", 100)
        bed_files.append(output_prefix+"_"+"methylkit_dmrs.bed")
    return bed_files

def simulate_allc_files(coverages, path, error):
    """
    Make a run simulation call for each coverage and get the file names in names array
    """
    names = [] #list of filenames that will be created by simulation
    for cov in coverages:
        print "Simulating " + str(cov) + "x coverage"
        names.append(run_simulation(cov, path, error=error, nonDMR_methyl=1-error))
    return names

def run_simulation(coverage, path, context="CGA", chr="1", strand="+", dmr_size=200, between_size=1000, dmr_diff=.5, num_dmr=1000,
                   error=.01, nonDMR_methyl=.7, num_between_pos=10, num_pos=8):
    """
    Run the simulation and make allc files
    """
    output = path + "grph" + str(coverage) + "output"
    make_allc(output, context, chr, strand, dmr_size, between_size, dmr_diff, num_dmr, error, coverage, nonDMR_methyl, 
              num_between_pos, num_pos=num_pos, num_files=2, balanced_mcs = True, replicates = True)
    return output

def collapse_methylkit(fname, dist):
    """
    Collapse a (methylkit) bed file into intervals based on positions within dist
    """
    f = open(fname, 'r')
    g = open(fname+"_collapsed", 'w')
    
    line = f.readline().rstrip().split('\t')
    prev_pos = int(line[1])
    start = int(line[1])
    chrom = line[0]
    for line in f:
        line = line.rstrip().split('\t')
        pos = int(line[1])
        if pos - prev_pos > dist:
            g.write('\t'.join([chrom, str(start), str(prev_pos)])+'\n')
            start = pos
            prev_pos = pos
        else:
            prev_pos = pos
    g.write('\t'.join([chrom, str(start), str(prev_pos)])+'\n')
    f.close()
    g.close()
    subprocess.check_call(['mv', fname+"_collapsed", fname])

def correct_name(path, names, coverages, replicates = False):
    """
    Adjust the simulation file names to match the DMRfind pattern:
        [ allc_sample_chr.tsv ]
        
    Each sample pair will be named: 
        [ A + coverage number, B + coverage number ]
    """
    samples = []
    for cov,name in zip(coverages,names):
        subprocess.check_call(["mv", name+"0", path + "allc_"+"A"+str(cov)+"_1.tsv"])
        subprocess.check_call(["mv", name+"1", path + "allc_"+"B"+str(cov)+"_1.tsv"])
        if replicates:
            subprocess.check_call(["mv", name+"0_replicate", path + "allc_"+"A"+str(cov)+"rep_1.tsv"])
            subprocess.check_call(["mv", name+"1_replicate", path + "allc_"+"B"+str(cov)+"rep_1.tsv"])
        samples.append(["A"+str(cov), "B"+str(cov)])
    return samples

def run_DMRfind(path, samples, coverages, mc_max_dist, procs):
    """
    Run DMRfind for each coverage pair and rename the results to:
        [ collapsed + coverage number ], [ results + coverage number ]
    """
    cutoffs = []
    dmrfindoutput = path + 'dmrfindpvalueoutput'
    for cov,sample_list in zip(coverages,samples):
        print "Running DMRfind for " + str(cov) + "x coverage"
        stdout = sys.stdout
        sys.stdout = open(dmrfindoutput, 'w')
        DMRfind("CG", {"chr1":[0,5000000]}, sample_list, path, num_procs=procs, 
                min_cov=0,mc_max_dist=mc_max_dist,dmr_max_dist=100,num_sims=1000,num_sig_tests=400, save_result=path+"temp")
        sys.stdout.close()
        sys.stdout = stdout
        subprocess.check_call(["mv", path + "temp_rms_results_collapsed.tsv", path + "collapsed" + str(cov)])
        subprocess.check_call(["mv", path + "temp_rms_results.tsv", path + "results" + str(cov)])
        cutoffs.append(get_pvalue_cutoff(dmrfindoutput))
        subprocess.check_call(["rm", dmrfindoutput])
    return cutoffs

def get_pvalue_cutoff(infile):
    """
    Parse the DMRfind output to get the pvalue cutoff
    """
    f = open(infile, 'r')
    for line in f:
        line = line.rstrip()
        if "The closest p-value cutoff for your desired FDR is" in line:
            line = line.split(' ')
            return line[9]

def judge(path, coverages, names, cutoffs):
    """
    Judge the correctness of DMRfind by looking through and finding whether the positions found were correct.
    
    Get the p-value cutoff from DMRfind, look through the results file and judge whether positions were below the cutoff.
    If the position was below the cutoff, then that site was judged as significant by DMRfind. 
    
    If a position is significant &:
        - exists in one of the simulated blocks -> true positive
        - not in one of the simulated blocks -> false positive
    
    If a position is insignificant &:
        - exists in one of the simulated blocks -> false negative 
        - not in one of the simulated blocks -> true negative
    """
    #actual files
    results = [path + "results" + str(x) for x in coverages]
    #simulated files
    sim_summaries = [str(x) + "_allc_sim_summary" for x in names]
    sim_positions = [path + "allc_A" + str(x) + "_1.tsv" for x in coverages]
    
    #get list of blocks from first simulated summary file (all simulations have same block coordinates)
    sim_blocks = []
    f = open(sim_summaries[0], 'r')
    f.readline()
    for line in f:
        line = line.rstrip().split('\t')
        sim_blocks.append((int(line[1]), int(line[2])))
    f.close()
    
    performance = [] #list of results from different coverages
    for result, sim_position, cutoff in zip(results, sim_positions, cutoffs):
        tp, fp, tn, fn = 0, 0, 0, 0
        sim_pos_dict = {} #dictionary of simulated position to 0 if insignificant or 1 if significant
        sim_block_index = 0
        fsim_positions = open(sim_position, 'r') #get simulated positions and find which were inside DMR blocks
        for simline in fsim_positions:
            simline = simline.rstrip().split('\t')
            simpos = int(simline[1])
            if simpos >= sim_blocks[sim_block_index][0]:
                if simpos <= sim_blocks[sim_block_index][1]: #was in between start and end of DMR block
                    sim_pos_dict[simpos] = 1 #inside DMR block
                elif sim_block_index + 1 >= len(sim_blocks):
                    sim_pos_dict[simpos] = 0 #very last edge of the file is outside DMR block
                elif simpos < sim_blocks[sim_block_index+1][0]:
                    sim_pos_dict[simpos] = 0 #outside current DMR block
                    sim_block_index += 1
            else:
                sim_pos_dict[simpos] = 0 #this position is before the next DMR block but after the previous
        fsim_positions.close()

        fresult = open(result, 'r')
        fresult.readline()
        for resline in fresult: #loop through each position found by dmrfind
            resline = resline.rstrip().split('\t')
            respos = int(resline[1])
            pvalue = float(resline[4])
            if pvalue <= float(cutoff): #significant by dmrfind
                if sim_pos_dict[respos] == 1: #also significant by simulation
                    tp += 1
                else: #true one was not significant
                    fp += 1
            else: #insignificant by dmrfind
                if sim_pos_dict[respos] == 1: #significant by simulation
                    fn += 1
                else: #was insignificant by simulation
                    tn += 1
            del sim_pos_dict[respos] #delete the positions from hash as they are found
        fresult.close()
        
        for val in sim_pos_dict.values():
            if val == 1: #simulated position was significant, not found by dmrfind
                fn += 1
            else: #simulated position was insignificant, thrown away by dmrfind
                tn += 1
        performance.append((tp,tp+fp,tn,tn+fn))
    return performance

def judge_interval(bed_files):
    """
    bed_files is a list of tuples with (correct_bed_filename, actual_bed_filename)
    
    Sensitivity = #true pos / #total pos
        # true pos -> number of bp in the overlaps
        # total pos -> number of bp in correct blocks
    
    Specificity = #true neg / #total neg
        #true neg -> invert both correct and actual blocks and count number of bp
        #total neg -> number of bp outside of correct blocks
    """
    performance = []
    for correct, actual in bed_files:
        try:
            correct = pbt.BedTool(correct)
            actual = pbt.BedTool(actual)
            intersect = correct.intersect(actual)
        except:
            print "DMRfind collapsed file is empty"
            performance.append((0,1,1,1))
            continue
        
        correct_neg = [] # 'negative' blocks from simulated data
        
        #get total pos and total neg
        total_pos, total_neg = 0, 0
        prev_end = 0
        for interval in correct:
            st = int(interval[1])
            end = int(interval[2])
            total_pos += end - st + 1
            total_neg += st - prev_end + 1
            correct_neg.append("chr1 " + str(prev_end) + " " + str(st))
            prev_end = end
        correct_neg = "\n".join(correct_neg)
        correct_neg = pbt.BedTool(correct_neg, from_string=True)

        actual_neg = [] # 'negative' blocks from real data (inverse)

        #get true positives, number of bp in overlaps
        true_pos = 0
        prev_end = 0
        for interval in intersect:
            st = int(interval[1])
            end = int(interval[2])
            true_pos += end - st + 1
            if st < prev_end:
                print "Warning, some blocks found with your method overlap!"
            actual_neg.append("chr1 " + str(prev_end) + " " + str(st))
            prev_end = end
        actual_neg = "\n".join(actual_neg)
        actual_neg = pbt.BedTool(actual_neg, from_string=True)
        
        #get true negatives by intersecting simulated negative and real negative blocks
        true_neg = 0
        for interval in correct_neg.intersect(actual_neg):
            st = int(interval[1])
            end = int(interval[2])
            true_neg += end - st + 1
        
        performance.append((true_pos, total_pos, true_neg, total_neg))
    
    return performance 
    
def summary_to_bed(names, coverages):
    """
    Convert simulated summary files to bed files
    """
    sim_summaries = [str(x) + "_allc_sim_summary" for x in names]
    for summary in sim_summaries:
        f = open(summary + ".bed", 'w')
        process1 = subprocess.Popen(["sed", "1d", summary], stdout=subprocess.PIPE)
        process2 = subprocess.Popen(["cut", "-f", "1-3"], stdin=process1.stdout, stdout=f)
        f.close()
    return [x + ".bed" for x in sim_summaries]

def collapsed_to_bed(path, coverages):
    """
    Convert collapsed files from DMRfind to bed files
    """    
    results = [path + "collapsed" + str(x) for x in coverages]
    for result in results:
        f = open(result + ".bed", 'w')
        process1 = subprocess.Popen(["sed", "1d", result], stdout=subprocess.PIPE)
        process2 = subprocess.Popen(["cut", "-f", "1-3"], stdin=process1.stdout, stdout=f)
        f.close()       
    return [x + ".bed" for x in results]
    
def print_per(c, per):
    """
    Print performance information
    """
    for p, cov in zip(per,c):
        print str(cov) + "x coverage: " + "(TP, TP+FP, TN, TN+FN): " + str(p)

def graph(coverages, results, output, interval):
    """
    Graph the performance results
    """
    #x-axis: coverages, [2, 5..., 30]
    #y-axis: performance values
    colors = ['b', 'r', 'g', 'k', 'c', 'm', 'y']
    titles = ["Sensitivity", "Specificity"]
    y_axis = [[], []]
    legend_titles = []
    for (param,option), color in zip(results, colors):
        y_axis[0].append((color,[x[0]/float(x[1]) for x in option]))
        y_axis[1].append((color,[x[2]/float(x[3]) for x in option]))
        legend_titles.append(param)
    
    plt.subplots_adjust(wspace=0.4, hspace=0.4)
    for sub, y, title in zip(range(211,222), y_axis, titles):
        plt.subplot(sub)
        plt.xlim(min(coverages), max(coverages))
        #y_min, y_max = min(y[0][1]), max(y[0][1])
        plt.xlabel("Coverage")
        if interval:
            plt.ylabel("Fraction of positions (interval)")
        else:
            plt.ylabel("Fraction of positions")
        plt.title(title)
        to_be_labeled = []
        for color, performance in y:
            to_be_labeled.append(plt.plot(coverages, performance, color, linewidth=2.0))
            #y_min = min(y_min, min(performance))
            #y_max = max(y_max, max(performance))
        #plt.ylim(max(y_min-10,0), y_max + 10)
    plt.figlegend(to_be_labeled, legend_titles, 'upper right', prop={'size':7})
    plt.subplots_adjust(bottom=0.1, right=0.8, top=0.9)
    plt.savefig(output + ".pdf")
    print "Graph result saved in " + output + ".pdf"

def clean(path, names):
    """
    Clean up files created by run()
    """
    try:
        subprocess.check_call("rm " + path + "allc_A*", shell=True)
        subprocess.check_call("rm " + path + "allc_B*", shell=True)
        for name in names:
            subprocess.check_call("rm " + name + "*", shell=True)
    except:
        pass
    try:
        subprocess.check_call("rm " + path + "collapsed*", shell=True)
        subprocess.check_call("rm " + path + "results*", shell=True)
    except:
        pass
    try:
        subprocess.check_call("rm " + path + "cov_*", shell=True)
    except:
        pass

if __name__ == '__main__':
    if len(sys.argv) > 2:
        path = sys.argv[1]
        procs = sys.argv[2]
    elif len(sys.argv) == 2:
        path = sys.argv[1]  
        procs = '1'      
    else:
        path = os.getcwd()
        procs = '1'
    main(path, procs)