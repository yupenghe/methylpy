'''
Created on May 9, 2012

@author: Matteo
'''
import sys
import random as rand
from subprocess import Popen, PIPE,check_call
import operator
import pdb
import warnings
import exceptions
try:
    import numpy as np
    #cimport numpy as np
except Exception,e:
    exc_type, exc_obj, exc_tb = sys.exc_info()
    print(exc_type, exc_tb.tb_lineno)
    print e
    sys.exit("methylpy.stats requires the numpy module")
try:
    import scipy.stats as sci
except Exception,e:
    exc_type, exc_obj, exc_tb = sys.exc_info()
    print(exc_type, exc_tb.tb_lineno)
    print e
    sys.exit("methylpy.stats requires the scipy module")
"""
def run_chisquare_tests(files,output,samples,min_cov = 0,buffer_limit=1):
    #NOTE IT IS ASSUMED THAT SAMPLES AND FILES ARE IN THE SAME ORDER!!!!
    #These buffer variables are to prevent too many writes to disk. This allows the
    #user to decide how many reads they'd like to hold in memory before making a write to the disk.
    try:
        buffer_limit = int(buffer_limit)
    except:
        sys.exit("buffer_limit must be an integer or something that can be case as an integer")
    buffer_string = ""
    results = open(output,'w')
    cdef int buffer_count = 0
    cdef int total_tests = 0
    cdef int index
    cdef char* position
    cdef char* chrom
    cdef char* strand
    cdef char* mc_class
    cdef bytes methylated
    cdef int sample_index
    
    lines = {}
    fields = {}
    file_handles = {}
    for file,sample in zip(files,samples):
        file_handles[sample] = open(file,"r")
    #for sample in samples:
    for sample_index in range(len(samples)):
        #lines[sample]=file_handles[sample].readline().rstrip()
        lines[samples[sample_index]]=file_handles[samples[sample_index]].readline().rstrip()
        #fields[sample]=lines[sample].split("\t")
        fields[samples[sample_index]]=lines[samples[sample_index]].split("\t")

    #while there are still lines to be read
    while [line for line in lines.values() if len(line) > 0]:
        min_pos = str(min([int(item[1]) for sample,item in fields.items() if len(lines[sample])>0]))
        counts = []
        current_samples = []
        mc = []
        h = []
        frac = []
		#for sample in samples:
        for sample_index in range(len(samples)):    
            #if len(lines[sample]) > 0:
            if len(lines[samples[sample_index]]) > 0:
                #value = fields[sample]
                value = fields[samples[sample_index]]
                position = value[1]
                if position == min_pos and int(value[5]) >= min_cov:
                    chrom = value[0]
                    strand = value[2]
                    mc_class = value[3]
                    mc.append(value[4])
                    h.append(value[5])
                    methylated = value[6]
                    if methylated == "1":
                        counts.append([int(value[4]),int(value[5])-int(value[4])])
                        frac.append(float(value[4])/float(value[5]))
                    else:
                        counts.append([0,int(value[5])])
                        frac.append(0.0)
                    current_samples.append(samples[sample_index])
                else:
                    mc.append("NA")
                    h.append("NA")
                    frac.append("NA")
            else:
                mc.append("NA")
                h.append("NA")
                frac.append("NA")
                
        #counts = np.array(counts)
        #make sure there's enough data to test
        if len(current_samples)> 1:
            #counts = counts.reshape(len(counts),2)
            #Make sure there's actually a potential difference to be tested
            measured_fractions = [i for i in frac if i != "NA"]
            if frac.count(measured_fractions[0]) != len(measured_fractions):
                stat, pvalue,residuals = chisquare_test(counts,True)
                if np.isnan(stat) == False:
                    #get order of methylation levels for samples
                    #sample_order = "NA"
                    sample_order = [i[1] for i in sorted([((float(count[0]) / float(count[0] + count[1])),sample) for sample,count in zip(current_samples,counts)])]
                    mc_resid = []
                    uc_resid = []
                    index = 0
                    #for sample in samples:
                    for sample_index in range(len(samples)):                  
                        #if sample in current_samples:
                        if samples[sample_index] in current_samples:                      
                            mc_resid.append(residuals[index][0])
                            uc_resid.append(residuals[index][1])
                            index += 1
                        else:
                            mc_resid.append("NA")
                            uc_resid.append("NA")
                    mc_resid = "\t".join(map(str,mc_resid))
                    uc_resid = "\t".join(map(str,uc_resid))
                    mc = "\t".join(map(str,mc))
                    h ="\t".join(map(str,h))
                    frac = "\t".join(map(str,frac))
                    sample_order = ",".join(sample_order)
                    #pdb.set_trace()
                    if min_pos == "84424376":
                        print map(str,[chrom,position,strand,mc_class,sample_order,pvalue,mc,h,frac,mc_resid,uc_resid])
                    buffer_string+="\t".join(map(str,[chrom,position,strand,mc_class,sample_order,pvalue,mc,h,frac,mc_resid,uc_resid]))+"\n"
                    buffer_count+=1
                    if buffer_count == buffer_limit:
                        results.write(buffer_string)
                        buffer_count = 0
                        buffer_string =""
                    total_tests += 1
        for sample,field in fields.items():
            if len(lines[sample]) > 0 and field[1] == min_pos:
                lines[sample]=file_handles[sample].readline().strip()
        for sample,line in lines.items():
            fields[sample]=lines[sample].split("\t")
    results.write(buffer_string)
    results.close()
    check_call("sort -k 6g,6g -o "+output+" "+output,shell=True)
    return total_tests
"""
def chisquare_test_general(obs,use_n_minus_one=False):
    """
    obs is a numpy array shaped the correct way (i.e., a 3x3 matrix for a 3x3 test)
    See this link for a discussion of the N-1 chisquare test and other helpful chisquare tips:
    https://sites.google.com/a/lakeheadu.ca/bweaver/Home/statistics/notes/chisqr_assumptions
    They are largely based on this paper:
    http://www.ncbi.nlm.nih.gov/pubmed/17315184
    """
    warnings.simplefilter("ignore",RuntimeWarning)
    if not(isinstance(obs,np.ndarray)):
        sys.exit("obs argument to chisquare_test must be a numpy array")
    exp = np.empty(obs.shape,dtype=float)
    #This is the denominator for computing the adjusted standardized residuals
    resid_denom = np.empty(obs.shape,dtype=float) 
    total = np.sum(obs,dtype=float)  
    for row_ind in xrange(0,obs.shape[0]):
        rowsum = np.sum(obs[row_ind,:],dtype=float)
        resid_std_denom = (1-(rowsum/total))
        for col_ind in xrange(0,obs.shape[1]):
            #exp[row_ind][col_ind] = (np.sum(obs[row_ind,:],dtype=float) / np.sum(obs,dtype=float)) * \
                                    #(np.sum(obs[:,col_ind],dtype=float) / np.sum(obs,dtype=float)) *np.sum(obs,dtype=float)
            colsum = np.sum(obs[:,col_ind],dtype=float)            
            exp[row_ind][col_ind] = (rowsum * colsum) / total
            resid_denom[row_ind][col_ind] = (exp[row_ind][col_ind] * resid_std_denom * (1-(colsum/total))) ** 0.5 
    #2x2 table
    if obs.shape[0] == 2 and obs.shape[1] == 2: 
        if use_n_minus_one == False:                   
            if len(exp[exp==0]) != 0:
                stat = np.nan
                pvalue = np.nan
                residuals = np.nan
            elif len(exp[exp<5]) / float(exp.size) > .2:
                stat = np.nan
                pvalue = np.nan
                residuals = np.nan
            else:
                stat = np.nansum(((obs - exp) * (obs - exp)) / exp)
                residuals = (obs-exp) / resid_denom
                df = reduce(operator.mul,np.array(obs.shape)-1)
                pvalue = 1 - sci.chi2.cdf(stat,df)
        else:
            #This is one of the assumptions of the N-1 chisquare test
            if len(exp[exp<1]) == 0:
                stat = np.nansum(((obs - exp) * (obs - exp)) / exp)
                stat = stat * (np.nansum(obs)-1)/np.nansum(obs)
                residuals = (obs-exp) / resid_denom
                df = reduce(operator.mul,np.array(obs.shape)-1)
                pvalue = 1 - sci.chi2.cdf(stat,df)
            else:
                #print("Warning: N-1 assumption failed! Using Fisher's exact test")
                stat, pvalue = sci.fisher_exact(obs)
                residuals = (obs-exp) / resid_denom
    #bigger than 2x2 table
    else:
        if len(exp[exp==0]) != 0:
            stat = np.nan
            pvalue = np.nan
            residuals = np.nan
        elif len(exp[exp<5]) / float(exp.size) > .2:
            stat = np.nan
            pvalue = np.nan
            residuals = np.nan
        else:
            stat = np.nansum(((obs - exp) * (obs - exp)) / exp)
            residuals = (obs-exp) / resid_denom
            df = reduce(operator.mul,np.array(obs.shape)-1)
            if np.isnan(stat):
                pvalue = np.nan
            else:
                pvalue = 1 - sci.chi2.cdf(stat,df)
    return stat,pvalue,residuals
"""
def chisquare_test(obs, use_n_minus_one=False):
    NOTE THIS FUNCTION ONLY WORKS FOR 2D chisquare tests
    obs is a 2D array shaped the correct way (i.e., a 3x3 matrix for a 3x3 test)
    See this link for a discussion of the N-1 chisquare test and other helpful chisquare tips:
    https://sites.google.com/a/lakeheadu.ca/bweaver/Home/statistics/notes/chisqr_assumptions
    They are largely based on this paper:
    http://www.ncbi.nlm.nih.gov/pubmed/17315184
    #warnings.simplefilter("ignore",RuntimeWarning)
    #if not(isinstance(obs,np.ndarray)):
    #    sys.exit("obs argument to chisquare_test must be a numpy array")
    cdef int nrow
    cdef int ncol
    cdef int row_ind
    cdef int col_ind
    cdef long rowsum
    cdef long colsum
    #Number of cells below 5 observations 
    cdef int low_obs_count = 0
    #The number of cells (rows*cols)
    cdef int size
    #The number of cells that have to be below 5 observations
    #for the test to be skipped
    cdef double low_obs_cutoff = 0
    
    stat = 0
    nrow = len(obs)
    ncol = len(obs[0])
    size = nrow*ncol
    low_obs_cutoff = 0.2 * size
    exp = [[0]*ncol for i in xrange(nrow)]
    #This is the denominator for computing the adjusted standardized residuals
    resid_denom = [[0]*ncol for i in xrange(nrow)]
    residuals =  [[0]*ncol for i in xrange(nrow)]
    total = float(sum([sum(i) for i in obs]))

    #2x2 table
    if nrow == 2 and ncol == 2: 
        if use_n_minus_one == False:
            for row_ind in xrange(0,nrow):
                rowsum = sum(obs[row_ind])
                row_std_resid = (1-(rowsum/total))
                for col_ind in xrange(0,ncol):
                    colsum = sum([i[col_ind] for i in obs])
                    exp[row_ind][col_ind] = (rowsum * colsum) / total
                    if exp[row_ind][col_ind] < 5:
                        low_obs_count+=1
                    if exp[row_ind][col_ind] == 0 or low_obs_count > low_obs_cutoff:
                        stat = np.nan
                        pvalue = np.nan
                        residuals = np.nan
                        return stat,pvalue,residuals
                    resid_denom[row_ind][col_ind] = (exp[row_ind][col_ind] * row_std_resid * (1-(colsum/total))) ** 0.5
                    stat += ((obs[row_ind][col_ind] - exp[row_ind][col_ind]) ** 2) / exp[row_ind][col_ind]
                    residuals[row_ind][col_ind] = (obs[row_ind][col_ind] - exp[row_ind][col_ind]) / resid_denom[row_ind][col_ind]

                                        
            df = (nrow - 1 ) * (ncol - 1)
            pvalue = 1 - sci.chi2.cdf(stat,df)
            return stat,pvalue,residuals
                    
        else:
            #This is one of the assumptions of the N-1 chisquare test
            for row_ind in xrange(0,nrow):
                rowsum = sum(obs[row_ind])
                row_std_resid = (1-(rowsum/total))
                for col_ind in xrange(0,ncol):
                    colsum = sum([i[col_ind] for i in obs])
                    exp[row_ind][col_ind] = (rowsum * colsum) / total
                    if exp[row_ind][col_ind] < 1:
                        low_obs_count += 1
                    if exp[row_ind][col_ind] == 0:
                        #To prevent division by 0
                        exp[row_ind][col_ind]= 0.0001
                    resid_denom[row_ind][col_ind] = (exp[row_ind][col_ind] * row_std_resid * (1-(colsum/total))) ** 0.5
                    stat += ((obs[row_ind][col_ind] - exp[row_ind][col_ind]) ** 2) / exp[row_ind][col_ind]
                    residuals[row_ind][col_ind] = (obs[row_ind][col_ind] - exp[row_ind][col_ind]) / resid_denom[row_ind][col_ind]
            if low_obs_count > 0:
                stat, pvalue = sci.fisher_exact(obs)
            else:
                stat = stat * (total-1)/total
                df = (nrow - 1 ) * (ncol - 1)
                pvalue = 1 - sci.chi2.cdf(stat,df)
            return stat,pvalue,residuals

    #bigger than 2x2 table
    else:
        for row_ind in xrange(0,nrow):
            rowsum = sum(obs[row_ind])
            row_std_resid = (1-(rowsum/total))
            for col_ind in xrange(0,ncol):
                colsum = sum([i[col_ind] for i in obs])
                exp[row_ind][col_ind] = (rowsum * colsum) / total
                if exp[row_ind][col_ind] < 5:
                    low_obs_count+=1
                if exp[row_ind][col_ind] == 0 or low_obs_count > low_obs_cutoff:
                    stat = np.nan
                    pvalue = np.nan
                    residuals = np.nan
                    return stat,pvalue,residuals
                resid_denom[row_ind][col_ind] = (exp[row_ind][col_ind] * row_std_resid * (1-(colsum/total))) ** 0.5
                stat += ((obs[row_ind][col_ind] - exp[row_ind][col_ind]) ** 2) / exp[row_ind][col_ind]
                residuals[row_ind][col_ind] = (obs[row_ind][col_ind] - exp[row_ind][col_ind]) / resid_denom[row_ind][col_ind]

        df = (nrow - 1 ) * (ncol - 1)
        if np.isnan(stat):
            pvalue = np.nan
        else:
            pvalue = 1 - sci.chi2.cdf(stat,df)
        return stat,pvalue,residuals
"""
def root_mean_square_test(obs,nsim=10000):
    stat = 0
    nrow = len(obs)
    ncol = len(obs[0])
    size = nrow*ncol
    exp = [[0]*ncol for i in xrange(nrow)]
    #This is the denominator for computing the adjusted standardized residuals
    resid_denom = [[0]*ncol for i in xrange(nrow)]
    residuals =  [[0]*ncol for i in xrange(nrow)]
    total = float(sum([sum(i) for i in obs]))
    exp_prob = [[0]*ncol for i in xrange(nrow)]
    for row_ind in xrange(0,nrow):
        rowsum = sum(obs[row_ind])
        row_std_resid = (1-(rowsum/total))
        for col_ind in xrange(0,ncol):
            colsum = sum([i[col_ind] for i in obs])
            exp[row_ind][col_ind] = (rowsum * colsum) / total
            exp_prob[row_ind][col_ind] = exp[row_ind][col_ind] / total
        resid_denom[row_ind][col_ind] = (exp[row_ind][col_ind] * row_std_resid * (1-(colsum/total))) ** 0.5
        stat += ((obs[row_ind][col_ind] - exp[row_ind][col_ind]) ** 2)
        residuals[row_ind][col_ind] = (obs[row_ind][col_ind] - exp[row_ind][col_ind]) / resid_denom[row_ind][col_ind]
    stat = (stat / total) ** 0.5
    
    #compute p-value
    pvalue = simulate_rms_test_pvalue(exp_prob,total,stat)
    return stat,pvalue,residuals

def simulate_rms_test_pvalue(exp, num_obs, obs_stat,num_sims=10000):
    sig_count = 0
    nrow = len(exp)
    ncol = len(exp[0])
    num_cells = nrow * ncol
    exp_flatten = [item for sublist in exp for item in sublist]
    for sim in xrange(num_sims):
        stat = 0
        sampling = np.random.multinomial(num_obs,exp_flatten)
        obs = [[sampling[i],sampling[i+1]] for i in xrange(0,num_cells,2)]
        for row_ind in xrange(0,nrow):
            rowsum = sum(obs[row_ind])
            row_std_resid = (1-(rowsum/num_obs))
            for col_ind in xrange(0,ncol):
                colsum = sum([i[col_ind] for i in obs])
                exp[row_ind][col_ind] = (rowsum * colsum) / num_obs
            stat += ((obs[row_ind][col_ind] - exp[row_ind][col_ind]) ** 2)
        stat = (stat / num_obs) ** 0.5
        if stat >= obs_stat:
            sig_count += 1
    return sig_count / float(num_sims)
    
def float_range(a, b, inc):
    """
    Associated with QVALUE function.
    Function to simulated the python range function, except for floats.
    """
    floats = []
    while a < b + inc/2: #add inc/2 to make sure that the last float is included even if it off by a few decimals bc
        floats.append(a) #of errors when handling floats
        a += inc
    return floats
  
def qvalue(p, lambd=float_range(0,0.90,0.05), pi0_method="bootstrap", robust=False, 
           smooth_df = 3, smooth_log_pi0 = False):
    """
    Translated from R function; 'smoothed' option NOT IMPLEMENTED
    Input
    =============================================================================
    p: a numpy array of p-values (only necessary input)
    lambda: the value of the tuning parameter to estimate pi0 (optional)
    pi0.method: either "smoother" or "bootstrap"; the method for automatically
               choosing tuning parameter in the estimation of pi0, the proportion
               of true null hypotheses
    robust: an indicator of whether it is desired to make the estimate more robust
            for small p-values and a direct finite sample estimate of pFDR (optional)
    smooth.df: degrees of freedom to use in smoother (optional)
    smooth.log.pi0: should smoothing be done on log scale? (optional)
    
    Output
    =============================================================================
    call: gives the function call
    pi0: an estimate of the proportion of null p-values
    qvalues: a vector of the estimated q-values (the main quantity of interest)
    pvalues: a vector of the original p-values
    """
    #preprocessing
    if min(p)<0 or max(p)>1:
        sys.exit("ERROR: p-values not in valid range.")
    p = np.array(p)
    
    lambda_length = len(lambd)
    if lambda_length>1:
        if lambda_length<4:
            sys.exit("ERROR: If length of lambda greater than 1, you need at least 4 values.")
        if min(lambd)<0 or max(lambd)>=1:
            sys.exit("ERROR: Lambda must be within [0,1).")
            
    m = len(p)
    #Ways to estimate pi0
    if lambda_length==1:
        if lambd<0 or lambd>=1:
            sys.exit("ERROR: MAbda must be within [0,1).")
        pi0 = np.mean(p>=lambd)/(1-lambd)
        pi0 = min(pi0, 1)
    else:
        pi0 = np.array([0.0]*lambda_length)
        for i in xrange(lambda_length):
            pi0[i] = np.mean(p>=lambd[i])/(1-lambd[i])

        if pi0_method=="smoother":
            sys.exit("Not implemented")
            if smooth_log_pi0:
                pi0 = np.log(pi0)
                
            #spi0 = smooth.spline(lambd, pi0, df=smooth_df)
            #pi0 = predict(spi0, x=max(lambd))$y
            
            if smooth_log_pi0:
                pi0 = np.exp(pi0)
            pi0 = min(min(pi0),1)
        elif pi0_method=="bootstrap":
            minpi0 = min(pi0)
            mse = np.array([0.0]*lambda_length)
            pi0_boot = np.array([0.0]*lambda_length)
            for i in xrange(100):
                p_boot = []
                for j in xrange(m):
                    p_boot.append(rand.choice(p))
                p_boot=np.array(p_boot)
                for i in xrange(lambda_length):
                    pi0_boot[i] = np.mean(p_boot > lambd[i])/(1-lambd[i])
                mse = mse + np.square(pi0_boot-minpi0)
            pi0 = min(pi0[mse==min(mse)])
            pi0 = min(pi0,1)
        else:
            sys.exit("ERROR: pi0_method must be either 'smoother' or 'bootstrap'.")
            
    if pi0<=0:
        sys.exit("ERROR: The estimated pi0 <= 0. Check that you have valid p-values or use another lambda method.")

    #The estimated q-values calculated here
    u = np.argsort(p)
    v = qvalue_rank(p)

    qvalue = pi0*m*p/v
    if robust:
        qvalue = pi0*m*p/(v*(1-np.power(1-p,[m]*len(p))))
    qvalue[u[m-1]] = min(qvalue[u[m-1]],1)
    for i in xrange((m-2),-1,-1):
        qvalue[u[i]] = min(qvalue[u[i]], qvalue[u[i+1]],1)
    
    #The results are returned
    return (qvalue, p)

def qvalue_rank(x):
    """ 
    Ranking function which returns number of observations less than or equal
    This function calculates the rank of each value in the list by assigning it a number based on the number of values smaller or equal to it.
    
    ex: [1,2,3] -> [1,2,3] 
    [1,1,1] -> [3,3,3]
    [1,1,2] -> [2,2,3]
    """
    orig_index_order = np.argsort(x) #keep track of original index ordering of the numbers before sorting
    sorted_list = np.sort(x) #sorted list
    
    cumsum = [0] #put temporary sums here
    prev_num = "NA"
    count = 0
    for num in sorted_list:
        if num!=prev_num and prev_num!="NA": #if the number does not match previous, then append its sum to the list
            for i in xrange(count - cumsum[-1]): #append the sum the number of times that it repeated
                cumsum.append(count)
        prev_num = num
        count+=1
    for i in xrange(count - cumsum[-1]):
        cumsum.append(count)
    cumsum.remove(0)
    
    result = [0]*len(sorted_list)
    
    index = 0
    for num in orig_index_order: #put all of the sum values back in original ordering
        result[num] = cumsum[index] 
        index += 1
    
    return np.array(result)

if __name__ == '__main__':
    pass