'''
Created on Dec 16, 2011

@author: Matt
'''
import sys
import subprocess
import os
import itertools
import time
import errno
import shlex

try:   
    import cPickle
except:
    sys.exit("methylpy requires the cPickle module")
try:
    from methylpy.gelist import gelist, sort_by_coord_gel, sort_by_coord_tuple, write_gelist_Line
except:
    sys.exit("methylpy.mCdata requires the methylpy.gelist module")
try:
    from methylpy.utilities import print_checkpoint, connect_dbs, parse_gff_line
except:
    sys.exit("methylpy.mCdata requires the methylpy.utilities module")
try:
    import MySQLdb
except:
    sys.exit("methylpy.mCdata requires the MySQLdb module")
try:
    from multiprocessing import Pool, Pipe, Process
except:
    sys.exit("methylpy.mCdata requires Pool")
try:
    import numpy as np
except:
    sys.exit("methylpy.mCdata requires the numpy module")
import pdb
def create_smp_matrix(samples, servers, databases, coverage, chromosomes=False, to_file=False):
    """
    create_smp_matrix generates a matrix of chromosome, strand, position, class for all of the samples in the sample list
    where that coordinate exists in the mc table and also has an h count greater than 4 in the allc tables.
    
    samples is a list of samples
    servers is a list of servers
    databases is a list of databases
    chromosomes is a list of chromosomes to use (False to use all)
    to_file writes the matrix to a file 
    coverage is the minimum h count that you want to use
    """
    #input check
    if not isinstance(samples, list):
        if isinstance(samples, basestring):
            samples = [samples]
        else:
            sys.exit("samples must be a list of strings")
    if not isinstance(servers, list):
        if isinstance(servers, basestring):
            servers = [servers]
        else:
            sys.exit("servers must be a list of strings")
    if not isinstance(databases, list):
        if isinstance(databases, basestring):
            databases = [databases]
        else:
            sys.exit("databases must be a list of strings")
    
    #First go to the mc tables and generate a dictionary of all positions that exist in each sample table;
    if samples!="all" and len(samples)!=len(servers)!=len(databases):
        sys.exit("Length of sample, server, and database list should be the same")
        
    db_connections, sample_dbs, cursors, samples = connect_dbs(samples, servers, databases, "mc")
    allc_db_connections, allc_sample_dbs, allc_cursors, samples = connect_dbs(samples, servers, databases, "allc")

    if chromosomes==False:
        (server, database) = sample_dbs[samples[0]]
        cursors[(server,database)].execute("select table_name from information_schema.tables where table_schema='" + \
                                           database + "' and table_name regexp 'mc_" + samples[0] + "_.{1,2}$';")
        results = cursors[(server,database)].fetchall()
        chromosomes = []
        for table in results:
            table = table[0]
            if table[-1]=="C" or table[-1]=="M":
                continue
            chromosomes.append(table[table.rfind("_")+1:])
    
    if to_file!=False:
        f = open(to_file, "w")
        f.write("Chromosome, Strand, Position, Class \t")
        for sample in samples:
            f.write(sample + "\t")
        f.write("\n")

    samples.sort()
    samplelist = {}
    samplelist["samples"] = samples
    print_checkpoint("Getting sample info")
    for chromosome in chromosomes:
        for sample in allc_sample_dbs:
            print_checkpoint("Getting allc info for: " + sample + " chromosome " + chromosome)
            (server, database) = allc_sample_dbs[sample]
            if allc_cursors[(server,database)].execute("select position, strand, class,methylated from allc_" + \
                                                       sample + "_" + chromosome + " where h>" + coverage + ";")==0:
                print "Warning: Sample(" + sample + ") chromosome " + chromosome + " allc table does not exist"
                continue
            rows = allc_cursors[(server,database)].fetchall()
            for row in rows:
                position = row[0]
                strand = row[1]
                clss = row[2]
                methylated = row[3]
                if clss[1]=="G":
                    clss = "CG"
                elif clss[2]=="G":
                    clss = "CHG"
                else:
                    clss = "CHH"
                if methylated == 1:
                    try: 
                        samplelist[(chromosome,strand,position,clss)][sample] = 1
                    except KeyError:
                        samplelist[(chromosome,strand,position,clss)] = {}
                        samplelist[(chromosome,strand,position,clss)][sample] = 1
                else:
                    try: 
                        samplelist[(chromosome,strand,position,clss)][sample] = 0
                    except KeyError:
                        samplelist[(chromosome,strand,position,clss)] = {}
                        samplelist[(chromosome,strand,position,clss)][sample] = 0
            allc_cursors[(server,database)].execute("RESET QUERY CACHE;")

        if "samples" in samplelist:
            del samplelist[("samples")]
        if to_file!=False:
            print_checkpoint("Sorting and writing to file chromosome " + chromosome)
            for chrom,strand,position,clss in sorted(samplelist.iterkeys()):
                f.write(chrom+"\t" + strand + "\t" + str(position) + "\t" + clss + "\t")
                for sample in samples:
                    if sample in samplelist[(chrom,strand,position,clss)]:
                        f.write(str(samplelist[(chrom,strand,position,clss)][sample]) + "\t")
                    else:
                        f.write("NA \t")
                f.write("\n")
        samplelist={}
        
    f.close()
    for server,database in db_connections:
        db_connections[(server,database)].close()
        cursors[(server,database)].close()
    for server,database in allc_db_connections:
        allc_db_connections[(server,database)].close()
        allc_cursors[(server,database)].close()
        
    print_checkpoint("Finished")
    return

def assign_mC_values(gel,servers,databases,samples,table_type, processors=1, save_result=False, by_element=False,
                     use_c_index=False, save_memory=False,store_level=False,store_method="weighted",store_reads=True):
    """
    This function extracts information from mc tables and saves it in a gff format
    The fields that can be accessed in the attributes column are:
    mc_positions: a list of genomic coordinates for the mc sites
    mc_classes: a list of the mc types of each position (CG, CHG, or CHH)
    mc_counts: a list of the number of methylated reads covering each position
    h_counts: a list of the number of reads covering each position
    strands: a list of the strands for each position
    NOTE: Each of these "lists" is actually a comma separated string, so they need to first be
    split to make proper python lists
    
    Arguments:
    gel is a gelist object
    server is a list of names of servers the tables are on
    database is a list of the databases the tables are in
    samples is a list of the sample names (NOTE MUST BE A LIST EVEN IF THERE'S ONLY ONE SAMPLE!)
    table_type indicates whether or not to extract the information from the allc tables or mc tables
    processors indicates how many cpus you'd like to parallelize this function to
    save_result indicates whether or not you'd like gff files saved to disk (REQUIRED IF # PROCESSORS>1)
    by_element indicates that you'd like the returned value to be a dictionary of the elements of the gel
    rather than a dictionary of the samples
    use_c_index indicates that you'd like to load a c_index to determine the C class (CG,CHG, or CHH) of 
    unmethylated cytosines.  This argument must be the path to a base index where the sample specific indexes
    are coded up to the chr at the very end. Example: base_index=c_index_chrN sample_index = c_index_sample_chrN
    /path/to/c_index.
    save_memory indicates whether or not you'd like to keep all the results in memory until the creation of the 
        final bedtool. In general you should set this to false, but if the data set you're working with is large, 
        set to true and each sample's results will be written to a temporary file, which will all be cat'ed 
        together at the end before bedtool creation.
    store_levels is a list of C contexts for which you want to precompute methylation levels. These will be stored
        in the attributes column as <ctype>_<method>_<use all/mc_reads>_level (e.g., CAA_mattia_True_level is CAA_m_T=xxx
        or CHG_weighted_False_level is CHG_w_F=xxx)
    store_method indicates which method (weighted, mattia, etc.) you'd like to use to precompute the mC levels
    store_reads is a boolean indicating whether or not you'd like to include all reads that cover an element in the 
    mC level calculation or just the ones covering mCs
    """
    gff_string = ""

    if not isinstance(samples, list):
        if isinstance(samples, basestring):
            samples = [samples]
        else:
            sys.exit("samples must be a list")
    if not isinstance(servers, list):
        if isinstance(servers, basestring):
            servers = [servers]
        else:
            sys.exit("servers must be a list")
    if not isinstance(databases, list):
        if isinstance(databases, basestring):
            databases = [databases]
        else:
            sys.exit("databases must be a list")

    if not(isinstance(gel,gelist)):
        sys.exit("assign_mC_values expects a gelist as the argument to gel.")
    if processors>1 and save_result==False:
        sys.exit("assign_mC_values requires a save_result parameter when using multiprocessing")
    if samples == "all":
        samples = []
        conn_list = zip(servers,databases)
        servers = []
        databases = []
        for server,database in conn_list:
            try:
                connection = MySQLdb.connect(server,"mysql","rekce",database)
                cursor = connection.cursor()
            except:
                sys.exit("Could not establish database connection with these parameters: server="+server+ \
                         " database="+database)
            cursor.execute("show tables like '"+table_type+"_%_1'")
            rows = cursor.fetchall()
            for row in rows:
                sample = row[0][len(table_type)+1:row[0].rfind("_")]
                samples.append(sample)
                databases.append(database)
                servers.append(server)
            connection.close()   

    if use_c_index != False:
        set_c_index(use_c_index)
        load_base_c_index(gel.chromosomes)
        
    results = {}
    if(processors>1):
        for sample,server,database in zip(samples,servers,databases):
            assign_mc_pipes(gel, sample, server, database, table_type, processors=processors, use_c_index=use_c_index,
                            save_memory=save_memory,store_level=store_level,store_method=store_method,
                            store_reads=store_reads)        
    else:
        for sample,server,database in zip(samples,servers,databases):
            gff_string += create_samplestring(gel, sample, server, database, table_type,use_c_index=use_c_index,
                                              save_memory=save_memory,store_level=store_level,store_method=store_method,
                                              store_reads=store_reads)
                
    if save_memory==False and processors==1:        
        gc = gelist(gff_string,from_string=True)
    else:
        gc = gelist(combine_mc_files(samples, save_result))
        
    if save_result != False and save_memory==False and processors==1:
        gc.saveas(save_result)

    return gc

def combine_mc_files(samples, save_result, path="tmp_results/"):
    """
    Concatenates the files containing lines that have mc values gotten from assign_mc_values. 
    This will happen if save_memory option is enabled or if multiprocessing was used. 
    You can specify your own folder if you want using the path parameter.
    Samples is the list of samples used.
    Save_result is what you want your file named.
    """
    print_checkpoint("Concatenating files..")

    save = shlex.split("find " + path +" -type f -print | xargs cat > allsamplestmp.gff")
    subprocess.check_call(save) #cat file file file > allsamplestmp.gff
    cmd = shlex.split("rm -r " + path)
    subprocess.check_call(cmd) #rm file file file
        
    temp_pathname = os.path.abspath(os.path.join(os.getcwd(), "allsamplestmp.gff"))
    new_pathname = os.path.abspath(os.path.join(os.getcwd(), save_result))
    os.rename(temp_pathname, new_pathname)
    return new_pathname
    
def assign_mc_pipes(gel, sample, server, database, table_type, processors=1, use_c_index=False,chromosomes=False,
                    save_memory=False,store_level=False,store_method="weighted",store_reads=True):
    """
    Modified version of create_sample_string that is used when multiprocessing using pipes is needed. This function should
    be called from assign_mc_values so users should not have to call it. It sends all required info to its children so they
    can call get_mysql_data and then write the results to file. Those temp files are combined at the end of the function.
    
    see assign_mC_values for specific parameter info
    """
    #start children
    all_conns = []
    sql_conns = []
    process = []
    for proc in xrange(processors): #create all the child/parent pipes for processors
        parent_conn, child_conn = Pipe()
        connection = MySQLdb.connect(server,"mysql","rekce",database)
        cursor = connection.cursor()
        process.append(Process(target=child_sql, args=(child_conn, parent_conn, cursor)))
        process[-1].start()
        child_conn.close()
        all_conns.append(parent_conn)
        sql_conns.append((connection,cursor))

    c_index = False
    current_index = 0
    gff_string = ""
    print_checkpoint("\tBegin sample "+sample)
    #used to keep track of which chromosome we're on so we know when to load the next c_index
    current_chrom = ""

    #need to sort the gelist so that you can smartly load the c_indexes (i.e., one chromosome at a time)
    for element in sorted(gel, cmp=sort_by_coord_gel):
        if element["_format"] == "gff":
            chrom = element[0]
            start = element[3]
            end = element[4]
            strand = element[6]
            element_type = element[2]
            attributes = element[8]
        else:
            chrom = element[0]
            #Tables are 1 based
            start = str(int(element[1])) #removed +1, used to be str(int(element[1]) +1 )
            end = str(int(element[2]))   #removed +1 here as well
            if element["_fields"] > 3:
                strand = element[3]
                if element["_fields"] > 4:
                    element_type = element[4]
                else:
                    element_type = "none"
            else:
                strand = "."
                element_type = "none"
            attributes = ""
        
        if use_c_index!=False:
            
            if chrom != current_chrom:
                current_chrom = chrom
                c_index = load_sample_c_index(chrom, sample)
                #Test to make sure c_index loaded
                if c_index == False:
                    return gff_string
        
        #send all required information to use get_mysql_data in child
        all_conns[current_index].send((chrom, start, end, strand, element_type, sample, table_type, 
                                       attributes, c_index, store_level, store_method, store_reads))
        current_index += 1
        if current_index >= processors:
            current_index = 0
    
    for parent_conn in all_conns:
        parent_conn.close()
    for p in process:
        p.join()
    for connection, cursor in sql_conns:
        cursor.close()
        connection.close()
    
    return 

def child_sql(child_conn, parent_conn, cursor):
    """
    Child processes that are piped line data from the gelist. They run get_mysql_data and write the results to a file.
    
    child_conn is the child end of the pipe
    parent_conn is the parent end, should not be used in this function
    cursor is the cursor to the database you will be querying
    """
    parent_conn.close() #close unused end of the pipe
    ensure_dir("tmp_results") #make temporary directory to put results in
    while True:
        try:
            chrom, start, end, strand, element_type, sample, table_type, \
            attributes, c_index, store_level, store_method, store_reads = child_conn.recv() 
        except EOFError:
            child_conn.close() #close pipe when you reach end of file
            break

        f = open("tmp_results/" + str(chrom) + "-" + str(start) + "-" + str(sample), 'w')    
        temp_string = get_mysql_data(chrom,start,end,strand,element_type,sample,table_type,cursor,attributes,c_index, 
                                 store_level=store_level, store_method=store_method, store_reads=store_reads)
        f.write(temp_string) #write results to temporary file
        f.close() #close everything  

def ensure_dir(d):
    """
    Create directory d; if this fails check if it already exists
    """
    try:
        os.makedirs(d)
    except OSError, e:
        if e.errno == errno.EEXIST:
            #it already exists, so continue on
            pass
        else:
            #some other error occurred, so raise it
            raise

def create_samplestring(gel, sample, server,database,table_type, use_c_index=False,chromosomes=False,save_memory=False,
                        store_level=False,store_method="weighted",store_reads=True):
    """
    This function is used with assign mc values in order to create a string of information for the given sample. The info
    is fetched from an SQL table using the function get_mysql_data. Depending on the file format, fields have to be taken from
    different columns of the genome file. If store_level is true, then the methylation levels will be put in the gff file 
    as well.
    
    Called from assign_mC_values or assign_mC_pipes
    see assign_mC_values for specific parameter info
    """
    c_index = False
    if save_memory == True:
        path = "tmp_results"
        ensure_dir(path) #make temporary directory to put results in
        f=open(path + "/" + str(sample)+"tmp.gff", "w")
    gff_string = ""
    print_checkpoint("\tBegin sample "+sample)
    #used to keep track of which chromosome we're on so we know when to load the next c_index
    current_chrom = ""
    connection = MySQLdb.connect(server,"mysql","rekce",database)
    cursor = connection.cursor()
    #need to sort the pybedtool so that you can smartly load the c_indexes (i.e., one chromosome at a time)
    for element in sorted(gel,cmp=sort_by_coord_gel):
        if element["_format"] == "gff": #gff format
            chrom = element[0]
            start = element[3]
            end = element[4]
            strand = element[6]
            element_type = element[2]
            attributes = element[8]
        else: #otherwise, it is assumed to be a bed format
            chrom = element[0]
            #Tables are 1 based
            start = str(int(element[1])) #removed +1, used to be str(int(element[1]) +1 )
            end = str(int(element[2]))   #removed +1 here as well
            if element["_fields"] > 3:
                strand = element[3]
                if element["_fields"] > 4:
                    element_type = element[4]
                else:
                    element_type = "none"
            else:
                strand = "."
                element_type = "none"
            attributes = ""
        
        if use_c_index!=False:
            
            if chrom != current_chrom:
                current_chrom = chrom
                c_index = load_sample_c_index(chrom, sample)
                #Test to make sure c_index loaded
                if c_index == False:
                    return gff_string
                
        temp_string = get_mysql_data(chrom,start,end,strand,element_type,sample,table_type,cursor,attributes, c_index, 
                                     store_level=store_level, store_method=store_method, store_reads=store_reads)
        if save_memory == False: #keep all the information in memory
            gff_string += temp_string
        else:
            f.write(temp_string) #write the string to file
            
    cursor.close()
    connection.close()
    if save_memory == False:
        return gff_string
    else: #if it was written to a file, return an empty string
        f.close()
        return ""

def get_mysql_data(chr,start,end,strand,dmr_type,sample,table_type,cursor,attributes,c_index=False,
                   store_level=False,store_method="weighted",store_reads=True):
    """
    This function pulls the methylation data for a specific range in a specific sample from MySQL.
    chr, start, end, and strand indicate which elements from the table to pull out
        NOTE strand is currently not used!
    dmr_type is the attribute from the feature field of a gff file
    sample is the attribute to be put in the source fields of a gff file
    table_type is the kind of table (mc and stacks or mc and allc) that the methylation data should be drawn from 
    cursor is a cursor object pointing to the appropriate MySQL table
    attributes are the attributes from the gff file (last, flexible column)
    c_index is the c_index object that is used to get the type of C for a particular position (e.g., CG). This
        is more or less deprecated so make allc tables
    store_levels is a list of C contexts for which you want to precompute methylation levels. These will be stored
        in the attributes column as <ctype>_<method>_<all/mc_reads>_level (e.g., CAA_mattia_True_level, or 
        CHG_weighted_False_level)
    store_method indicates which method (weighted, mattia, etc.) you'd like to use to precompute the mC levels
    store_reads is a boolean indicating whether or not you'd like to include all reads that cover an element in the 
        mC level calculation or just the ones covering mCs
    """

    mc_positions = []     #position in the mc table
    mc_classes = []       #class in the mc table
    mc_counts = []        #mc in the mc table
    h_counts = []         #h in the mc table
    mc_strands = []       #strand in the mc table
    
    stacks_positions = [] #position in the stacks table
    stacks_strands = []   #strand is based on whether the C in the stacks table was from the fbase or rbase column
    stacks_call = []      #fcall or rcall in the stacks table
    stacks_total = []     #ftotal or rtotal in the stacks table depending on whether the C was from the fbase/rbase col
    stacks_classes = []   #Based on sequence from the actual
    
    if table_type == "allc":
        query = "select position,strand,class,mc,h,methylated from allc_"+sample+"_"+chr.replace("chr","") \
                + " where position >="+str(start)+" and position <= "+str(end)
        cursor.execute(query)
        rows = cursor.fetchall()
        if len(rows) > 0:
            for row in rows:
                if row[5] == 1:
                    mc_positions.append(str(row[0]))
                    mc_strands.append(row[1])
                    mc_classes.append(row[2])
                    mc_counts.append(str(row[3]))
                    h_counts.append(str(row[4]))
                stacks_positions.append(str(row[0]))
                stacks_strands.append(row[1])
                stacks_classes.append(row[2])
                stacks_total.append(str(row[4]))
    else:       
        query = "select position,strand,class,mc,h from mc_"+sample+"_"+chr.replace("chr","")+" where position >=" + \
                str(start)+" and position <= "+str(end)
        
        cursor.execute(query)
        rows = cursor.fetchall()
        if len(rows) > 0:
            for row in rows:
                mc_positions.append(str(row[0]))
                mc_strands.append(row[1])
                mc_classes.append(row[2])
                mc_counts.append(str(row[3]))
                h_counts.append(str(row[4]))
            """
            Try not using C_INDEX switch -- not implemented yet
            """
            if c_index == False and table_type != "allc":
                query = "select position,fcall,ftotal,fbase from stacks_bs_"+sample+"_"+chr.replace("chr","") + \
                        " where position >="+str(int(start)-2)+" and position <= "+str(int(end)+2)
                cursor.execute(query)
                rows = cursor.fetchall()
                
                if len(rows) > 0:
                    #We're going to get the mC classes from the stacks table so ignore whatever we saw in
                    #the mC tables. This allows us to expand beyond the three contexts
                    mc_classes = []
                    for index in xrange(2,len(rows)-2):
                        row = rows[index]
                        if row[2] > 0 and row[3]=="C":
                            stacks_positions.append(str(row[0]))
                            stacks_call.append(row[1])
                            stacks_strands.append("+")
                            stacks_total.append(str(row[2]))
                            try:
                                c_class = "C"
                                if rows[index+1][0] == rows[index][0] + 1:
                                    c_class += rows[index+1][3]
                                else:
                                    c_class += "N"
                                if rows[index+2][0] == rows[index][0] + 2:
                                    c_class += rows[index+2][3]
                                else:
                                    c_class += "N"                  
                            except:
                                c_class = "N"
                            stacks_classes.append(c_class)
                            if stacks_positions[-1] in mc_positions:
                                #Since we'll go through everything at this stage in positional order
                                #we can simply append the mc_class to the end of the list
                                mc_classes.append(c_class)
        
                query = "select position,rcall,rtotal,rbase from stacks_bs_"+sample+"_"+chr.replace("chr","") + \
                        " where position >="+str(int(start)-2)+" and position <= "+str(int(end)+2)
                cursor.execute(query)
                rows = cursor.fetchall()
                
                if len(rows) > 0:
                    for index in xrange(2,len(rows)-2):
                        row = rows[index]
                        if row[2] > 0 and row[3]=="C":
                            stacks_positions.append(str(row[0]))
                            stacks_call.append(row[1])
                            stacks_strands.append("-")
                            stacks_total.append(str(row[2]))
                            try:
                                c_class = "C"
                                if rows[index-1][0] == rows[index][0] - 1:
                                    c_class += rows[index-1][3]
                                else:
                                    c_class += "N"
                                if rows[index-2][0] == rows[index][0] - 2:
                                    c_class += rows[index-2][3]
                                else:
                                    c_class += "N"                  
                            except:
                                c_class = "N"
                            stacks_classes.append(c_class)
                            if stacks_positions[-1] in mc_positions:
                                #Since we'll go through everything at this stage in positional order
                                #we can simply append the mc_class to the end of the list
                                mc_classes.append(c_class)
    
            else: #should go here if using c_index
                query = "select position,fcall,ftotal from stacks_bs_"+sample+"_"+chr.replace("chr","") + \
                        " where position >="+str(start)+" and position <= "+str(end)+" and fbase = 'C' and ftotal > 0"
                cursor.execute(query)
                rows = cursor.fetchall()
                
                if len(rows) > 0:
                    for row in rows:
                        stacks_positions.append(str(row[0]))
                        stacks_call.append(row[1])
                        stacks_strands.append("+")
                        stacks_total.append(str(row[2]))
                        #Look in the base_c_index for the class
                        try:
                            c_class = c_index[(int(row[0]),"+")]
                        except:
                            #Look in the sample specific c_index for the class
                            try:
                                c_class = base_c_index[chr][(int(row[0]),"+")]
                            except:
                                c_class = "N"
                        stacks_classes.append(c_class)
        
                query = "select position,rcall,rtotal from stacks_bs_"+sample+"_"+chr.replace("chr","") + \
                        " where position >="+str(start)+" and position <= "+str(end)+" and rbase = 'C' and rtotal > 0"
                cursor.execute(query)
                rows = cursor.fetchall()
                
                if len(rows) > 0:
                    for row in rows:
                        stacks_positions.append(str(row[0]))
                        stacks_call.append(row[1])
                        stacks_strands.append("-")
                        stacks_total.append(str(row[2]))
                        #Look in the sample specific c_index for the class
                        try:
                            c_class = c_index[(int(row[0]),"-")]
                        except:
                            #Look in the base_c_index for the class
                            try:
                                c_class = base_c_index[chr][(int(row[0]),"-")]
                            except:
                                c_class = "N"
                        stacks_classes.append(c_class)
    
    features_list = ";".join(["mc_positions="+",".join(mc_positions),
                                "mc_strands="+",".join(mc_strands),
                                "mc_classes="+",".join(mc_classes),
                                "mc_counts="+",".join(mc_counts),
                                "h_counts="+",".join(h_counts),
                                "stacks_positions="+",".join(stacks_positions),
                                "stacks_call="+",".join(stacks_call),
                                "stacks_strands="+",".join(stacks_strands),
                                "stacks_total="+",".join(stacks_total),
                                "stacks_classes="+",".join(stacks_classes)])
    gff_string = "\t".join([chr,sample,dmr_type,str(start),str(end),".",strand,".",attributes+";"+features_list])
    
    if store_level != False:
        if store_reads==True:
            all_reads = "T"
        else:
            all_reads = "F"
        level_field = ','.join(store_level) + "_" + store_method[0] + "_" + all_reads #field name = type,type_method_T/F
        level = get_methylation_level(gff_string, store_level,method=store_method,use_all_reads=store_reads)
        gff_string+=";"+ level_field +"=" + str(level)
            
    return gff_string+"\n"

def parallel_get_methylation_level_mysql(bins,mc_type, samples,databases, servers, method = "weighted", use_all_reads=True, 
                                   processors=2,write_to_file=False,min_reads=0,as_matrix=False):
    
    """
    This function takes a list of lists [[chr,start,end],[chr,start,end],...]
    Remember that the list is a list of the methylation types you want returned so ["CHG","CHH"] will return
    the levels for CHG AND SEPARATELY CHH. So if you wanted the combined level of CHG and CHH you'd have to do CHN
    as_matrix means the output will be written in the form chr, start, end, type, sample1, sample2, sample3,...
    """
    if not isinstance(mc_type, list):
        if isinstance(mc_type, basestring):
            mc_type = [mc_type]
        else:
            sys.exit("mc_type must be a list")
    if not isinstance(samples, list):
        if isinstance(samples, basestring):
            samples = [samples]
        else:
            sys.exit("samples must be a list")
    if not isinstance(servers, list):
        if isinstance(servers, basestring):
            servers = [servers]
        else:
            sys.exit("servers must be a list")
    if not isinstance(databases, list):
        if isinstance(databases, basestring):
            databases = [databases]
        else:
            sys.exit("databases must be a list")
    
    if processors > 1:
        try:
            pool = Pool(processors)
            results = []
            for sample,database,server in zip(samples,databases,servers):
                conn = MySQLdb.connect(server,"mysql","rekce",database)
                cursor = conn.cursor()
                cursor.execute("select max(position) from allc_"+sample+"_"+bins[0][0])
                max_end = cursor.fetchone()[0]
                chunk_size = max(len(bins) / processors,1)
                for chunk_start in xrange(0,len(bins),chunk_size):
                    #print_checkpoint("\t\t"+str(bin_start))
                    results.append(pool.apply_async(get_mc_level_mysql,[bins[chunk_start:(chunk_start+chunk_size)],sample,database,server,mc_type,method,use_all_reads,min_reads]))
                cursor.close()        
                conn.close()
            pool.close()
            pool.join()
    
        except Exception, e:
            pdb.set_trace()
            print e
            pool.terminate()
            sys.exit("Killed Pool")

        result_dict = []
        for result in results:
            for bin in result.get():
                #numpy arrays may be more flexible for indexing
                result_dict.append(bin)
                """
                if bin[0] not in result_dict:
                    result_dict[bin[0]] = {}
                result_dict[bin[0]][(bin[1],bin[2],bin[3])] = bin[5]
            """
    else:
        results = []
        for sample,database,server in zip(samples,databases,servers):
            conn = MySQLdb.connect(server,"mysql","rekce",database)
            cursor = conn.cursor()
            cursor.execute("select max(position) from allc_"+sample+"_"+bins[0][0])
            max_end = cursor.fetchone()[0]
            chunk_size = len(bins) / processors
            for bin in bins:
                #print_checkpoint("\t\t"+str(bin_start))
                try:
                    results.append(get_mc_level_mysql([bin],sample,database,server,mc_type,method,use_all_reads,min_reads))
                except Exception, e:
                    print e
                    pdb.set_trace()
            cursor.close()        
            conn.close()
            result_dict = []
            for bin in results:
                result_dict.append(bin)
                """
                if bin[0] not in result_dict:
                    result_dict[bin[0]] = {}
                result_dict[bin[0]][(bin[1],bin[2],bin[3])] = bin[5]
                """
                
    array = np.array(result_dict,{"names":("chr","start","end","sample","type","mc","h","level"),"formats":("a2","i8","i8","a20","a3","i8","i8","a21")})            
    
    if write_to_file != False:
        g = open(write_to_file,'w')
        if as_matrix == True:
            coords = ""
            g.write("chr\tstart\tend\ttype\t"+"\t".join(samples)+"\n")
            for bin in np.sort(array,order=["chr","start","end","sample"]):
                if [bin["chr"],bin["start"],bin["end"],bin["type"]] != coords:
                    g.write("\n")
                    g.write("\t".join(map(str,[bin["chr"],bin["start"],bin["end"],bin["type"],bin["level"]])))
                    coords = [bin["chr"],bin["start"],bin["end"],bin["type"]]
                else:
                    g.write("\t"+bin["level"])
            
        else:
            for name in array.dtype.names:
                g.write(name+"\t")
            g.write("\n")
            for bin in np.sort(array,order=["chr","start","end","sample"]):
                for name in array.dtype.names:
                    g.write(str(bin[name])+"\t")
                g.write("\n")
                #g.write(sample+"\t"+"\t".join(map(str,bin))+"\t"+str(result_dict[sample][bin])+"\n")
        g.close()
    return array

def get_mc_level_mysql(bins,sample, database, server,mc_type,method="weighted",use_all_reads=True,min_reads = 0):
        conn = MySQLdb.connect(server,"mysql","rekce",database)
        cursor=conn.cursor()
        levels = []
        #pdb.set_trace()
        for chr,bin_start,bin_end in bins:
            #REMEMBER EACH TYPE IN THE LIST IS EVALUATED SEPARATELY AND RETURNED ON IT'S OWN. THIS SAVES TIME
            #BECAUSE THE LONG PART IS GETTING THE ROWS FROM MYSQL
            cursor.execute("select * from allc_"+sample+"_"+chr+" where position >= "+str(bin_start)+" and position <= "+str(bin_end))
            rows=cursor.fetchall()
            for type in mc_type:
                current_mc_types = [type]
                iub_dict = {"N":["A","C","G","T"],"H":["A","C","T"],"C":["C"],"G":["G"],"T":["T"],"A":["A"]}
                type += "N" * (3 - len(type))
                current_mc_types.extend(["".join(i) for i in itertools.product(*[iub_dict[nuc] for nuc in type])])
                
                if type == "C":
                    current_mc_types.extend(["CG","CHG","CHH","N"])
                numerator = 0
                denominator = 0
                for row in rows:
                    if row[5] < min_reads:
                        continue
                    if row[3] not in current_mc_types:
                        continue
                    if method == "weighted":
                        if row[6] == 0 and use_all_reads == True:
                            denominator += row[5]
                        if row[6] == 1:
                            numerator += row[4]
                            denominator += row[5]
                    elif method == "mattia":
                        sys.exit("MATTIA'S METRIC NOT IMPLEMENTED")
                    elif method == "average":
                        sys.exit("AVERAGE METRIC NOT IMPLEMENTED")
                    elif method== "count":
                        denominator = 1
                        if row[6] == 1:
                            numerator += 1
                    elif method == "fraction":
                        if row[6] == 1:
                            numerator += 1
                        denominator += 1
                try:
                    levels.append((chr,bin_start,bin_end,sample,type,numerator,denominator,float(numerator) / float(denominator)))
                except:
                    levels.append((chr,bin_start,bin_end,sample,type,numerator,denominator,"NA"))
        cursor.close()        
        conn.close()
        return levels
        
    
def parallel_get_methylation_level(files, mc_type, method="weighted", use_all_reads=True, start=False, end=False, 
                                   processors=2, receive_size=False, write_to_file=False):
    """
    Parallelized get_methylation_level. Input takes a list of gff files OR a single filename OR a list of gelists and 
    writes the results of get_methylation_level on each file in a new file named filename_methlevels.gff.
    
    Receive_size is the number of lines that you want to be fetched at a time (in a block). If it is too large, 
    then the program will warn you. If this is left false, then it tries to estimate the best receive size based on 
    the first line of your gff file.
    mc_type is the type of methylated C's you want to use in get_methylation_level
    method is the way used to get the methylation levels
    use_all_reads // start/end are whether you want to use every read
    write_to_file is the filename of where you want to save the gff result of get_methylation_level; if none->False.
    """
    MAX_PACKET_SIZE = 33900 #empirical, could not find actual size in documentation
    meth_level_results = []
    if not isinstance(files, list):
        files = [files]
        
    for filename in files:
        #start children
        all_conns = [] #list of child,parent connections for the pipes
        process = [] #list of pipe processes
        for proc in xrange(processors): #create all the child/parent pipes for processors
            parent_conn, child_conn = Pipe()
            process.append(Process(target=child_methylation_levels, 
                                   args=(child_conn, parent_conn, method, use_all_reads, start, end)))
            process[-1].start()
            child_conn.close()
            all_conns.append(parent_conn)
        
        if not isinstance(filename, gelist):
            gel = gelist(filename)
        else:
            gel = filename
        
        #estimate optimal packet size to fetch
        if receive_size==False:
            receive_size = MAX_PACKET_SIZE/sys.getsizeof(write_gelist_Line(gel[0])) 
        
        gel_index = 0       #used to keep track of where we are in the gelist
        processor_index = 0 #used to increment through the processors sending info
        parent_index = 0    #keeps track of the index for the processors that the parent is using to receive information
        line_counter = 0    #count # of lines sent (equal to # of lines you should receive back)
        packet_size_sum = 0 #keeps track of the size of the current number of packets that are currently being processed
        gelist_length = len(gel)
        result = ""
        
        #loop while there is gelist or while there are lines that havent been returned
        while gel_index < gelist_length or line_counter > 0: 

            packet_size_sum += sys.getsizeof(write_gelist_Line(gel[0]))  #add this line to the packet size

            if packet_size_sum > MAX_PACKET_SIZE: #if the packet size exceeds the limit, terminate
                for p in process:
                    p.terminate()
                sys.exit("Warning! The lines of your gff file were too long. You should decrease the packet receive size.")

            if gel_index < gelist_length: #IF there are lines to send, send them
                line = gel[gel_index] #get current line
                all_conns[processor_index].send((line, mc_type)) #send info to use get_mysql_data in child
                gel_index +=1
                line_counter += 1
            
            #go receive levels if within receive size or no more gelist
            if line_counter%receive_size==0 or gel_index>=gelist_length: 
                packet_size_sum = 0 #reset the packet size that is currently being processed
                while line_counter > 0: #while there are still lines to receive
                    line = all_conns[parent_index].recv() #get the line (string)
                    result += line 
                    line_counter -= 1 
                    parent_index += 1 #increment the processor so that the next loop receives from the correct processor
                    if parent_index >= processors:
                        parent_index = 0
            
            processor_index += 1 #increment index to next child
            if processor_index >= processors:
                processor_index = 0
        
        for parent_conn in all_conns:
            parent_conn.close()
        for p in process:
            p.join()
        subprocess.check_call(["rm", "-r", "tmp_results/"]) #rm temp folder
        if write_to_file!=False:
            f=open(write_to_file, 'w')
            f.write(result)
            f.close()
        else:    
            meth_level_results.append(result)
    
    return meth_level_results

def child_methylation_levels(child_conn, parent_conn, method, use_all_reads, start, end):
    """
    Child processes that are piped line data for get_meth_level. They run get_methylation_level and write the results to 
    a temporary file that is combined again in parallel_get_methylation level after the process is finished.
    """
    parent_conn.close() #close unused end of the pipe
    ensure_dir("tmp_results") #make temporary directory to put results in
    while True:
        try:
            line, mc_type = child_conn.recv() 
        except EOFError:
            child_conn.close() #close pipe when you reach end of file or other end is closed
            break
        line = write_gelist_Line(line, num_fields=9) #reproduce the line string from dictionary
        if use_all_reads==True:
            all_reads = "T"
        else:
            all_reads = "F"
        level_field = ','.join(mc_type) + "_" + method[0] + "_" + all_reads #field name = type,type_method_T/F
        level = get_methylation_level(line, mc_type, method, use_all_reads, start, end)
        line = line[:-1] + ";" + level_field + "=" + str(level) + "\n" #trim off newline of old string and add newline at end
        child_conn.send(line)

def get_methylation_level(element,mc_type,method="weighted",use_all_reads=True,start=False,end=False):
    """
    Expects an interval from a BedTool and returns the methylation level.
    mc_type is a list of the C contexts you want the level for (e.g., CHG). For example ["CG","CHG"] will only
    consider CG and CHG methylation for the methylation level. You may use Ns and Hs anywhere in the list.  
    method is one of weighted, simple, mattia
    -average indicates that you want to weight each mc site equally
    -weighted takes coverage into account by weighting the methylation level by the number of reads at the site
    -mattia calculates methylation level as (mc/h + mc/h + mc/h) / (basepairs)
    -count indicates that you want to use the number of mc sites in the window as the methylation level
    -fraction indicates that you want to return the fraction of Cs in the window that are methylated
    use_all_reads indicates whether or not you want to include unmethylated Cs in the calculation of the methylation level
    start indicates the beginning of a subregion of the element that you'd like the methylation level for
    end indicates the end of a subregion of the element that you'd like the methylation level for 
    NOTE: STACKS TABLES ARE ALWAYS USED TO CHECK FOR COVERAGE. IF THERE ARE NO STACKS POSITIONS IN A REGION AN NA 
    WILL GET RETURNED!
    """
    if method not in ["weighted", "average", "mattia", "count", "fraction"]:
        sys.exit("In get_methylation_level, method must be one of 'weighted', 'simple', 'average', 'mattia', \
                 'count', or 'fraction', not "+method)
    
    if not isinstance(mc_type, list):
        if isinstance(mc_type, basestring):
            mc_type = [mc_type]
        else:
            sys.exit("mc_type must be a list")
    
    if isinstance(element,str):
        element = parse_gff_line(element)
    
    if use_all_reads:
        all_reads = "T"
    else:
        all_reads = "F"
    level_field = ','.join(mc_type) + "_" + method[0] + "_" + all_reads
    if level_field in element:
        if element[level_field] == "NA":
            return "NA"
        else:
            return float(element[level_field])
    
    iub_dict = {"N":["A","C","G","T"],"H":["A","C","T"],"C":["C"],"G":["G"],"T":["T"],"A":["A"]}
    expanded_mc_type = mc_type[:]
    for type in mc_type:
        type += "N" * (3 - len(type))
        expanded_mc_type.extend(["".join(i) for i in itertools.product(*[iub_dict[nuc] for nuc in type])])
            
    if "C" in mc_type:
        expanded_mc_type.extend(["CG","CHG","CHH","N"])
    mc_type = expanded_mc_type
    if len(element['mc_classes']) == 0:
        mc_classes_list = []
    else:
        mc_classes_list = element['mc_classes'].split(',')
        
    if len(element['stacks_classes']) == 0:
        stacks_classes_list = []
    else:
        stacks_classes_list = element['stacks_classes'].split(',')
        
    if len(element['mc_positions']) == 0:
        mc_positions_list = []
    else:
        mc_positions_list = element['mc_positions'].split(',')
        
    if len(element['stacks_positions']) == 0:
        stacks_positions_list = []
    else:    
        stacks_positions_list = element['stacks_positions'].split(',')
    #Check if there are any elements 
    #These are lists of indexes that meet the start/end requirements as well as the mc class requirements
    mc_indexes = [index for index, value in enumerate(zip(mc_positions_list,mc_classes_list)) if
                  (start == False or int(value[0]) >= start) and (end == False or int(value[0]) <= end) and 
                  value[1] in mc_type]
    stacks_indexes = [index for index, value in enumerate(zip(stacks_positions_list,stacks_classes_list)) if 
                      (start == False or int(value[0]) >= start) and (end == False or int(value[0]) <= end) and 
                      value[1] in mc_type]
    #Check that there is coverage in the block at all
    if len(stacks_indexes) == 0:
        return "NA"
    if len(mc_indexes) == 0:
        return 0.0
    
    mc_counts_list = np.array([float(value) for value in element['mc_counts'].split(',')])[mc_indexes]
    h_counts_list = np.array([float(value) for value in element['h_counts'].split(',')])[mc_indexes]
    mc_strands_list = np.array(element['mc_strands'].split(','))[mc_indexes]    
    stacks_total_list = np.array([float(value) for value in element['stacks_total'].split(',')])[stacks_indexes]
    if method == "average":
        values = mc_counts_list / h_counts_list
        if use_all_reads == False:
            num_c = len(mc_counts_list)
        else:
            num_c = len(stacks_total_list)
        meth_level = sum(values) / float(num_c)
        return meth_level
    elif method == "mattia":
        #I implemented this just for legacy's sake. I really don't think anyone should use it...
        if element["_format"] == "gff":
            meth_level = sum((mc_counts_list / h_counts_list) / (int(element[4]) - int(element[3])))
        else:
            meth_level = sum((mc_counts_list / h_counts_list) / (int(element[2]) - int(element[1])))
        return meth_level
    elif method == "count":
        return len(mc_counts_list)
    elif method == "fraction":
        return float(len(mc_counts_list)) / float(len(stacks_total_list))
    #This is the weighted method
    else:
        if use_all_reads == False:
            meth_level = (sum(mc_counts_list) / sum(h_counts_list))
        else:
            meth_level = (sum(mc_counts_list) / sum(stacks_total_list))
        return meth_level
    
def create_matrix(gel, method, mc_type, file,use_all_reads=True):
    """
    This function will write out a tab separated file of the elements in the gelist
    as well as their methylation levels.
    method indicates which method of methylation level computation you'd like to use
    mc_type indicaes which methylation context (C,CG,CHG,CHH) you'd like the level for
    file indicates the file where you'd like the result saved
    """
    f = open(file,'w')
    #This helps reorder the matrix so that elements can be rows and samples columns
    temp_dict = {}
    for element in gel:
        sample = element[1]
        level = get_methylation_level(element,mc_type,method,use_all_reads)
        coord = (element[0],element[3],element[4])
        if coord not in temp_dict:
            temp_dict[coord] = {}
        temp_dict[coord][sample] = level
    
    f.write("chr\tstart\tend\t"+"\t".join(gel.samples)+"\n")
    
    for coord in sorted(temp_dict.keys(),cmp=sort_by_coord_tuple):
        f.write("\t".join([coord[0],coord[1],coord[2]]))
        for sample in gel.samples:
            try:
                value = temp_dict[coord][sample]
            except:
                value = 0
            f.write("\t"+str(value))
        f.write("\n")
    f.close()
    
def set_c_index(path):
    #This path indicates the directory where the c_index files are stored.  
    #They are expected to be in the format <sample>_c_index_chr<chrom num>
    global c_index_path
    
    if path[-1] == "/":
        sys.exit("Your c_index path ends in a /, which means you haven't entered this variable correctly. " \
                  + "It should include the entire c_index name up to the part specifiying the chromosome")
    else:
        c_index_path = path

def load_base_c_index(chromosomes):
    global base_c_index

    base_c_index = {}
    print_checkpoint("Begin loading base_c_index")
    for chr in chromosomes:
        try:
            file = open(c_index_path+"_"+chr,'r')
        except:
            print "Could not open c_index: "+c_index_path+"_"+chr
            return False
        base_c_index[chr] = cPickle.load(file)
        file.close()
    print_checkpoint("Done loading base_c_index")
    
def load_sample_c_index(chr,sample):

    sample = sample.replace("_bud","")
    try:
        file = open(c_index_path+"_"+sample+"_"+chr,'r')
    except:
        print "Could not open c_index: "+c_index_path+"_"+sample+"_"+chr
        return False
    c_index = cPickle.load(file)
    file.close()
    return c_index

if __name__ == '__main__':
    pass