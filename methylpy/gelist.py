'''
Created on Dec 16, 2011

@author: Matteo
'''

import sys
import pdb
try:
    import itertools as it
except:
    sys.exit("methylpy.gelist requires the itertools module")
try:
    from methylpy.utilities import print_checkpoint, parse_gff_line
except:
    sys.exit("methylpy.gelist requires the methylpy.utilities module")
try:
    import MySQLdb
except:
    sys.exit("methylpy.gelist requires the MySQLdb module")
try:
    from collections import defaultdict
except:
    sys.exit("methylpy.gelist requires the collections module")    
try:
    import linecache
except:
    sys.exit("methylpy.gelist requires the linecache module")
try:
    import cStringIO
except:
    sys.exit("methylpy.gelist requires the cStringIO module")

class gelist:
        
    def __init__(self, filename, from_string=False, check_mc_values=False, load_values=False, format=False):
        """
        Initiates a gelist object by filling in all of the appropriate dictionaries of information on where to
        find certain coordinates, samples, features, or lines. A gelist also has a string indicating the format of 
        the file backing the gelist, a set of samples, a set of chromosomes, and a set of features that are present 
        in this gelist.
        """
        if from_string==True:
            self.file = cStringIO.StringIO(filename)        
        else:
            try:
                self.file = open(filename, 'r')
            except:
                sys.exit("Cannot open file "+filename)
        
        self.chromosomes = set()       #A set of the chromosomes in the pybedtool
        self.samples = set()           #A set of the samples in the pybedtool
        self.coord_to_index = {}       #A hash of coords to their indexes in the pybedtool
        self.index_to_coord = {}       #hash of indexes in the pybedtool to the coordinates
        self.features = set()          #features present in this gelist
        self.format = ""               #format of the file backing the gelist
        self.values = {}               #hash of coords -> a list of lines that have those coords
        self.line_to_byte = {}         #hash of line number to bytes
        self.sample_to_byte = {}       #hash of sample to bytes
        self.feature_to_byte = {}      #hash of feature to bytes
        
        #find file type by first using extension, then number of fields. gff files have 9 fields
        if format!=False:
            self.format=format
        elif (filename.rfind(".")==-1 and from_string==False) or from_string==True:
            line = self.file.readline()
            if line.find("\t") != -1: #check to see if it is tab-separated (most often)
                line_length = len(line.split("\t"))
            else:
                line_length = len(line.split(" ")) #if it is not tab-separated, split by spaces
            if line_length==9:
                self.format = "gff"
            else:
                self.format = "bed"
        else:
            self.format = filename.split(".")[-1]
        byte = 0    #byte index
        linenum = 0  #line number
        self.file.seek(0) #go back to beginning
        while True: #do not use file iterator, f.tell() doesn't work with iterator
            line = self.file.readline()
            if not line:
                break
            
            line = parse_gff_line(line, format=self.format) #parse line
            self.chromosomes.add(line[0])
            self.line_to_byte[linenum] = byte
            #find format of file
            if self.format=="gff":
                self.samples.add(line[1])
                if line[1] not in self.sample_to_byte:
                    self.sample_to_byte[line[1]] = []
                self.sample_to_byte[line[1]].append(byte)
                
                self.features.add(line[2])
                if line[2] not in self.feature_to_byte:
                    self.feature_to_byte[line[2]] = []
                self.feature_to_byte[line[2]].append(byte)
                
                coord = (line[0], line[3], line[4])
            elif self.format=="bed":
                coord = (line[0], line[1], line[2])
            else:
                sys.exit("Unrecognized file argument of the gelist constructor")
            if coord not in self.coord_to_index:
                self.coord_to_index[coord] = []
            #add values to hashes
            self.coord_to_index[coord].append(byte)
            self.index_to_coord[byte]= coord
            
            if check_mc_values == True:
                #When you create a gelist double check that all the mC information is within a features coordinates
                #delete positions that are outside of it                 
                #Check for mc_positions field
                if "mc_positions" not in line:
                    sys.exit("Missing mc_positions field in: " + str(line))
                #This variable indicates whether or not mc_positions or stacks_positions are wrong
                changed = False
                mc_positions = line["mc_positions"].split(',')
                for index,position in enumerate(mc_positions[:]):
                    if int(position) < int(coord[1]) or int(position) > int(coord[2]):
                        changed = True
                if changed == True:
                    sys.exit("mc values are incorrect on this line: " + str(line))

                if "stacks_positions" not in line:
                    sys.exit("Missing stacks_positions field in: " + str(line))
                changed = False
                stacks_positions = line["stacks_positions"].split(',')
                for index,position in enumerate(stacks_positions[:]):
                    if int(position) < int(coord[1]) or int(position) > int(coord[2]):
                        changed = True
                if changed == True:
                    sys.exit("stacks values are incorrect on this line: " + str(line))
                            
            #put in memory if necessary
            if load_values == True:
                if coord not in self.values:
                    self.values[coord] = []
                self.values[coord].append(line)
            byte = self.file.tell() 
            linenum += 1        
        
    def __iter__(self):
        """
        Returns an iterator pointed at the beginning of the gelist file.
        """
        self.file.seek(0)
        return self
    
    def next(self):
        """
        Parses the next line of the gelist and returns it. Used for the iterator.
        """
        line = self.file.readline()
        if not line:
            raise StopIteration
        else:
            return parse_gff_line(line, format=self.format)
    
    def __contains__(self, i):
        """
        Tests whether a coordinate 'i' exists in the gelist.
        """
        if len(i)!=3:
            sys.exit("To check if a coordinate exists, you must provide chromosome, start, and stop")
        elif isinstance(i[0], int): #check if they input chromosome as integer only
            i = ('chr'+str(i[0]), i[1], i[2])
        elif i[0].find('chr')==-1: #check if they input chromosome as string without 'chr' in it
            i = ('chr'+i[0], i[1], i[2])

        return ((str(i[0]), str(int(i[1])), str(int(i[2]))) in self.coord_to_index)
    
    def __getitem__(self, i):
        """
        Index into the gelist using the bracket operator. Finds the byte value of the appropriate line and returns the parsed
        dictionary. Slicing of the gelist is also supported.
        
        You can filter by specifying your arguments in the brackets gelist[("filter_type", arguments)]
        Sample filter requires a list of samples.
        Chromosome filter requires a list of chromosomes.
        Feature filter requires a list of features.
        Coordinate filter requires chrom, start, stop.
        Range filter requires at least chrom, start, stop, and there are 3 optional arguments: contained, fraction_query, or
            fraction_subject. If fraction_query or fraction_subject is chosen, then the next argument must be the decimal that
            you want to limit the filter by.
            ie -
            gelist[("range", chromosome, start, stop, [contained, fraction_query, fraction_subject], [fraction decimal])]
        """
        if isinstance(i, slice):
            indices = i.indices(len(self)) #returns (start, end, step) tuple
            return self.slice(indices)
        elif isinstance(i, tuple):
            if i[0]=="sample":
                if len(i)!=2 or not isinstance(i[1],list):
                    sys.exit("Filtering by sample requires 1 argument: a list of samples")
                else:
                    return self.sample_filter(i[1])
            elif i[0]=="chromosome":
                if len(i)!=2 or not isinstance(i[1],list):
                    sys.exit("Filtering by chromosome requires 1 argument: a list of chromosomes")
                else:
                    return self.chr_filter(i[1])                
            elif i[0]=="feature":
                if len(i)!=2 or not isinstance(i[1],list):
                    sys.exit("Filtering by feature requires 1 argument: a list of features")
                else:
                    return self.feature_filter(i[1])                
            elif i[0]=="coordinate":
                if len(i)!=4:
                    sys.exit("Filtering by coordinate requires 3 arguments: chromosome, start, stop")
                else:
                    return self.coord_filter(i[1], i[2], i[3])                
            elif i[0]=="range":
                if len(i)<4 or len(i)>6:
                    sys.exit("Filtering by range requires at least 3 arguments: chromosome, start, stop.")
                else:
                    if len(i)>4:
                        if i[4]=="contained":
                            return self.range_filter(i[1], i[2], i[3], contained=True)  
                        elif i[4]=="fraction_query":
                             return self.range_filter(i[1], i[2], i[3], fraction_query=i[5])  
                        elif i[4]=="fraction_subject":
                            return self.range_filter(i[1], i[2], i[3], fraction_subject=i[5])  
                    else:
                        return self.range_filter(i[1], i[2], i[3])  
        else:
            byte = self.line_to_byte[i]
            self.file.seek(byte)
            line = self.file.readline()
            if not line:
                raise IndexError
            else:
                return parse_gff_line(line, format=self.format) 
    
    def slice(self, indices):
        """
        Helper function for the slicing of the gelist. Indices is a list of (start, end, step). 
        """
        for index in xrange(*indices): #loop over range of indices
            yield self[index]
    
    def __setitem__(self, i):
        sys.exit("You cannot set items of a gelist")
        
    def __delitem__(self, i):
        sys.exit("You cannot delete items of a gelist") 
    
    def __len__(self):
        """
        Returns the length of the gelist (number of lines present in the original file).
        """
        return len(self.line_to_byte.keys())
    
    def sample_filter(self, samples):
        """
        Filtering functions:
        
        Sample filter returns parsed dictionaries of the relevant lines from the samples query.
        Samples is a list of strings that contain the sample names that you want to look for.
        """
        sample_list = []
        for sample in samples:
            if sample not in self.samples:
                print "Warning: " + str(sample) + " does not exist here"
            else:
                sample_list.append(sample)
        for sample in sample_list:
            bytes = self.sample_to_byte[sample] #list of byte values where the lines are
            for byte in sorted(bytes):
                self.file.seek(byte)
                line = self.file.readline()
                if not line:
                    raise IndexError
                else:
                    yield parse_gff_line(line, format=self.format) 
    
    def chr_filter(self, chrom):
        """
        chr_filter filters the gelist by the list of strings containing chromosomes that is given to it.
        The list chrom must contain strings in the same format as was in the gelist, ie 'chr1' or 'ctg1', not just the number.
        """
        chr_list = []
        for chr in chrom:
            if chr not in self.chromosomes:
                print "Warning: " + str(chr) + " does not exist here"
            else:
                chr_list.append(chr)
        bytes = [self.coord_to_index[a,b,c] for a,b,c in self.coord_to_index if a in chr_list]
        #bytes is a list of lists; each list probably contains one byte inside, but just in case the function can 
        #handle more than 1 byte in each list
        for byte in sorted(bytes): 
            if len(byte)==1:
                self.file.seek(byte[0])
                line = self.file.readline()
                if not line:
                    raise IndexError
                else:
                    yield parse_gff_line(line, format=self.format) 
            else:
                for b in byte:
                    self.file.seek(b)
                    line = self.file.readline()
                    if not line:
                        raise IndexError
                    else:
                        yield parse_gff_line(line, format=self.format) 
    
    def feature_filter(self, features):
        """
        feature_filter filters the gelist by the third field in gff format and returns relevant dictionaries.
        """
        feature_list = []
        for feature in features:
            if feature not in self.features:
                print "Warning: " + str(feature) + " does not exist here"
            else:
                feature_list.append(feature) #list of wanted features that appear in this gelist
        for feature in feature_list:
            bytes = self.feature_to_byte[feature] #list of bytes where the wanted features can be found
            for byte in sorted(bytes):
                self.file.seek(byte)
                line = self.file.readline()
                if not line:
                    raise IndexError
                else:
                    yield parse_gff_line(line, format=self.format) 
    
    def coord_filter(self, chrom, start, end):
        """
        coord_filter returns subsets of the gelist that contain the exact match of chromosome, start, and end that is 
        given to the function. Start and end can be integers or strings.
        """
        try:
            bytes = self.coord_to_index[(chrom,str(start),str(end))] #see if coordinate exists
        except KeyError:                
            print "Warning: coordinate " + str(chrom) + ", start: " + str(start) + ", end: " + str(end) + \
                  ", does not exist here"
            return
        for byte in sorted(bytes): #go to line where coordinate exists
            self.file.seek(byte)
            line = self.file.readline()
            if not line:
                raise IndexError
            else:
                yield parse_gff_line(line, format=self.format) 
    
    def range_filter(self, chrom, start, end, contained=False, fraction_query=False, fraction_subject=False):
        """
        range_filter returns subsets of the gelist that match whose interval start and stop fit in the specified parameters.
        
        The default option requires only 1bp to overlap with the interval for it to be returned.
        The contained option forces that start and end of any returned subset to be completely contained within the query range.
        The fraction_query option calculates the fraction of overlapping bp / total number of bp in the query to return or not.
        The fraction_subject option calculated the fraction of overlapping bp/ total number of bp in the interval it wants to
            return to determine if that interval should be included. 
        """
        if contained!=False: #the interval must be completely contained within the query
            bytes = [self.coord_to_index[a,b,c] for a,b,c in self.coord_to_index if 
                     a==chrom and int(b)>=int(start) and int(c)<=int(end)]
        elif fraction_query!=False: #the fraction of overlap/query length must be greater than parameter
            bytes = [self.coord_to_index[a,b,c] for a,b,c in self.coord_to_index if a==chrom and 
                     (getOverlap((float(start),float(end)), 
                                 (float(b),float(c))) / (float(end)-float(start))) >= float(fraction_query)]
        elif fraction_subject!=False: #the fraction of overlap/gelist interval length must be greater than parameter
            bytes = [self.coord_to_index[a,b,c] for a,b,c in self.coord_to_index if a==chrom and 
                     (getOverlap((float(start),float(end)), 
                                 (float(b),float(c))) / (float(c)-float(b))) >= float(fraction_subject)]            
        else: #default; there must be at least 1 bp overlap
            bytes = [self.coord_to_index[a,b,c] for a,b,c in self.coord_to_index if 
                     a==chrom and (int(b)>=int(start) or int(c)<=int(end))]

        #bytes is a lists of lists of bytes. Each list should contain 1 byte, but there may be more     
        for byte in sorted(bytes):
            if len(byte)==1:
                self.file.seek(byte[0])
                line = self.file.readline()
                if not line:
                    raise IndexError
                else:
                    yield parse_gff_line(line, format=self.format) 
            else:
                for b in byte:
                    self.file.seek(b)
                    line = self.file.readline()
                    if not line:
                        raise IndexError
                    else:
                        yield parse_gff_line(line, format=self.format) 
            
    def saveas(self, name):
        """
        saveas writes the information from the gelist to the current directory
        """   
        name = str(name) + "." + self.format
        f = open(name, 'w')
        f.write(self.file.getvalue())
        f.close()

def write_gelist_Line(line, num_fields=9):
    """
    Formats a gelist interval dictionary as a string in gff format
    If you only want a certain number of fields, specify in num_fields variable.
    """
    if not isinstance(line, dict):
        sys.exit("write_gelist_Line takes a parsed gelist line as input")
        
    result_string = ""
    for field in xrange(num_fields):
        result_string += line[field] + "\t"
    result_string = result_string[:-1] #trim extra tab at end of string
    result_string += "\n" #add newline

    return result_string

def getOverlap(a, b):
    """
    getOverlap returns the overlap between the intervals a and b. Inclusive.
    """
    return max(0, 1 + min(a[1], b[1]) - max(a[0], b[0]))

def sort_by_coord_tuple(coord1,coord2):
    """
    Sorts the two coordinate tuples.
    """
    #expects tuple of (chrom,start,stop)
    #if the chromosomes are different return is based on chromosome
    lettertonum = {}
    lettertonum["X"] = "23"
    lettertonum["Y"] = "24"
    lettertonum["M"] = "25"
    lettertonum["C"] = "26"
    lettertonum["L"] = "27"
    if coord1[0] != coord2[0]:
        if coord1[0][3:] in lettertonum:
            coord1_chr = lettertonum[coord1[0][3:]]
        else:
            coord1_chr = coord1[0][3:]
        if coord2[0][3:] in lettertonum:
            coord2_chr = lettertonum[coord2[0][3:]]
        else:
            coord2_chr = coord2[0][3:]
        return int(coord1_chr) - int(coord2_chr)
    elif coord1[1] != coord2[1]:
        return int(coord1[1]) - int(coord2[1])
    else:
        return int(coord1[2]) - int(coord2[2])

def sort_by_coord_gel(interval1, interval2):
    """
    sort_by_coord_gel is a function used in the python comparator functions in order for it to sort the intervals
    from the gelist.
    """
    try:
        #expects coord1 and coord2 to be intervals from pybed tool
        lettertonum = {}
        lettertonum["X"] = "23"
        lettertonum["Y"] = "24"
        lettertonum["M"] = "25"
        lettertonum["C"] = "26"
        lettertonum["L"] = "27"
        if interval1["_format"] == "gff":
            coord1 = [interval1[0],interval1[3], interval1[4]]
        else:
            coord1 = [interval1[0],interval1[1], interval1[2]]
        if interval2["_format"] == "gff":
            coord2 = [interval2[0],interval2[3], interval2[4]]
        else:
            coord2 = [interval2[0],interval2[1], interval2[2]]
        
        if coord1[0] != coord2[0]:
            if coord1[0][-1] in lettertonum:
                coord1[0] = coord1[0][:-1] + lettertonum[coord1[0][-1]]
            if coord2[0][-1] in lettertonum:
                coord2[0] = coord2[0][:-1] + lettertonum[coord2[0][-1]]
                
            return int(coord1[0][3:]) - int(coord2[0][3:])
        elif coord1[1] != coord2[1]:
            return int(coord1[1]) - int(coord2[1])
        else:
            return int(coord1[2]) - int(coord2[2])
    except:
        pdb.set_trace()

def sort_by_sample_gel(interval1, interval2):
    """
    sort_by_sample_gel is a function used in the python comparator functions in order to sort intervals from gelist 
    first by sample, then by start of the interval, then by end interval
    
    APPLIES TO GFF FORMAT ONLY
    """
    lettertonum = {}
    lettertonum["X"] = "23"
    lettertonum["Y"] = "24"
    lettertonum["M"] = "25"
    lettertonum["C"] = "26"
    lettertonum["L"] = "27"
    if interval1["_format"] != "gff" or interval2["_format"] != "gff":
        sys.exit("sort_by_sample is only usable on GFF format intervals")
        
    sample1 = interval1[1]
    sample2 = interval2[1]
    coord1 = [interval1[0],interval1[3], interval1[4]]
    coord2 = [interval2[0],interval2[3], interval2[4]]
    
    if sample1 != sample2:
        if sample1 > sample2:
            return 1
        return -1
    elif coord1[0] != coord2[0]:
        if coord1[0][-1] in lettertonum:
            coord1[0] = coord1[0][:-1] + lettertonum[coord1[0][-1]]
        if coord2[0][-1] in lettertonum:
            coord2[0] = coord2[0][:-1] + lettertonum[coord2[0][-1]]
            
        return int(coord1[0][3:]) - int(coord2[0][3:])
    elif coord1[1] != coord2[1]:
        return int(coord1[1]) - int(coord2[1])
    else:
        return int(coord1[2]) - int(coord2[2])
        
#THE UPDATE FUNCTIONS ALLOW YOU TO CHANGE ASPECTS OF THE GFF FILEDS.
#PYBEDTOOLS DON'T SEEM TO NORMALLY ALLOW FOR ITEM ASSIGNMENT
#FOR EXAMPLE I CAN'T UPDATE THE START POSITION OF AN INTERVAL EXCEPT
#THROUGH METHODS LIKE THESE
def update_source(feature):
    """
    NEEDS TO BE COMPLETED DON'T USE!!
    This function changes the source (2nd column, index 1) of a pybedtool from a gff 
    """
    if feature.file_type != "gff":
        return feature
    
if __name__ == '__main__':
    pass