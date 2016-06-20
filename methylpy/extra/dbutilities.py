import sys
try:
    import MySQLdb
except:
    sys.exit("methylpy.dbutilities requires the MySQLdb module")
try:
    from utilities import print_checkpoint
except:
    sys.exit("methylpy.dbutilities requires the methylpy.utilities module")
    
def create_allc_table(sample,chr,server,database,import_tables=False):
    """
    create_allc_table goes into the stacks tables and makes a new allc table that contains all positions where there is a "C".
    It fills in the type of C it is by checking the next two bases. It also determines whether the C is methylated or not and fills
    in the methylated column with a 0 or 1.
    
    """
    connection = MySQLdb.connect(server,"mysql","rekce",database)
    cursor = connection.cursor(cursorclass=MySQLdb.cursors.SSCursor)
    if import_tables == True:
        cursor.execute("drop table if exists allc_"+sample+"_"+chr.replace("chr",""))
        cursor.execute("CREATE TABLE allc_"+sample+"_"+chr.replace("chr","")+"( assembly char(2) NOT NULL," + \
                   " position INT(10) UNSIGNED NOT NULL DEFAULT 0, strand ENUM('+', '-', '.') NOT NULL DEFAULT '.'," + \
                   " class CHAR(3) NOT NULL DEFAULT 'NNN', mc SMALLINT(5) UNSIGNED NOT NULL DEFAULT 0, h SMALLINT(5)" + \
                   " UNSIGNED NOT NULL DEFAULT 0, methylated BOOL, PRIMARY KEY (position,strand));")
        if cursor.fetchall():
            sys.exit("Could not create mysql table")
    
    #go through all mc positions in mc table and add position,strand to set if it exists
    print_checkpoint("Getting mc information")
    mc_positions = set()
    mcquery = "select * from "+database+".mc_"+sample+"_"+chr.replace("chr","")
    g = open("mysql_mc_output_"+sample+"_"+chr+".tsv",'w')
    subprocess.check_call(shlex.split("mysql -u mysql --password=rekce -h "+server+" -e "+"'"+mcquery+"'"),stdout=g)
    g.close()
    f = open("mysql_mc_output_"+sample+"_"+chr+".tsv",'r')
    for line in f:
        fields = line.split("\t")
        mc_positions.add((str(fields[1]), str(fields[2])))
    f.close()
    #cursor.execute(mcquery)
    #rows = cursor.fetchall()
    #for line in rows:
    #    mc_positions.add((str(line[1]), str(line[2])))
    
    #A file to temporarily store results to be loaded into mysql
    print_checkpoint("Getting Rows")
    query = "select position,fcall,ftotal,fbase,rcall,rtotal,rbase,fc,rc from "+database+".stacks_bs_"+sample+"_"+chr.replace("chr","")
    #cursor.execute(query)
    #cursor.arraysize = 5000000
    #rows = cursor.fetchmany()
    g = open("mysql_stacks_output_"+sample+"_"+chr+".tsv",'w')
    subprocess.check_call(shlex.split("mysql -u mysql --password=rekce -h "+server+" -e "+"'"+query+"'"),stdout=g)
    g.close()
    f = open("mysql_stacks_output_"+sample+"_"+chr+".tsv",'r')
    g = open("allc_"+sample+"_"+chr+".tsv",'w')
    print_checkpoint("Finding C Contexts")
    line1 = f.readline().rstrip()
    line2 = f.readline().rstrip()
    line3 = f.readline().rstrip()
    line4 = f.readline().rstrip()
    line5 = f.readline().rstrip()
    while line3:
        fields1 = line1.split("\t")
        fields2 = line2.split("\t")
        fields3 = line3.split("\t")
        fields4 = line4.split("\t")
        fields5 = line5.split("\t")
        
        if int(fields3[2]) > 0 and fields3[3]=="C":
            allc_row = []
            allc_row.append(chr.replace("chr",""))  #append chromosome
            allc_row.append(str(fields3[0]))            #append position
            allc_row.append("+")                    #append strand
            try:
                c_class = "C"
                if int(fields4[0]) == int(fields3[0]) + 1:
                    c_class += fields4[3]
                else:
                    c_class += "N"
                if int(fields5[0]) == int(fields3[0]) + 2:
                    c_class += fields5[3]
                else:
                    c_class += "N"                  
            except:
                c_class = "N"
            allc_row.append(c_class)
            allc_row.append(str(fields3[7]))    #append mc
            allc_row.append(str(fields3[2]))    #append h
            if (allc_row[1], allc_row[2]) in mc_positions:
                allc_row.append("1")
            else:
                allc_row.append("0")
            g.write("\t".join(allc_row)+"\n")
            
        elif int(fields3[5]) > 0 and fields3[6]=="C":
            allc_row = []
            allc_row.append(chr.replace("chr",""))  #append chromosome
            allc_row.append(str(fields3[0]))            #append position
            allc_row.append("-")                    #append strand
            try:
                c_class = "C"
                if int(fields2[0]) == int(fields3[0]) - 1:
                    c_class += fields2[6]
                else:
                    c_class += "N"
                if int(fields1[0]) == int(fields3[0]) - 2:
                    c_class += fields1[6]
                else:
                    c_class += "N"                  
            except:
                c_class = "N"
            allc_row.append(c_class)        #append class
            allc_row.append(str(fields3[8]))    #append mc
            allc_row.append(str(fields3[5]))    #append h
            if (allc_row[1], allc_row[2]) in mc_positions:
                allc_row.append("1")
            else:
                allc_row.append("0")
            g.write("\t".join(allc_row)+"\n")
        
        line1 = line2
        line2 = line3
        line3 = line4
        line4 = line5
        line5 = f.readline().rstrip()
    f.close()
    g.close()
    subprocess.check_call(["rm", "mysql_stacks_output_"+sample+"_"+chr+".tsv"])
    subprocess.check_call(["rm", "mysql_mc_output_"+sample+"_"+chr+".tsv"])
    print_checkpoint("Done.")


    """
    code for mysql
    while rows:    
        #Duplicating this for loop helps avoid more checks in the if statements below
        for index in range(0,2):
            row = rows[index]
            if row[2] > 0 and row[3]=="C":
                allc_row = []
                allc_row.append(chr.replace("chr",""))  #append chromosome
                allc_row.append(str(row[0]))            #append position
                allc_row.append("+")                    #append strand
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
                allc_row.append(c_class)
                allc_row.append(str(row[7]))    #append mc
                allc_row.append(str(row[2]))    #append h
                if (allc_row[1], allc_row[2]) in mc_positions:
                    allc_row.append("1")
                else:
                    allc_row.append("0")
                g.write("\t".join(allc_row)+"\n")
           
        for index in range(2,len(rows)-2):
            row = rows[index]
            if row[2] > 0 and row[3]=="C":
                allc_row = []
                allc_row.append(chr.replace("chr",""))  #append chromosome
                allc_row.append(str(row[0]))            #append position
                allc_row.append("+")                    #append strand
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
                allc_row.append(c_class)        #append class
                allc_row.append(str(row[7]))    #append mc
                allc_row.append(str(row[2]))    #append h

    
            elif row[5] > 0 and row[6]=="C":
                allc_row = []
                allc_row.append(chr.replace("chr",""))  #append chromosome
                allc_row.append(str(row[0]))            #append position
                allc_row.append("-")                    #append strand
                try:
                    c_class = "C"
                    if rows[index-1][0] == rows[index][0] - 1:
                        c_class += rows[index-1][6]
                    else:
                        c_class += "N"
                    if rows[index-2][0] == rows[index][0] - 2:
                        c_class += rows[index-2][6]
                    else:
                        c_class += "N"                  
                except:
                    c_class = "N"
                allc_row.append(c_class)        #append class
                allc_row.append(str(row[8]))    #append mc
                allc_row.append(str(row[5]))    #append h
            else:
                continue
            if (allc_row[1], allc_row[2]) in mc_positions:
                allc_row.append("1")
            else:
                allc_row.append("0")
            g.write("\t".join(allc_row)+"\n")
    
        for index in range(len(rows)-2,len(rows)):
            row = rows[index]
            if row[5] > 0 and row[6]=="C":
                allc_row = []
                allc_row.append(chr.replace("chr",""))      #append chromosome
                allc_row.append(str(row[0]))                #append position
                allc_row.append("-")                        #append strand
                try:
                    c_class = "C"
                    if rows[index-1][0] == rows[index][0] - 1:
                        c_class += rows[index-1][6]
                    else:
                        c_class += "N"
                    if rows[index-2][0] == rows[index][0] - 2:
                        c_class += rows[index-2][6]
                    else:
                        c_class += "N"                  
                except:
                    c_class = "N"
                allc_row.append(c_class)        #append class
                allc_row.append(str(row[8]))    #append mc
                allc_row.append(str(row[5]))    #append h
                if (allc_row[1], allc_row[2]) in mc_positions:
                    allc_row.append("1")
                else:
                    allc_row.append("0")
                g.write("\t".join(allc_row)+"\n")
                
        next_rows = cursor.fetchmany()
        #If fetchmany returs rows than we have more things to process
        #otherwise rows will be set to an empty list and the while loop above will exit
        if len(next_rows) > 0:
            rows = [rows[-1],rows[-2]]
            rows.extend(next_rows)
        else:
            rows = []

    
    g.close()
    if import_tables == True:
        cursor.execute("LOAD DATA LOCAL INFILE 'allc_"+sample+"_"+chr+".tsv' into table "+"allc_"+sample+"_"+chr.replace("chr",""))
    cursor.close()
    connection.close()
    """

if __name__ == '__main__':
    pass