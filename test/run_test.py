#/usr/bin/env python
import sys
import subprocess
import shlex
verbose=True

def print_error(msg):
    sys.stdout.write(msg)
    sys.stdout.flush()
    exit()

def print_msg(message,print_out=True):
    if print_out:
        sys.stdout.write(message)
        sys.stdout.flush()

def test_methylpy_import(function_name,verbose):
    print_msg("- importing %s: " %(function_name),
              verbose)
    try:
        exec("from methylpy import "+function_name)
        print_msg("pass\n",verbose)
    except:
        print_msg("failed\nExit\n",verbose)
        exit()

def get_executable_version(exec_name):
    try:
        sys.stdout.write("- find %s: " %(exec_name))
        if exec_name != "java":
            out = subprocess.check_output(shlex.split(exec_name+" --version"))
        else:
            out = subprocess.check_output(shlex.split(exec_name+" -version"),
                                          stderr=subprocess.STDOUT)
            out = out.replace("\"","")
        first_line = out.split("\n")[0]
        fields = first_line.split(" ")
        print("found version %s" %(fields[-1]))
        return(True)
    except:
        sys.stdout.write("failed\n")
        return(False)


f_stdout = open("test_output_msg.txt",'w')
f_stderr = open("test_error_msg.txt",'w')
    
# 1 - Check methylpy installation
print("Test importing methylpy module.")
test_methylpy_import("call_mc_se",verbose)
test_methylpy_import("call_mc_pe",verbose)
test_methylpy_import("DMRfind",verbose)
test_methylpy_import("parser",verbose)
# executable
sys.stdout.write("Check methylpy executable:")
try:
    out = subprocess.check_output(shlex.split("methylpy -h"))
    sys.stdout.write("pass\n")
except:
    print_error("failed\n"
                +"Please check whether methylpy/bin/ is included in PATH.\n"
                +"Exit\n")
# successful
print("\nmethylpy is successfully installed!\n")


# 2 - Check dependencies
print("Check whether dependencies are available in PATH")
bowtie = get_executable_version("bowtie")
bowtie2 = get_executable_version("bowtie2")
samtools = get_executable_version("samtools")
cutadapt = get_executable_version("cutadapt")
java = get_executable_version("java")

path_to_samtools = '""'
if not samtools:
    print("")
    while not samtools:
        print("samtools is required but it is not detected.")
        path_to_samtools = str(raw_input("Please enter path to samtools: "))
        samtools = get_executable_version(path_to_samtools+"samtools")

path_to_bowtie, path_to_bowtie2 = '""','""'
if not bowtie and not bowtie2:
    print("")
    while not bowtie and not bowtie2:
        print("Neither bowtie nor bowtie2 is detected.")
        path_to_bowtie = str(raw_input("Please enter path to bowtie: "))
        path_to_bowtie2 = str(raw_input("Please enter path to bowtie2: "))
        bowtie = get_executable_version(path_to_bowtie+"/bowtie")
        bowtie2 = get_executable_version(path_to_bowtie2+"/bowtie2")
        
if bowtie and bowtie2:
    print("Both bowtie and bowtie2 are available and will be tested.")
elif bowtie:
    print("Only bowtie is available. Skip tests for bowtie2 related code")
elif bowtie2:
    print("Only bowtie2 is available. Skip tests for bowtie related code")

# 3 - Build reference
print("")
if bowtie:
    sys.stdout.write("Test build-reference with bowtie: ")
    subprocess.check_call(
        shlex.split("methylpy build-reference --input-files data/chrL.fa "
                    + "--output-prefix chrL/chrL --bowtie2 False "
                    +"--path-to-aligner "+path_to_bowtie),
        stdout=f_stdout,
        stderr=f_stderr)
    sys.stdout.write("pass\n")
if bowtie2:
    sys.stdout.write("Test build-reference with bowtie2: ")
    subprocess.check_call(
        shlex.split("methylpy build-reference --input-files data/chrL.fa "
                    + "--output-prefix chrL/chrL --bowtie2 True "
                    +"--path-to-aligner "+path_to_bowtie2),
        stdout=f_stdout,
        stderr=f_stderr)
    sys.stdout.write("pass\n")


# 4 - Single-end pipeline
print("")
if bowtie:
    sys.stdout.write("Test single-end-pipeline with bowtie: ")
    subprocess.check_call(
        shlex.split("methylpy single-end-pipeline "
                    +"--read-files data/test_data_R1.fastq "
                    +"--sample se_bt "
                    +"--path-to-output results/ "
                    +"--forward-ref chrL/chrL_f "
                    +"--reverse-ref chrL/chrL_r "
                    +"--ref-fasta data/chrL.fa "
                    +"--num-procs 1 "
                    +"--path-to-aligner "+path_to_bowtie+" "
                    +"--path-to-samtools "+path_to_samtools+" "
                    +"--remove-clonal False "
                    +"--binom-test True "
                    +"--bowtie2 False "
                    +"--unmethylated-control 0.005"),
    stdout=f_stdout,
    stderr=f_stderr)
    sys.stdout.write("pass\n")
if bowtie2:
    sys.stdout.write("Test single-end-pipeline with bowtie: ")
    subprocess.check_call(
        shlex.split("methylpy single-end-pipeline "
                    +"--read-files data/test_data_R1.fastq "
                    +"--sample se_bt2 "
                    +"--path-to-output results/ "
                    +"--forward-ref chrL/chrL_f "
                    +"--reverse-ref chrL/chrL_r "
                    +"--ref-fasta data/chrL.fa "
                    +"--num-procs 1 "
                    +"--path-to-aligner "+path_to_bowtie2+" "
                    +"--path-to-samtools "+path_to_samtools+" "
                    +"--remove-clonal False "
                    +"--binom-test True "
                    +"--bowtie2 True "
                    +"--unmethylated-control 0.005"),
        stdout=f_stdout,
        stderr=f_stderr)
    sys.stdout.write("pass\n")

# 5 - bam-quality-filter
sys.stdout.write("\nTest quality filter for BAM file of single-end data: ")
if bowtie:
    subprocess.check_call(
        shlex.split("methylpy bam-quality-filter "
                    +"--input-file results/se_bt_processed_reads.bam "
                    +"--output-file results/se_bt_processed_reads.filtered.bam "
                    +"--ref-fasta data/chrL.fa "
                    +"--path-to-samtools "+path_to_samtools+" "
                    +"--quality-cutoff 30 "
                    +"--min-num-ch 3 "
                    +"--max-mch-level 0.7 "
                    +"--buffer-line-number 100"),
        stdout=f_stdout,
        stderr=f_stderr)
if bowtie2:
    subprocess.check_call(
        shlex.split("methylpy bam-quality-filter "
                    +"--input-file results/se_bt2_processed_reads.bam "
                    +"--output-file results/se_bt2_processed_reads.filtered.bam "
                    +"--ref-fasta data/chrL.fa "
                    +"--path-to-samtools "+path_to_samtools+" "
                    +"--quality-cutoff 30 "
                    +"--min-num-ch 3 "
                    +"--max-mch-level 0.7 "
                    +"--buffer-line-number 100"),
        stdout=f_stdout,
        stderr=f_stderr)
sys.stdout.write("pass\n")

# 6 - Paired-end pipeline
print("")
if bowtie:
    sys.stdout.write("Test paired-end-pipeline with bowtie: ")
    subprocess.check_call(
        shlex.split("methylpy paired-end-pipeline "
                    +"--read1-files data/test_data_R1.fastq "
                    +"--read2-files data/test_data_R2.fastq "
                    +"--sample se_bt "
                    +"--path-to-output results/ "
                    +"--forward-ref chrL/chrL_f "
                    +"--reverse-ref chrL/chrL_r "
                    +"--ref-fasta data/chrL.fa "
                    +"--num-procs 1 "
                    +"--path-to-aligner "+path_to_bowtie+" "
                    +"--path-to-samtools "+path_to_samtools+" "
                    +"--remove-clonal False "
                    +"--binom-test True "
                    +"--bowtie2 False "
                    +"--unmethylated-control 0.005"),
    stdout=f_stdout,
    stderr=f_stderr)
    sys.stdout.write("pass\n")
if bowtie2:
    sys.stdout.write("Test paired-end-pipeline with bowtie: ")
    subprocess.check_call(
        shlex.split("methylpy paired-end-pipeline "
                    +"--read1-files data/test_data_R1.fastq "
                    +"--read2-files data/test_data_R2.fastq "
                    +"--sample se_bt2 "
                    +"--path-to-output results/ "
                    +"--forward-ref chrL/chrL_f "
                    +"--reverse-ref chrL/chrL_r "
                    +"--ref-fasta data/chrL.fa "
                    +"--num-procs 1 "
                    +"--path-to-aligner "+path_to_bowtie2+" "
                    +"--path-to-samtools "+path_to_samtools+" "
                    +"--remove-clonal False "
                    +"--binom-test True "
                    +"--bowtie2 True "
                    +"--unmethylated-control 0.005"),
        stdout=f_stdout,
        stderr=f_stderr)
    sys.stdout.write("pass\n")

# successful
print("\nAll tests are done!\n")
