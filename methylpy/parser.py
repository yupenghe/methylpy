try:
    from argparse import ArgumentParser
except Exception,e:
    exc_type, exc_obj, exc_tb = exc_info()
    print(exc_type, exc_tb.tb_lineno)
    print e
    exit("methylpy.parser requires ArgumentParser from the argparse module")

from methylpy.DMRfind import DMRfind
from methylpy.call_mc_se import run_methylation_pipeline
from methylpy.call_mc_pe import run_methylation_pipeline_pe

def parse_args():
     # create the top-level parser
     parser = ArgumentParser(prog='PROG')
     subparsers = parser.add_subparsers(help='Process all commands', dest='command')

     add_DMRfind_subparser(subparsers)
     add_se_pipeline_subparser(subparsers)
     #add_DMRfind_parser(parser)
     #add_DMRfind_parser(parser)
     args = parser.parse_args()
     
     return(0)
 
     if args.command == "DMRfind":


         DMRfind(args.mc_type, reg_dict, args.samples, args.path_to_allc, args.num_procs, args.save_result,
                 args.min_cov, args.keep_temp_files, args.mc_max_dist, args.dmr_max_dist, args.resid_cutoff,
                 args.sig_cutoff,
                 args.num_sims,
                 args.min_tests, args.seed, args.min_num_dms, args.collapse_samples, 
                 args.sample_category, args.min_cluster)

     args = parser.parse_args()

     if args.command == "run_methylation_pipeline":
         if not args.aligner_options:
              args.aligner_options = ["-S","-k 1","-m 1","--chunkmbs 3072","--best","--strata","-o 4","-e 80","-l 20","-n 0"]
              
         run_methylation_pipeline(args.files,args.libraries,args.sample,args.forward_ref,args.reverse_ref,args.ref_fasta,
                                  args.unmethylated_control,args.path_to_samtools,args.path_to_aligner,
                                  args.aligner_options,args.num_procs,args.trim_reads,
                                  args.path_to_cutadapt,args.adapter_seq,args.max_adapter_removal,
                                  args.overlap_length,args.zero_cap,args.error_rate,args.min_qual_score,args.min_read_len,
                                  args.sig_cutoff,args.min_cov,args.binom_test,args.keep_temp_files,
                                  args.bowtie2,args.sort_mem,args.path_to_output)
         
     elif args.command == "call_methylated_sites":
          call_methylated_sites(args.inputf, args.sample, args.reference, args.control, args.casava_version, args.sig_cutoff,
                                args.num_procs, args.min_cov, args.binom_test, args.min_mc, args.path_to_samtools, args.sort_mem,
                                args.bh, args.path_to_files)

         

def add_DMRfind_subparser(subparsers):
     # create the parser for the "DMRfind" command
     parser_dmrfind = subparsers.add_parser('DMRfind', help='Use to run the DMRfind function')
     parser_dmrfind.add_argument('--mc_type', type=str, nargs="+", required=True, help="List of space separated mc nucleotide contexts for which you want to look for DMRs. These classifications may use the wildcards H (indicating anything but a G) and N (indicating any nucleotide)")
     parser_dmrfind.add_argument('--region', type=str, nargs="+", required=True, help="Space separated listing of chromosome, start, and end. The list elements should be positive integers. chr1 start end chr2 start end")
     parser_dmrfind.add_argument('--samples', type=str, nargs="+", required=True, help='List of space separated samples')
     parser_dmrfind.add_argument('--path_to_allc', type=str, required=True, help="String indicating the beginning of tab separated files containing methylation information for all C nucleotides in the genome")
     parser_dmrfind.add_argument('--num_procs', type=int, default=1, help='Number of processors you wish to use to parallelize this function')
     parser_dmrfind.add_argument('--save_result', type=str, default="temp", help='String indicating the prefix for result files')
     parser_dmrfind.add_argument('--min_cov', type=int, default=0, help='Minimum number of reads that must cover a site for it to be considered')
     parser_dmrfind.add_argument('--keep_temp_files', type=bool, default=False, help='Boolean; keep intermediate files?')
     parser_dmrfind.add_argument('--mc_max_dist', type=int, default=0, help='Integer indicating the maximum distance two sites can be from one another for their methylation counts to be combined. This option helps with low coverage experiments where you may want to leverage the correlation of methylation between sites to get more statistical power.')
     parser_dmrfind.add_argument('--dmr_max_dist', type=int, default=100, help='Maximum distance two significant sites can be to be included in the same block')
     parser_dmrfind.add_argument('--resid_cutoff', type=int, default=2, help='Results will have to show deviations in the contingency table in the same direction as the rest of the window')
     parser_dmrfind.add_argument('--sig_cutoff', type=float, default=.01, help='Float indicating at what FDR you want to consider a result significant')
     parser_dmrfind.add_argument('--num_sims', type=int, default=10000, help="Number of permutation tests you'd like to run to estimate the p-values of the differential methylation tests")
     parser_dmrfind.add_argument('--min_tests', type=int, default=100, help="Minimum number of permuation tests you'd like to run for each mC")
     parser_dmrfind.add_argument('--seed', type=int, default=-1, help='A seed to provide to the random number generator for permutation testing. Only change this if you are debugging and want to make sure the permutation output is consistent')
     parser_dmrfind.add_argument('--min_num_dms', type=int, default=0, help='The minimum number of differentially methylated sites that a differentially methylated region needs to contain to be reported')
     parser_dmrfind.add_argument('--collapse_samples', type=str, nargs='+', default=False, help='A list of samples for collapsing blocks')
     parser_dmrfind.add_argument('--sample_category', type=int, nargs='+', default=False, help='A list of categories that each respective sample belongs to; the categories must begin at 0 and increase by 1 for each category added. ex: samples [A,B,C] categories [0,1,2] or categories [0, 1, 0] ')
     parser_dmrfind.add_argument('--min_cluster', type=int, default=0, help='The minimum number of each sample category that must be present in every block that is output.')

def add_se_pipeline_subparser(subparsers):
    
     # create the parser for the "call_mc" command
     parser_pipeline = subparsers.add_parser('run_methylation_pipeline', help='Use to run the methylation pipeline')
     parser_pipeline.add_argument('--files', type=str, nargs="+", required=True, help="list of all the fastq files you'd like to run \
        through the pipeline. Note that globbing is supported here (i.e., you can use * in your paths)")
     parser_pipeline.add_argument('--libraries', type=str, nargs="+", required=True, help="list of library IDs (in the same order as \
        the files list) indiciating which libraries each set of fastq files belong to. If you use a glob, you only need to indicate \
        the library ID for those fastqs once (i.e., the length of files and libraries should be the same)")
     parser_pipeline.add_argument('--sample', type=str, required=True, help="String indicating the name of the sample you're processing. \
        It will be included in the output files.")
     parser_pipeline.add_argument('--forward_ref', type=str, required=True, help="string indicating the path to the forward strand \
        reference created by build_ref")
     parser_pipeline.add_argument('--reverse_ref', type=str, required=True, help="string indicating the path to the reverse strand \
        reference created by build_ref")
     parser_pipeline.add_argument('--ref_fasta', type=str, required=True, help="string indicating the path to a fasta file containing \
        the sequences you used for mapping")
     parser_pipeline.add_argument('--unmethylated_control', type=str, required=True, help="name of the chromosome/region that you want \
        to use to estimate the non-conversion rate of your sample, or the non-conversion rate you'd like to use. Consequently, control \
        is either a string, or a decimal. If control is a string then it should be in the following format: 'chrom:start-end'. \
        If you'd like to specify an entire chromosome simply use 'chrom:'")
     parser_pipeline.add_argument('--path_to_samtools', type=str, default="", help='Path to samtools installation (default is current dir)')     
     parser_pipeline.add_argument('--path_to_aligner', type=str, default="", help='Path to bowtie installation (default is current dir)')
     parser_pipeline.add_argument('--aligner_options', type=str, nargs='+', help="list of strings indicating options you'd like passed to bowtie \
        (e.g., ['-k 1','-l 2']")
     parser_pipeline.add_argument('--num_procs', type=int, default=1, help='Number of processors you wish to use to \
        parallelize this function')
     parser_pipeline.add_argument('--trim_reads', type=bool, default=True, help='Whether to trim reads using cutadapt (default is True)')
     parser_pipeline.add_argument('--path_to_cutadapt', type=str, default="", help='Path to cutadapt installation (default is current dir)')
     parser_pipeline.add_argument('--adapter_seq', type=str, default="AGATCGGAAGAGCACACGTCTG", help="sequence of an adapter that was ligated \
        to the 3' end. The adapter itself and anything that follows is trimmed.")
     parser_pipeline.add_argument('--max_adapter_removal', type=int, help="Indicates the maximum number of times to try to remove adapters. \
        Useful when an adapter gets appended multiple times.")
     parser_pipeline.add_argument('--overlap_length', type=int, help="Minimum overlap length. If the overlap between the read and the adapter \
        is shorter than LENGTH, the read is not modified. This reduces the no. of bases trimmed purely due to short random adapter matches.")
     parser_pipeline.add_argument('--zero_cap', type=bool, help="Flag that causes negative quality values to be set to zero (workaround to avoid \
        segmentation faults in BWA)")
     parser_pipeline.add_argument('--error_rate', type=float, help="maximum allowed error rate (no. of errors divided by the length \
        of the matching region)")
     parser_pipeline.add_argument('--min_qual_score', type=int, default=10, help="allows you to trim low-quality ends from reads before \
        adapter removal. The algorithm is the same as the one used by BWA (Subtract CUTOFF from all qualities; compute partial sums from \
        all indices to the end of the sequence; cut sequence at the index at which the sum is minimal).")
     parser_pipeline.add_argument('--min_read_len', type=int, default=30, help="indicates the minimum length a read must be to be kept. \
        Reads that are too short even before adapter removal are also discarded. In colorspace, an initial primer is not counted.")
     parser_pipeline.add_argument('--sig_cutoff', type=float, default=.01, help="float indicating the adjusted p-value cutoff you wish to \
        use for determining whether or not a site is methylated")
     parser_pipeline.add_argument('--min_cov', type=int, default=0, help="integer indicating the minimum number of reads for a site to be tested.")
     parser_pipeline.add_argument('--binom_test', type=bool, default=False, help="Indicates that you'd like to use a binomial test, rather than the \
        alternative method outlined here: https://bitbucket.org/schultzmattd/methylpy/wiki/Methylation%20Calling")
     parser_pipeline.add_argument('--keep_temp_files', type=bool, default=False, help="Boolean indicating that you'd like to keep the intermediate \
        files generated by this function. This can be useful for debugging, but in general should be left False.")
     parser_pipeline.add_argument('--bowtie2', type=bool, default=False, help="Specifies whether to use the bowtie2 aligner instead of bowtie")
     parser_pipeline.add_argument('--sort_mem', type=str, help="Parameter to pass to unix sort with -S/--buffer-size command")
     parser_pipeline.add_argument('--path_to_output', type=str, default="", help="Path to a directory where you would like the output to be stored. \
        The default is the same directory as the input fastqs.")
     parser_pipeline.add_argument('--path_to_picard', type=str, default=False, help="The path to MarkDuplicates jar from picard. Default is false indicating that you don't want to use this jar for duplication removal")
     parser_pipeline.add_argument('--remove_clonal', type=bool, default=True, help="Remove clonal reads or not")

     
     #create the parser for the "call_methylated_sites" command
     parser_call = subparsers.add_parser('call_methylated_sites', help='Use to run the call_methylated_sites function')
     parser_call.add_argument('inputf', type=str, help='inputf is the path to a bam file that contains mapped bisulfite sequencing reads')
     parser_call.add_argument('sample', type=str, help="output is the name you'd like for the allc files. The files will be named like so: allc_<sample>_<chrom>.tsv")
     parser_call.add_argument('reference', type=str, help="reference is the path to a samtools indexed fasta file")
     parser_call.add_argument('control', type=str, help="control is the name of the chromosome/region that you want to use to \
        estimate the non-conversion rate of your sample, or the non-conversion rate you'd like to use. Consequently, control \
        is either a string, or a decimal. If control is a string then it should be in the following format: 'chrom:start-end'. \
        If you'd like to specify an entire chromosome simply use 'chrom:'")
     parser_call.add_argument('casava_version', type=float, help="casava_version is a float indicating which version of casava was used to generate the fastq files.")
     parser_call.add_argument('--sig_cutoff', type=float, default=0.01, help="sig_cutoff is a float indicating the adjusted \
        p-value cutoff you wish to use for determining whether or not a site is methylated")
     parser_call.add_argument('--num_procs', type=int, default=1, help="processers is an integer indicating how many processors you'd like to run this function over") 
     parser_call.add_argument('--min_cov', type=int, default=1, help="min_cov is an integer indicating the minimum number of reads for a site to be tested")
     parser_call.add_argument('--binom_test', type=bool, default=False, help="Boolean indicating if you want to run binomial tests")
     parser_call.add_argument('--min_mc', type=int, default=0, help="Minimum number of mCs that must be observed")
     parser_call.add_argument('--path_to_samtools', type=str, default="", help='Path to samtools installation (default is current dir)')
     parser_call.add_argument('--sort_mem', type=str, default=False, help="Parameter to pass to unix sort with -S/--buffer-size command")
     parser_call.add_argument('--bh', type=bool, default=False, help="Boolean flag indicating whether or not you'd like to use the benjamini-hochberg FDR \
        instead of an FDR calculated from the control reference")
     parser_call.add_argument('--path_to_files', type=str, default="", help="string indicating the path for the output and the input bam, mpileup, or allc files \
        for methylation calling.")                                                                                                                                                   
     
          
    
def add_pe_pipeline_subparser(subparsers):
     # create the parser for the "call_mc" command
     parser_pipeline = subparsers.add_parser('run_methylation_pipeline', help='Use to run the methylation pipeline')
     parser_pipeline.add_argument('--files', type=str, nargs="+", required=True, help="list of all the fastq files you'd like to run \
        through the pipeline. Note that globbing is supported here (i.e., you can use * in your paths)")
     parser_pipeline.add_argument('--libraries', type=str, nargs="+", required=True, help="list of library IDs (in the same order as \
        the files list) indiciating which libraries each set of fastq files belong to. If you use a glob, you only need to indicate \
        the library ID for those fastqs once (i.e., the length of files and libraries should be the same)")
     parser_pipeline.add_argument('--sample', type=str, required=True, help="String indicating the name of the sample you're processing. \
        It will be included in the output files.")
     parser_pipeline.add_argument('--forward_ref', type=str, required=True, help="string indicating the path to the forward strand \
        reference created by build_ref")
     parser_pipeline.add_argument('--reverse_ref', type=str, required=True, help="string indicating the path to the reverse strand \
        reference created by build_ref")
     parser_pipeline.add_argument('--ref_fasta', type=str, required=True, help="string indicating the path to a fasta file containing \
        the sequences you used for mapping")
     parser_pipeline.add_argument('--unmethylated_control', type=str, required=True, help="name of the chromosome/region that you want \
        to use to estimate the non-conversion rate of your sample, or the non-conversion rate you'd like to use. Consequently, control \
        is either a string, or a decimal. If control is a string then it should be in the following format: 'chrom:start-end'. \
        If you'd like to specify an entire chromosome simply use 'chrom:'")
     parser_pipeline.add_argument('--path_to_samtools', type=str, default="", help='Path to samtools installation (default is current dir)')     
     parser_pipeline.add_argument('--path_to_aligner', type=str, default="", help='Path to bowtie installation (default is current dir)')
     parser_pipeline.add_argument('--aligner_options', type=str, nargs='+', help="list of strings indicating options you'd like passed to bowtie \
        (e.g., ['-k 1','-l 2']")
     parser_pipeline.add_argument('--num_procs', type=int, default=1, help='Number of processors you wish to use to \
        parallelize this function')
     parser_pipeline.add_argument('--trim_reads', type=bool, default=True, help='Whether to trim reads using cutadapt (default is True)')
     parser_pipeline.add_argument('--path_to_cutadapt', type=str, default="", help='Path to cutadapt installation (default is current dir)')
     parser_pipeline.add_argument('--adapter_seq', type=str, default="AGATCGGAAGAGCACACGTCTG", help="sequence of an adapter that was ligated \
        to the 3' end. The adapter itself and anything that follows is trimmed.")
     parser_pipeline.add_argument('--max_adapter_removal', type=int, help="Indicates the maximum number of times to try to remove adapters. \
        Useful when an adapter gets appended multiple times.")
     parser_pipeline.add_argument('--overlap_length', type=int, help="Minimum overlap length. If the overlap between the read and the adapter \
        is shorter than LENGTH, the read is not modified. This reduces the no. of bases trimmed purely due to short random adapter matches.")
     parser_pipeline.add_argument('--zero_cap', type=bool, help="Flag that causes negative quality values to be set to zero (workaround to avoid \
        segmentation faults in BWA)")
     parser_pipeline.add_argument('--error_rate', type=float, help="maximum allowed error rate (no. of errors divided by the length \
        of the matching region)")
     parser_pipeline.add_argument('--min_qual_score', type=int, default=10, help="allows you to trim low-quality ends from reads before \
        adapter removal. The algorithm is the same as the one used by BWA (Subtract CUTOFF from all qualities; compute partial sums from \
        all indices to the end of the sequence; cut sequence at the index at which the sum is minimal).")
     parser_pipeline.add_argument('--min_read_len', type=int, default=30, help="indicates the minimum length a read must be to be kept. \
        Reads that are too short even before adapter removal are also discarded. In colorspace, an initial primer is not counted.")
     parser_pipeline.add_argument('--sig_cutoff', type=float, default=.01, help="float indicating the adjusted p-value cutoff you wish to \
        use for determining whether or not a site is methylated")
     parser_pipeline.add_argument('--min_cov', type=int, default=0, help="integer indicating the minimum number of reads for a site to be tested.")
     parser_pipeline.add_argument('--binom_test', type=bool, default=False, help="Indicates that you'd like to use a binomial test, rather than the \
        alternative method outlined here: https://bitbucket.org/schultzmattd/methylpy/wiki/Methylation%20Calling")
     parser_pipeline.add_argument('--keep_temp_files', type=bool, default=False, help="Boolean indicating that you'd like to keep the intermediate \
        files generated by this function. This can be useful for debugging, but in general should be left False.")
     parser_pipeline.add_argument('--save_space', type=bool, default=True, help="indicates whether or not you'd like to perform read collapsing \
        right after mapping or once all the libraries have been mapped. If you wait until after everything has been mapped, the collapsing can be \
        parallelized. Otherwise the collapsing will have to be done serially. The trade-off is that you must keep all the mapped files around, rather \
        than deleting them as they are processed, which can take up a considerable amount of space. It's safest to set this to True.")
     parser_pipeline.add_argument('--bowtie2', type=bool, default=False, help="Specifies whether to use the bowtie2 aligner instead of bowtie")
     parser_pipeline.add_argument('--sort_mem', type=str, help="Parameter to pass to unix sort with -S/--buffer-size command")
     parser_pipeline.add_argument('--path_to_output', type=str, default="", help="Path to a directory where you would like the output to be stored. \
        The default is the same directory as the input fastqs.")
     parser_pipeline.add_argument('--path_to_picard', type=str, default=False, help="The path to picard.jar. Default is false indicating that you don't want to use this jar for duplication removal")
     parser_pipeline.add_argument('--remove_clonal', type=bool, default=True, help="Remove clonal reads or not")

     
     #create the parser for the "call_methylated_sites" command
     parser_call = subparsers.add_parser('call_methylated_sites', help='Use to run the call_methylated_sites function')
     parser_call.add_argument('inputf', type=str, help='inputf is the path to a bam file that contains mapped bisulfite sequencing reads')
     parser_call.add_argument('sample', type=str, help="output is the name you'd like for the allc files. The files will be named like so: allc_<sample>_<chrom>.tsv")
     parser_call.add_argument('reference', type=str, help="reference is the path to a samtools indexed fasta file")
     parser_call.add_argument('control', type=str, help="control is the name of the chromosome/region that you want to use to \
        estimate the non-conversion rate of your sample, or the non-conversion rate you'd like to use. Consequently, control \
        is either a string, or a decimal. If control is a string then it should be in the following format: 'chrom:start-end'. \
        If you'd like to specify an entire chromosome simply use 'chrom:'")
     parser_call.add_argument('casava_version', type=float, help="casava_version is a float indicating which version of casava was used to generate the fastq files.")
     parser_call.add_argument('--sig_cutoff', type=float, default=0.01, help="sig_cutoff is a float indicating the adjusted \
        p-value cutoff you wish to use for determining whether or not a site is methylated")
     parser_call.add_argument('--num_procs', type=int, default=1, help="processers is an integer indicating how many processors you'd like to run this function over") 
     parser_call.add_argument('--min_cov', type=int, default=1, help="min_cov is an integer indicating the minimum number of reads for a site to be tested")
     parser_call.add_argument('--binom_test', type=bool, default=False, help="Boolean indicating if you want to run binomial tests")
     parser_call.add_argument('--min_mc', type=int, default=0, help="Minimum number of mCs that must be observed")
     parser_call.add_argument('--path_to_samtools', type=str, default="", help='Path to samtools installation (default is current dir)')
     parser_call.add_argument('--sort_mem', type=str, default=False, help="Parameter to pass to unix sort with -S/--buffer-size command")
     parser_call.add_argument('--bh', type=bool, default=False, help="Boolean flag indicating whether or not you'd like to use the benjamini-hochberg FDR \
        instead of an FDR calculated from the control reference")
     parser_call.add_argument('--path_to_files', type=str, default="", help="string indicating the path for the output and the input bam, mpileup, or allc files \
        for methylation calling.")                                                                                                                                                   
     
     args = parser.parse_args()

     if args.command == "run_methylation_pipeline":
         if not args.aligner_options:
              args.aligner_options = ["-S","-k 1","-m 1","--chunkmbs 3072","--best","--strata","-o 4","-e 80","-l 20","-n 0"]
              
         run_methylation_pipeline_pe(args.files,args.libraries,args.sample,args.forward_ref,args.reverse_ref,args.ref_fasta,
                                     args.unmethylated_control,args.path_to_samtools,args.path_to_aligner,
                                     args.aligner_options,args.num_procs,args.trim_reads,
                                     args.path_to_cutadapt,args.adapter_seq,args.max_adapter_removal,
                                     args.overlap_length,args.zero_cap,args.error_rate,args.min_qual_score,args.min_read_len,
                                     args.sig_cutoff,args.min_cov,args.binom_test,args.keep_temp_files,
                                     args.save_space,
                                     args.bowtie2,args.sort_mem,args.path_to_output)
         
     elif args.command == "call_methylated_sites":
          call_methylated_sites(args.inputf, args.sample, args.reference, args.control, args.casava_version, args.sig_cutoff,
                                args.num_procs, args.min_cov, args.binom_test, args.min_mc, args.path_to_samtools, args.sort_mem,
                                args.bh, args.path_to_files)
          
if __name__ == '__main__':
     parse_args()
