import sys
import argparse
import methylpy

def parse_args():
     # create the top-level parser
     parser = argparse.ArgumentParser(
          description = "You are using methylpy "
          + methylpy.__version__
          + " version (" 
          + methylpy.__file__[:methylpy.__file__.rfind("/")]+"/"
          + ")"
     )
     subparsers = parser.add_subparsers(
          title="functions",
          dest="command",metavar="")
     
     add_build_ref_subparser(subparsers)
     add_se_pipeline_subparser(subparsers)
     add_pe_pipeline_subparser(subparsers)
     add_DMRfind_subparser(subparsers)
     add_merge_DMS_subparser(subparsers)
     add_get_methylation_level_subparser(subparsers)
     add_bam_filter_subparser(subparsers)
     add_call_mc_subparser(subparsers)
     add_allc2bw_subparser(subparsers)
     add_merge_allc_subparser(subparsers)
     add_index_allc_subparser(subparsers)
     add_filter_allc_subparser(subparsers)
     add_test_allc_subparser(subparsers)

     if len(sys.argv) > 1:
          ## print out version
          if (sys.argv[1] == '--version' or sys.argv[1] == '-v'):
               print(methylpy.__version__)
               exit()
          ## all functions
          args = parser.parse_args()
     else:
          args = parser.parse_args(["-h"])
          exit()

     if args.command == "build-reference":
          from methylpy.call_mc_se import build_ref
          build_ref(input_files=args.input_files,
                    output=args.output_prefix,
                    aligner=args.aligner,
                    path_to_aligner=args.path_to_aligner,
                    num_procs=args.num_procs,
                    buffsize=args.buffsize)

     elif args.command == "DMRfind":
          from methylpy.DMRfind import DMRfind
          DMRfind(allc_files = args.allc_files,
                  samples=args.samples,
                  mc_type=args.mc_type,
                  chroms=args.chroms,
                  num_procs=args.num_procs,
                  output_prefix=args.output_prefix,
                  min_cov=args.min_cov,
                  dmr_max_dist=args.dmr_max_dist,
                  sig_cutoff=args.sig_cutoff,
                  num_sims=args.num_sims,
                  num_sig_tests=args.min_tests,
                  min_num_dms=args.min_num_dms,
                  sample_category=args.sample_category,
                  mc_max_dist=args.mc_max_dist,
                  resid_cutoff=args.resid_cutoff,
                  keep_temp_files=args.keep_temp_files,
                  min_cluster=args.min_cluster,
                  seed=args.seed)

     elif args.command == "reidentify-DMR":
          from methylpy.DMRfind import merge_DMS_to_DMR
          merge_DMS_to_DMR(input_rms_file=args.input_rms_file,
                           output_file=args.output_file,
                           collapse_samples=args.collapse_samples,
                           sample_category=args.sample_category,
                           min_cluster=args.min_cluster,
                           sig_cutoff=args.sig_cutoff,
                           dmr_max_dist=args.dmr_max_dist,
                           min_num_dms=args.min_num_dms,
                           resid_cutoff=args.resid_cutoff,
                           num_sims=args.num_sims,
                           num_sig_tests=args.min_tests)

     elif args.command == "single-end-pipeline":
          from methylpy.call_mc_se import run_methylation_pipeline
          run_methylation_pipeline(read_files=args.read_files,
                                   sample=args.sample,
                                   forward_reference=args.forward_ref,
                                   reverse_reference=args.reverse_ref,
                                   reference_fasta=args.ref_fasta,
                                   libraries=args.libraries,                                   
                                   path_to_output=args.path_to_output,
                                   pbat=args.pbat,
                                   check_dependency=args.check_dependency,
                                   num_procs=args.num_procs,
                                   sort_mem=args.sort_mem,
                                   num_upstr_bases=args.num_upstream_bases,
                                   num_downstr_bases=args.num_downstream_bases,
                                   generate_allc_file=args.generate_allc_file,
                                   generate_mpileup_file=args.generate_mpileup_file,
                                   compress_output=args.compress_output,
                                   bgzip=args.bgzip,
                                   path_to_bgzip=args.path_to_bgzip,
                                   path_to_tabix=args.path_to_tabix,
                                   trim_reads=args.trim_reads,
                                   path_to_cutadapt=args.path_to_cutadapt,
                                   path_to_aligner=args.path_to_aligner,
                                   aligner=args.aligner,
                                   aligner_options=args.aligner_options,
                                   merge_by_max_mapq=args.merge_by_max_mapq,
                                   min_mapq=args.min_mapq,
                                   remove_clonal=args.remove_clonal,
                                   path_to_picard=args.path_to_picard,
                                   keep_clonal_stats=args.keep_clonal_stats,
                                   java_options=args.java_options,
                                   path_to_samtools=args.path_to_samtools,
                                   remove_chr_prefix=args.remove_chr_prefix,
                                   add_snp_info=args.add_snp_info,
                                   adapter_seq=args.adapter_seq,
                                   unmethylated_control=args.unmethylated_control,
                                   binom_test=args.binom_test,
                                   sig_cutoff=args.sig_cutoff,
                                   min_cov=args.min_cov,
                                   max_adapter_removal=args.max_adapter_removal,
                                   overlap_length=args.overlap_length,
                                   zero_cap=args.zero_cap,
                                   error_rate=args.error_rate,
                                   min_qual_score=args.min_qual_score,
                                   min_read_len=args.min_read_len,
                                   min_base_quality=args.min_base_quality,
                                   keep_temp_files=args.keep_temp_files)
     
     elif args.command == "paired-end-pipeline":
          from methylpy.call_mc_pe import run_methylation_pipeline_pe
          run_methylation_pipeline_pe(read1_files=args.read1_files,
                                      read2_files=args.read2_files,
                                      sample=args.sample,
                                      forward_reference=args.forward_ref,
                                      reverse_reference=args.reverse_ref,
                                      reference_fasta=args.ref_fasta,
                                      libraries=args.libraries,
                                      path_to_output=args.path_to_output,
                                      pbat=args.pbat,
                                      check_dependency=args.check_dependency,
                                      num_procs=args.num_procs,
                                      sort_mem=args.sort_mem,
                                      num_upstr_bases=args.num_upstream_bases,
                                      num_downstr_bases=args.num_downstream_bases,
                                      generate_allc_file=args.generate_allc_file,
                                      generate_mpileup_file=args.generate_mpileup_file,
                                      compress_output=args.compress_output,
                                      bgzip=args.bgzip,
                                      path_to_bgzip=args.path_to_bgzip,
                                      path_to_tabix=args.path_to_tabix,
                                      trim_reads=args.trim_reads,
                                      path_to_cutadapt=args.path_to_cutadapt,
                                      path_to_aligner=args.path_to_aligner,
                                      aligner=args.aligner,
                                      aligner_options=args.aligner_options,
                                      merge_by_max_mapq=args.merge_by_max_mapq,
                                      min_mapq=args.min_mapq,
                                      remove_clonal=args.remove_clonal,
                                      path_to_picard=args.path_to_picard,
                                      keep_clonal_stats=args.keep_clonal_stats,
                                      java_options=args.java_options,
                                      path_to_samtools=args.path_to_samtools,
                                      remove_chr_prefix=args.remove_chr_prefix,
                                      add_snp_info=args.add_snp_info,
                                      adapter_seq_read1=args.adapter_seq_read1,
                                      adapter_seq_read2=args.adapter_seq_read2,
                                      unmethylated_control=args.unmethylated_control,                                   
                                      binom_test=args.binom_test,
                                      sig_cutoff=args.sig_cutoff,
                                      min_cov=args.min_cov,
                                      max_adapter_removal=args.max_adapter_removal,
                                      overlap_length=args.overlap_length,
                                      zero_cap=args.zero_cap,
                                      error_rate=args.error_rate,
                                      min_qual_score=args.min_qual_score,
                                      min_read_len=args.min_read_len,
                                      min_base_quality=args.min_base_quality,
                                      keep_temp_files=args.keep_temp_files)

     elif args.command == "bam-quality-filter":          
          from methylpy.call_mc_se import bam_quality_mch_filter
          bam_quality_mch_filter(inputf=args.input_file,
                                 outputf=args.output_file,
                                 reference_fasta=args.ref_fasta,
                                 min_mapq=args.min_mapq,
                                 min_ch=args.min_num_ch,
                                 max_mch_level=args.max_mch_level,
                                 buffer_line_number=args.buffer_line_number,
                                 path_to_samtools=args.path_to_samtools)

     elif args.command == "call-methylation-state":
          if args.paired_end:
               from methylpy.call_mc_pe import call_methylated_sites_pe
               call_methylated_sites_pe(inputf=args.input_file,
                                        sample=args.sample,
                                        reference_fasta=args.ref_fasta,
                                        unmethylated_control=args.unmethylated_control,
                                        sig_cutoff=args.sig_cutoff,
                                        num_procs=args.num_procs,                                     
                                        num_upstr_bases=args.num_upstream_bases,
                                        num_downstr_bases=args.num_downstream_bases,
                                        generate_mpileup_file=args.generate_mpileup_file,
                                        compress_output=args.compress_output,
                                        bgzip=args.bgzip,
                                        path_to_bgzip=args.path_to_bgzip,
                                        path_to_tabix=args.path_to_tabix,
                                        min_mapq=args.min_mapq,
                                        min_cov=args.min_cov,
                                        binom_test=args.binom_test,
                                        path_to_samtools=args.path_to_samtools,
                                        remove_chr_prefix=args.remove_chr_prefix,
                                        add_snp_info=args.add_snp_info,
                                        path_to_files=args.path_to_output,
                                        min_base_quality=args.min_base_quality)
          else:
               from methylpy.call_mc_se import call_methylated_sites
               call_methylated_sites(inputf=args.input_file,
                                     sample=args.sample,
                                     reference_fasta=args.ref_fasta,
                                     unmethylated_control=args.unmethylated_control,
                                     sig_cutoff=args.sig_cutoff,
                                     num_procs=args.num_procs,                                     
                                     num_upstr_bases=args.num_upstream_bases,
                                     num_downstr_bases=args.num_downstream_bases,
                                     generate_mpileup_file=args.generate_mpileup_file,
                                     compress_output=args.compress_output,
                                     bgzip=args.bgzip,
                                     path_to_bgzip=args.path_to_bgzip,
                                     path_to_tabix=args.path_to_tabix,
                                     min_mapq=args.min_mapq,
                                     min_cov=args.min_cov,
                                     binom_test=args.binom_test,
                                     path_to_samtools=args.path_to_samtools,
                                     remove_chr_prefix=args.remove_chr_prefix,
                                     add_snp_info=args.add_snp_info,
                                     path_to_files=args.path_to_output,
                                     min_base_quality=args.min_base_quality)

     elif args.command == "add-methylation-level":
          if args.extra_info:
               from methylpy.DMRfind import get_c_info_DMRfind
               get_c_info_DMRfind(input_tsv_file=args.input_tsv_file,
                                  output=args.output_file,
                                  input_allc_files=args.allc_files,
                                  samples=args.samples,
                                  mc_type=args.mc_type,
                                  num_procs=args.num_procs,
                                  min_cov=args.min_cov,
                                  max_cov=args.max_cov,
                                  buffer_line_number=args.buffer_line_number,
                                  input_no_header=args.input_no_header)
          else:
               from methylpy.DMRfind import get_methylation_levels_DMRfind
               get_methylation_levels_DMRfind(input_tsv_file=args.input_tsv_file,
                                              output=args.output_file,
                                              input_allc_files=args.allc_files,
                                              samples=args.samples,
                                              mc_type=args.mc_type,
                                              num_procs=args.num_procs,
                                              min_cov=args.min_cov,
                                              max_cov=args.max_cov,
                                              buffer_line_number=args.buffer_line_number,
                                              input_no_header=args.input_no_header)

     elif args.command == "merge-allc":
          from methylpy.utilities import merge_allc_files
          merge_allc_files(allc_files=args.allc_files,
                           output_file=args.output_file,
                           num_procs=args.num_procs,
                           mini_batch=args.mini_batch,
                           compress_output=args.compress_output,
                           skip_snp_info=args.skip_snp_info)

     elif args.command == "index-allc":
          from methylpy.utilities import index_allc_file_batch
          index_allc_file_batch(allc_files=args.allc_files,
                                num_procs=args.num_procs,
                                reindex=args.reindex)

     elif args.command == "filter-allc":
          from methylpy.utilities import filter_allc_files
          filter_allc_files(allc_files=args.allc_files,
                            output_files=args.output_files,
                            num_procs=args.num_procs,
                            mc_type=args.mc_type,
                            chroms=args.chroms,
                            compress_output=args.compress_output,
                            max_mismatch=args.max_mismatch,
                            max_mismatch_frac=args.max_mismatch_frac,
                            min_cov=args.min_cov,
                            max_cov=args.max_cov)

     elif args.command == "test-allc":
          from methylpy.call_mc_se import perform_binomial_test
          perform_binomial_test(allc_file=args.allc_file,
                                sample=args.sample,
                                path_to_output=args.path_to_output,
                                unmethylated_control=args.unmethylated_control,
                                min_cov=args.min_cov,
                                sig_cutoff=args.sig_cutoff,
                                num_procs=args.num_procs,
                                sort_mem=args.sort_mem,
                                compress_output=args.compress_output,
                                remove_chr_prefix=args.remove_chr_prefix)

     elif args.command == "allc-to-bigwig":
          from methylpy.utilities import convert_allc_to_bigwig
          convert_allc_to_bigwig(input_allc_file=args.allc_file,
                                 output_file=args.output_file,
                                 reference_fasta=args.ref_fasta,
                                 mc_type=args.mc_type,
                                 bin_size=args.bin_size,
                                 path_to_wigToBigWig=args.path_to_wigToBigWig,
                                 path_to_samtools=args.path_to_samtools,
                                 min_bin_sites=args.min_bin_sites,
                                 min_bin_cov=args.min_bin_cov,
                                 min_site_cov=args.min_site_cov,
                                 max_site_cov=args.max_site_cov,
                                 remove_chr_prefix=args.remove_chr_prefix,
                                 add_chr_prefix=args.add_chr_prefix)

def add_DMRfind_subparser(subparsers):
     # create the parser for the "DMRfind" command
     parser_dmrfind = subparsers.add_parser(
          "DMRfind",
          formatter_class=argparse.ArgumentDefaultsHelpFormatter,
          help="Identify differentially methylated regions")
     
     # add options
     parser_dmrfind_req = parser_dmrfind.add_argument_group("required inputs")
     parser_dmrfind_req.add_argument("--allc-files",
                                     type=str,
                                     nargs="+",
                                     required=True,
                                     help="List of allc files.")

     parser_dmrfind_req.add_argument("--output-prefix",
                                     type=str,
                                     required=True,
                                     help="String indicating the prefix for output files")
     
     parser_dmrfind_opt = parser_dmrfind.add_argument_group("optional inputs")

     parser_dmrfind_opt.add_argument("--samples",
                                     type=str,
                                     nargs="+",
                                     default=None,
                                     help="List of space separated samples matching allc files. By default "
                                     +"sample names will be inferred from allc filenames")
     
     parser_dmrfind_opt.add_argument("--chroms",
                                     type=str,
                                     nargs="+",
                                     required=False,
                                     default=None,
                                     help="Space separated listing of chromosomes where DMRs will "
                                     +"be called. If not specified, DMRs will be called across the chromosomes/contigs "
                                     +"that contained any data in all allc files.")

     parser_dmrfind_opt.add_argument("--mc-type",
                                     type=str,
                                     nargs="+",
                                     default=["CGN"],
                                     help="List of space separated mc nucleotide contexts for "
                                     + "which you want to look for DMRs. These classifications "
                                     + "may use the wildcards H (indicating anything but a G) and "
                                     + "N (indicating any nucleotide).")
     
     parser_dmrfind_opt.add_argument("--num-procs",
                                     type=int,
                                     default=1,
                                     help="Number of processors you wish to use to parallelize this function")
     
     parser_dmrfind_opt.add_argument("--min-cov",
                                     type=int,
                                     default=0,
                                     help="Minimum number of reads that must cover a site for it to be "
                                     + "considered.")

     parser_dmrfind_opt.add_argument("--dmr-max-dist",
                                     type=int,
                                     default=250,
                                     help="Maximum distance two significant sites can be to be included "
                                     + "in the same block.")
     
     parser_dmrfind_opt.add_argument("--sig-cutoff",
                                     type=float,
                                     default=.01,
                                     help="Float indicating at what FDR you want to consider a result "
                                     + "significant.")

     parser_dmrfind_opt.add_argument("--num-sims",
                                     type=int,
                                     default=3000,
                                     help="Number of permutation tests you would like to run to estimate "
                                     + "the p-values of the differential methylation tests")
     
     parser_dmrfind_opt.add_argument("--min-tests",
                                     type=int,
                                     default=100,
                                     help="Minimum number of permuation tests you\ would d like to run "
                                     + "for each mC")     
     
     parser_dmrfind_opt.add_argument("--min-num-dms",
                                     type=int,
                                     default=0,
                                     help="The minimum number of differentially methylated sites "
                                     + "that a differentially methylated region needs to contain to be "
                                     + "reported")
     
     parser_dmrfind_opt.add_argument("--sample-category",
                                     type=str,
                                     nargs="+",
                                     default=False,
                                     help="A list of categories that each respective sample belongs "
                                     + "to; the categories must begin at 0 and increase by 1 for "
                                     + "each category added. ex: samples [A,B,C] categories [0,1,2] "
                                     + "or categories [0, 1, 0] ")
     
     parser_dmrfind_opt.add_argument("--mc-max-dist",
                                     type=int,
                                     default=0,
                                     help="Integer indicating the maximum distance two sites can be "
                                     + "from one another for their methylation counts to be combined. "
                                     + "This option helps with low coverage experiments where you may "
                                     + "want to leverage the correlation of methylation between sites "
                                     + "to get more statistical power.")
     
     parser_dmrfind_opt.add_argument("--resid-cutoff",
                                     type=int,
                                     default=0.01,
                                     help="Results will have to show deviations in the contingency "
                                     + "table in the same direction as the rest of the window")
     
     parser_dmrfind_opt.add_argument("--keep-temp-files",
                                     type=str2bool,
                                     default=False,
                                     help="Boolean indicating that you would like to keep the intermediate files "
                                     +"generated by this function. This can be useful for debugging, but in general "
                                     +"should be left False.")

     
     parser_dmrfind_opt.add_argument("--min-cluster",
                                     type=int,
                                     default=2,
                                     help="The minimum number of each sample category that must be "
                                     + "present in every block that is output.")
     
     parser_dmrfind_opt.add_argument("--seed",
                                     type=int,
                                     default=-1,
                                     help="A seed to provide to the random number generator for "
                                     + "permutation testing. Only change this if you are debugging "
                                     + "and want to make sure the permutation output is consistent")

def add_merge_DMS_subparser(subparsers):
     # create the parser for the "merge_DMS" command
     parser_mergedms = subparsers.add_parser(
          "reidentify-DMR",
          formatter_class=argparse.ArgumentDefaultsHelpFormatter,
          help="Re-call DMRs from existing DMRfind result")
     
     # add options
     parser_mergedms_req = parser_mergedms.add_argument_group("required inputs")
     parser_mergedms_req.add_argument("--input-rms-file",
                                      type=str,
                                      required=True,
                                      help="File storing the results of RMS tests (from DMRfind function.")

     parser_mergedms_req.add_argument("--output-file",
                                      type=str,
                                      required=True,
                                      help="String indicating the name of output file")
     
     parser_mergedms_opt = parser_mergedms.add_argument_group("optional inputs")

     parser_mergedms_opt.add_argument("--collapse-samples",
                                      type=str,
                                      nargs="+",
                                      default=False,
                                      help="A list of samples for collapsing blocks")     
     
     parser_mergedms_opt.add_argument("--sample-category",
                                      type=str,
                                      nargs="+",
                                      default=False,
                                      help="A list of categories that each respective sample belongs "
                                      + "to; the categories must begin at 0 and increase by 1 for "
                                      + "each category added. ex: samples [A,B,C] categories [0,1,2] "
                                      + "or categories [0, 1, 0] ")
     
     parser_mergedms_opt.add_argument("--min-cluster",
                                      type=int,
                                      default=2,
                                      help="The minimum number of each sample category that must be "
                                      + "present in every block that is output.")
     
     parser_mergedms_opt.add_argument("--sig-cutoff",
                                      type=float,
                                      default=.01,
                                      help="Float indicating at what FDR you want to consider a result "
                                      + "significant.")

     parser_mergedms_opt.add_argument("--dmr-max-dist",
                                      type=int,
                                      default=250,
                                      help="Maximum distance two significant sites can be to be included "
                                      + "in the same block.")
     
     parser_mergedms_opt.add_argument("--min-num-dms",
                                      type=int,
                                      default=0,
                                      help="The minimum number of differentially methylated sites "
                                      + "that a differentially methylated region needs to contain to be "
                                      + "reported")

     parser_mergedms_opt.add_argument("--resid-cutoff",
                                      type=int,
                                      default=0.01,
                                      help="Results will have to show deviations in the contingency "
                                      + "table in the same direction as the rest of the window")

     parser_mergedms_opt.add_argument("--num-sims",
                                      type=int,
                                      default=3000,
                                      help="Number of permutation tests you would like to run to estimate "
                                      + "the p-values of the differential methylation tests")
     
     parser_mergedms_opt.add_argument("--min-tests",
                                      type=int,
                                      default=100,
                                      help="Minimum number of permuation tests you\ would d like to run "
                                      + "for each mC")     

def add_se_pipeline_subparser(subparsers):
     parser_se = subparsers.add_parser("single-end-pipeline",
                                       formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                       help="Methylation pipeline for single-end data")

     parser_se_req = parser_se.add_argument_group("required inputs")
          
     parser_se_req.add_argument("--read-files",
                                  type=str,
                                  nargs="+",
                                  required=True,
                                  help="list of all the fastq files you would like to run through "
                                  + "the pipeline. Note that globbing is supported here (i.e., you "
                                  + "can use * in your paths)")
     
     parser_se_req.add_argument("--sample",
                                  type=str,
                                  required=True,
                                  help="String indicating the name of the sample you are processing. "
                                  + "It will be included in the output files.")
     
     parser_se_req.add_argument("--forward-ref",
                                  type=str, required=True, help="string indicating the path to the "
                                  + "forward strand reference created by build_ref")
     
     parser_se_req.add_argument("--reverse-ref",
                                  type=str,
                                  required=True,
                                  help="string indicating the path to the reverse strand reference "
                                  + "created by build_ref")
     
     parser_se_req.add_argument("--ref-fasta",
                                  type=str,
                                  required=True,
                                  help="string indicating the path to a fasta file containing the "
                                  + "sequences you used for mapping")
     
     parser_se_opt = parser_se.add_argument_group("optional inputs")

     parser_se_opt.add_argument("--libraries",
                                type=str,
                                nargs="+",
                                default=["libA"],
                                help="list of library IDs (in the same order as the files list) "
                                + "indiciating which libraries each set of fastq files belong to. "
                                + "If you use a glob, you only need to indicate the library ID for "
                                + "those fastqs once (i.e., the length of files and libraries should "
                                + "be the same)")
     
     parser_se_opt.add_argument("--path-to-output",
                                type=str,
                                default="",
                                help="Path to a directory where you would like the output to be stored. "
                                + "The default is the same directory as the input fastqs.")

     parser_se_opt.add_argument("--pbat",
                                type=str2bool,
                                default=False,
                                help="Boolean indicating whether to process data in PBAT (Post-Bisulfite "
                                +"Adaptor Tagging) mode, in which reads will be mapped to opposite strand "
                                +"of C-T converted genome and the forward strand of G-A converted genome.")

     parser_se_opt.add_argument("--check-dependency",
                                type=str2bool,
                                default=False,
                                help="Boolean indicating whether to check dependency requirements are met.")

     parser_se_opt.add_argument("--num-procs",
                                type=int,
                                default=1,
                                help="Number of processors you wish to use to parallelize this function")     

     parser_se_opt.add_argument("--sort-mem",
                                type=str,
                                default="500M",
                                help="Parameter to pass to unix sort with -S/--buffer-size command")

     parser_se_opt.add_argument("--num-upstream-bases",
                                type=int,
                                default=0,
                                help="Number of base(s) upstream of each cytosine that you wish to include "
                                + "in output file. Recommend value 1 for NOMe-seq processing since the "
                                + "upstream base is required to tell apart cytosine at GC context.")

     parser_se_opt.add_argument("--num-downstream-bases",
                                type=int,
                                default=2,
                                help="Number of base(s) downstream of each cytosine that you wish to include "
                                + "in output file. Recommend value to be at least 1 to separate cytosines at "
                                + "different sequence context.")

     parser_se_opt.add_argument("--generate-allc-file",
                                type=str2bool,
                                default=True,
                                help="Boolean indicating whether to generate the final output file that "
                                +" contains the methylation state of each cytosine. If set to be false, "
                                +"only alignment file (in BAM format) will be generated.")

     parser_se_opt.add_argument("--generate-mpileup-file",
                                type=str2bool,
                                default=True,
                                help="Boolean indicating whether to generate intermediate mpileup file to save "
                                +"space. However, skipping mpileup step may cause problem due to the nature of "
                                +"python. Not skipping this step is recommended.")
     
     parser_se_opt.add_argument("--compress-output",
                                type=str2bool,
                                default=True,
                                help="Boolean indicating whether to compress (by gzip) the final output "
                                + "(allc file(s)).")     

     parser_se_opt.add_argument("--bgzip",
                                type=str2bool,
                                default=False,
                                help="Boolean indicating whether to bgzip compressed allc files and tabix index.")
     
     parser_se_opt.add_argument("--path-to-bgzip",
                                type=str,
                                default="",
                                help="Path to bgzip installation")

     parser_se_opt.add_argument("--path-to-tabix",
                                type=str,
                                default="",
                                help="Path to tabix installation")

     parser_se_opt.add_argument("--trim-reads",
                                type=str2bool,
                                default=True,
                                help="Boolean indicating whether to trim reads using cutadapt.")
     
     parser_se_opt.add_argument("--path-to-cutadapt",
                                type=str,
                                default="",
                                help="Path to cutadapt installation")
     
     parser_se_opt.add_argument("--path-to-aligner",
                                type=str,
                                default="",
                                help="Path to bowtie/bowtie2 installation")
     
     parser_se_opt.add_argument("--aligner",
                                type=str,
                                default="bowtie2",
                                help="Aligner to use. Currently, methylpy supports bowtie, bowtie2 and minimap2. ")

     parser_se_opt.add_argument("--aligner-options",
                                type=str,
                                nargs="+",
                                help="list of strings indicating options you would like passed "
                                + "to bowtie (e.g., \"-k 1 -l 2\")")

     parser_se_opt.add_argument("--merge-by-max-mapq",
                                type=str2bool,
                                default=False,
                                help="Boolean indicates whether to merge alignment results from two "
                                +"converted genomes by MAPQ score. Be default, we only keep reads that "
                                +"are mapped to only one of the two converted genomes. If this option "
                                +"is set to True, for a read that could be mapped to both converted "
                                +"genomes, the alignment that achieves larger MAPQ score will be kept.")

     parser_se_opt.add_argument("--remove-clonal",
                                type=str2bool,
                                default=False,
                                help="Boolean indicates whether to remove clonal reads or not")
     
     parser_se_opt.add_argument("--path-to-picard",
                                type=str,
                                default="",
                                help="The path to the picard.jar in picard tools. The jar file can "
                                + "be downloaded from https://broadinstitute.github.io/picard/index.html "
                                + "(default is current dir)")
          
     parser_se_opt.add_argument("--keep-clonal-stats",
                                type=str2bool,
                                default=True,
                                help="Boolean indicates whether to store the metric file from picard.")
     
     parser_se_opt.add_argument("--java-options",
                                type=str,
                                default="-Xmx20g",
                                help="String indicating the option pass the java when running picard.")

     parser_se_opt.add_argument("--path-to-samtools",
                                type=str,
                                default="",
                                help="Path to samtools installation")
     
     parser_se_opt.add_argument("--adapter-seq",
                                type=str,
                                default="AGATCGGAAGAGCACACGTCTG",
                                help="sequence of an adapter that was ligated to the 3\' end. The "
                                +"adapter itself and anything that follows is trimmed.")

     parser_se_opt.add_argument("--remove-chr-prefix",
                                type=str2bool,
                                default=True,
                                help="Boolean indicates whether to remove in the final output the \"chr\" prefix "
                                +"in the chromosome name")

     parser_se_opt.add_argument("--add-snp-info",
                                type=str2bool,
                                default=False,
                                help="Boolean indicates whether to add extra two columns in the output (allc) file "
                                +"regarding the genotype information of each site. The first (second) column contain "
                                +"the number of basecalls that support the reference gentype (variant) for nucleotides "
                                "in the sequence context.")
     
     parser_se_opt.add_argument("--unmethylated-control",
                                type=str,
                                default=None,
                                help="name of the chromosome/region that you want to use to estimate "
                                + "the non-conversion rate of your sample, or the non-conversion rate "
                                + "you would like to use. Consequently, control is either a string, or "
                                + "a decimal. If control is a string then it should be in the following "
                                + "format: \"chrom:start-end\". If you would like to specify an entire "
                                + "chromosome simply use \"chrom:\"")
     
     parser_se_opt.add_argument("--binom-test",
                                type=str2bool,
                                default=False,
                                help="Indicates that you would like to perform a binomial test on each cytosine "
                                +"to delineate cytosines that are significantly methylated than noise due to "
                                +"the failure of bisulfite conversion.")
     
     parser_se_opt.add_argument("--sig-cutoff",
                                type=float,
                                default=.01,
                                help="float indicating the adjusted p-value cutoff you wish to use for "
                                + "determining whether or not a site is methylated")

     parser_se_opt.add_argument("--min-mapq",
                                type=int,
                                default=30,
                                help="Minimum MAPQ for reads to be included.")

     parser_se_opt.add_argument("--min-cov",
                                type=int,
                                default=0,
                                help="Integer indicating the minimum number of reads for a site to be tested.")
     
     parser_se_opt.add_argument("--max-adapter-removal",
                                type=int,
                                help="Indicates the maximum number of times to try to remove adapters. Useful "
                                +"when an adapter gets appended multiple times.")
     
     parser_se_opt.add_argument("--overlap-length",
                                type=int,
                                help="Minimum overlap length. If the overlap between the read and the adapter "
                                +"is shorter than LENGTH, the read is not modified. This reduces the no. of "
                                +"bases trimmed purely due to short random adapter matches.")
     
     parser_se_opt.add_argument("--zero-cap",
                                type=str2bool,
                                help="Flag that causes negative quality values to be set to zero (workaround "
                                +"to avoid segmentation faults in BWA)")
     
     parser_se_opt.add_argument("--error-rate",
                                type=float,
                                help="maximum allowed error rate (no. of errors divided by the length of "
                                +"the matching region)")
     
     parser_se_opt.add_argument("--min-qual-score",
                                type=int,
                                default=10,
                                help="allows you to trim low-quality ends from reads before adapter removal. "
                                +"The algorithm is the same as the one used by BWA (Subtract CUTOFF from all "
                                +"qualities; compute partial sums from all indices to the end of the sequence; "
                                +" cut sequence at the index at which the sum is minimal).")
     
     parser_se_opt.add_argument("--min-read-len",
                                type=int,
                                default=30,
                                help="indicates the minimum length a read must be to be kept. Reads that "
                                +"are too short even before adapter removal are also discarded. In colorspace, "
                                "an initial primer is not counted.")

     parser_se_opt.add_argument("--min-base-quality",
                                type=int,
                                default=1,
                                help="Integer indicating the minimum PHRED quality score for a base to be "
                                +"included in the mpileup file (and subsequently to be considered for "
                                +"methylation calling).")
     
     parser_se_opt.add_argument("--keep-temp-files",
                                type=str2bool,
                                default=False,
                                help="Boolean indicating that you would like to keep the intermediate files "
                                +"generated by this function. This can be useful for debugging, but in general "
                                +"should be left False.")

def add_pe_pipeline_subparser(subparsers):
     parser_pe = subparsers.add_parser("paired-end-pipeline",
                                       formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                       help="Methylation pipeline for paired-end data")

     parser_pe_req = parser_pe.add_argument_group("required inputs")
          
     parser_pe_req.add_argument("--read1-files",
                                  type=str,
                                  nargs="+",
                                  required=True,
                                  help="list of all the read 1 fastq files you would like to run through "
                                  + "the pipeline. Note that globbing is supported here (i.e., you "
                                  + "can use * in your paths)")
     
     parser_pe_req.add_argument("--read2-files",
                                  type=str,
                                  nargs="+",
                                  required=True,
                                  help="list of all the read 2 fastq files you would like to run through "
                                  + "the pipeline. Note that globbing is supported here (i.e., you "
                                  + "can use * in your paths)")
     
     parser_pe_req.add_argument("--sample",
                                  type=str,
                                  required=True,
                                  help="String indicating the name of the sample you are processing. "
                                  + "It will be included in the output files.")
     
     parser_pe_req.add_argument("--forward-ref",
                                  type=str, required=True, help="string indicating the path to the "
                                  + "forward strand reference created by build_ref")
     
     parser_pe_req.add_argument("--reverse-ref",
                                  type=str,
                                  required=True,
                                  help="string indicating the path to the reverse strand reference "
                                  + "created by build_ref")
     
     parser_pe_req.add_argument("--ref-fasta",
                                  type=str,
                                  required=True,
                                  help="string indicating the path to a fasta file containing the "
                                  + "sequences you used for mapping")
     
     parser_pe_opt = parser_pe.add_argument_group("optional inputs")

     parser_pe_opt.add_argument("--libraries",
                                type=str,
                                nargs="+",
                                default=["libA"],
                                help="list of library IDs (in the same order as the files list) "
                                + "indiciating which libraries each set of fastq files belong to. "
                                + "If you use a glob, you only need to indicate the library ID for "
                                + "those fastqs once (i.e., the length of files and libraries should "
                                + "be the same)")
     
     parser_pe_opt.add_argument("--path-to-output",
                                type=str,
                                default="",
                                help="Path to a directory where you would like the output to be stored. "
                                + "The default is the same directory as the input fastqs.")

     parser_pe_opt.add_argument("--pbat",
                                type=str2bool,
                                default=False,
                                help="Boolean indicating whether to process data in PBAT (Post-Bisulfite "
                                +"Adaptor Tagging) mode, in which reads will be mapped to opposite strand "
                                +"of C-T converted genome and the forward strand of G-A converted genome.")

     parser_pe_opt.add_argument("--check-dependency",
                                type=str2bool,
                                default=False,
                                help="Boolean indicating whether to check dependency requirements are met.")

     parser_pe_opt.add_argument("--num-procs",
                                type=int,
                                default=1,
                                help="Number of processors you wish to use to parallelize this function")     

     parser_pe_opt.add_argument("--sort-mem",
                                type=str,
                                default="500M",
                                help="Parameter to pass to unix sort with -S/--buffer-size command")

     parser_pe_opt.add_argument("--num-upstream-bases",
                                type=int,
                                default=0,
                                help="Number of base(s) upstream of each cytosine that you wish to include "
                                + "in output file. Recommend value 1 for NOMe-seq processing since the "
                                + "upstream base is required to tell apart cytosine at GC context.")

     parser_pe_opt.add_argument("--num-downstream-bases",
                                type=int,
                                default=2,
                                help="Number of base(s) downstream of each cytosine that you wish to include "
                                + "in output file. Recommend value to be at least 1 to separate cytosines at "
                                + "different sequence contexts.")

     parser_pe_opt.add_argument("--generate-allc-file",
                                type=str2bool,
                                default=True,
                                help="Boolean indicating whether to generate the final output file that "
                                +" contains the methylation state of each cytosine. If set to be false, "
                                +"only alignment file (in BAM format) will be generated.")

     parser_pe_opt.add_argument("--generate-mpileup-file",
                                type=str2bool,
                                default=True,
                                help="Boolean indicating whether to generate intermediate mpileup file to save "
                                +"space. However, skipping mpileup step may cause problem due to the nature of "
                                +"python. Not skipping this step is recommended.")
     
     parser_pe_opt.add_argument("--compress-output",
                                type=str2bool,
                                default=True,
                                help="Boolean indicating whether to compress (by gzip) the final output "
                                + "(allc file(s)).")     

     parser_pe_opt.add_argument("--bgzip",
                                type=str2bool,
                                default=False,
                                help="Boolean indicating whether to bgzip compressed allc files and tabix index.")
     
     parser_pe_opt.add_argument("--path-to-bgzip",
                                type=str,
                                default="",
                                help="Path to bgzip installation")

     parser_pe_opt.add_argument("--path-to-tabix",
                                type=str,
                                default="",
                                help="Path to tabix installation")

     parser_pe_opt.add_argument("--trim-reads",
                                type=str2bool,
                                default=True,
                                help="Boolean indicating whether to trim reads using cutadapt.")
     
     parser_pe_opt.add_argument("--path-to-cutadapt",
                                type=str,
                                default="",
                                help="Path to cutadapt installation")
     
     parser_pe_opt.add_argument("--path-to-aligner",
                                type=str,
                                default="",
                                help="Path to bowtie/bowtie2 installation")
     
     parser_pe_opt.add_argument("--aligner",
                                type=str,
                                default="bowtie2",
                                help="Aligner to use. Currently, methylpy supports bowtie, bowtie2 and minimap2. ")
     
     parser_pe_opt.add_argument("--aligner-options",
                                type=str,
                                nargs="+",
                                help="list of strings indicating options you would like passed "
                                + "to bowtie (e.g., \"-k 1 -l 2\")")

     parser_pe_opt.add_argument("--merge-by-max-mapq",
                                type=str2bool,
                                default=False,
                                help="Boolean indicates whether to merge alignment results from two "
                                +"converted genomes by MAPQ score. Be default, we only keep read pairs "
                                "that are mapped to only one of the two converted genomes. If this option "
                                +"is set to True, for a read pair that could be mapped to both converted "
                                +"genomes, the alignment that achieves larger MAPQ score will be kept.")

     parser_pe_opt.add_argument("--remove-clonal",
                                type=str2bool,
                                default=False,
                                help="Boolean indicates whether to remove clonal reads or not")
     
     parser_pe_opt.add_argument("--path-to-picard",
                                type=str,
                                default="",
                                help="The path to the picard.jar in picard tools. The jar file can "
                                + "be downloaded from https://broadinstitute.github.io/picard/index.html "
                                + "(default is current dir)")
          
     parser_pe_opt.add_argument("--keep-clonal-stats",
                                type=str2bool,
                                default=True,
                                help="Boolean indicates whether to store the metric file from picard.")
     
     parser_pe_opt.add_argument("--java-options",
                                type=str,
                                default="-Xmx20g",
                                help="String indicating the option pass the java when running picard.")

     parser_pe_opt.add_argument("--path-to-samtools",
                                type=str,
                                default="",
                                help="Path to samtools installation")
     
     parser_pe_opt.add_argument("--adapter-seq-read1",
                                type=str,
                                default="AGATCGGAAGAGCACACGTCTGAAC",
                                help="sequence of an adapter that was ligated to the 3\' end of read 1. The "
                                +"adapter itself and anything that follows is trimmed.")
     
     parser_pe_opt.add_argument("--adapter-seq-read2",
                                type=str,
                                default="AGATCGGAAGAGCGTCGTGTAGGGA",
                                help="sequence of an adapter that was ligated to the 3\' end of read 2. The "
                                +"adapter itself and anything that follows is trimmed.")

     parser_pe_opt.add_argument("--remove-chr-prefix",
                                type=str2bool,
                                default=True,
                                help="Boolean indicates whether to remove in the final output the \"chr\" prefix "
                                +"in the chromosome name")

     parser_pe_opt.add_argument("--add-snp-info",
                                type=str2bool,
                                default=False,
                                help="Boolean indicates whether to add extra two columns in the output (allc) file "
                                +"regarding the genotype information of each site. The first (second) column contain "
                                +"the number of basecalls that support the reference gentype (variant) for nucleotides "
                                "in the sequence context.")

     parser_pe_opt.add_argument("--unmethylated-control",
                                type=str,
                                default=None,
                                help="name of the chromosome/region that you want to use to estimate "
                                + "the non-conversion rate of your sample, or the non-conversion rate "
                                + "you would like to use. Consequently, control is either a string, or "
                                + "a decimal. If control is a string then it should be in the following "
                                + "format: \"chrom:start-end\". If you would like to specify an entire "
                                + "chromosome simply use \"chrom:\"")
     
     parser_pe_opt.add_argument("--binom-test",
                                type=str2bool,
                                default=False,
                                help="Indicates that you would like to perform a binomial test on each cytosine "
                                +"to delineate cytosines that are significantly methylated than noise due to "
                                +"the failure of bisulfite conversion.")
     
     parser_pe_opt.add_argument("--sig-cutoff",
                                type=float,
                                default=.01,
                                help="float indicating the adjusted p-value cutoff you wish to use for "
                                + "determining whether or not a site is methylated")

     parser_pe_opt.add_argument("--min-mapq",
                                type=int,
                                default=30,
                                help="Minimum MAPQ for reads to be included.")

     parser_pe_opt.add_argument("--min-cov",
                                type=int,
                                default=0,
                                help="Integer indicating the minimum number of reads for a site to be tested.")
     
     parser_pe_opt.add_argument("--max-adapter-removal",
                                type=int,
                                help="Indicates the maximum number of times to try to remove adapters. Useful "
                                +"when an adapter gets appended multiple times.")
     
     parser_pe_opt.add_argument("--overlap-length",
                                type=int,
                                help="Minimum overlap length. If the overlap between the read and the adapter "
                                +"is shorter than LENGTH, the read is not modified. This reduces the no. of "
                                +"bases trimmed purely due to short random adapter matches.")
     
     parser_pe_opt.add_argument("--zero-cap",
                                type=str2bool,
                                help="Flag that causes negative quality values to be set to zero (workaround "
                                +"to avoid segmentation faults in BWA)")
     
     parser_pe_opt.add_argument("--error-rate",
                                type=float,
                                help="maximum allowed error rate (no. of errors divided by the length of "
                                +"the matching region)")
     
     parser_pe_opt.add_argument("--min-qual-score",
                                type=int,
                                default=10,
                                help="allows you to trim low-quality ends from reads before adapter removal. "
                                +"The algorithm is the same as the one used by BWA (Subtract CUTOFF from all "
                                +"qualities; compute partial sums from all indices to the end of the sequence; "
                                +" cut sequence at the index at which the sum is minimal).")
     
     parser_pe_opt.add_argument("--min-read-len",
                                type=int,
                                default=30,
                                help="indicates the minimum length a read must be to be kept. Reads that "
                                +"are too short even before adapter removal are also discarded. In colorspace, "
                                "an initial primer is not counted.")

     parser_pe_opt.add_argument("--min-base-quality",
                                type=int,
                                default=1,
                                help="Integer indicating the minimum PHRED quality score for a base to be "
                                +"included in the mpileup file (and subsequently to be considered for "
                                +"methylation calling).")
     
     parser_pe_opt.add_argument("--keep-temp-files",
                                type=str2bool,
                                default=False,
                                help="Boolean indicating that you would like to keep the intermediate files "
                                +"generated by this function. This can be useful for debugging, but in general "
                                +"should be left False.")

def add_build_ref_subparser(subparsers):
     # create the parser for the "DMRfind" command
     parser_build = subparsers.add_parser(
          "build-reference",
          formatter_class=argparse.ArgumentDefaultsHelpFormatter,
          help="Building reference for bisulfite sequencing data")
     
     # add options
     parser_build_req = parser_build.add_argument_group("required inputs")
     parser_build_req.add_argument("--input-files",
                                   type=str,
                                   nargs="+",
                                   required=True,
                                   help="List of genome fasta files to build a reference from.")

     parser_build_req.add_argument("--output-prefix",
                                   type=str,
                                   required=True,
                                   help="the prefix of the two output reference files that will be created.")
     
     parser_build_opt = parser_build.add_argument_group("optional inputs")
     parser_build_opt.add_argument("--num-procs",
                                   type=int,
                                   default=1,
                                   help="Number of processors you wish to use to parallelize this function")

     parser_build_opt.add_argument("--aligner",
                                   type=str,
                                   default="bowtie2",
                                   help="Aligner to use. Currently, methylpy supports bowtie, bowtie2 and minimap2. ")

     parser_build_opt.add_argument("--path-to-aligner",
                                   type=str,
                                   default="",
                                   help="Path to bowtie/bowtie2 installation")

     parser_build_opt.add_argument("--buffsize",
                                   type=int,
                                   default=100,
                                   help="The number of bytes that will be read in from the reference at once.")


def add_bam_filter_subparser(subparsers):
     parser_filter = subparsers.add_parser(
          "bam-quality-filter",
          formatter_class=argparse.ArgumentDefaultsHelpFormatter,
          help="Filter out single-end reads by mapping quality and mCH level")
     
     # add options
     parser_filter_req = parser_filter.add_argument_group("required inputs")
     parser_filter_req.add_argument("--input-file",
                                    type=str,
                                    required=True,
                                    help="BAM file to filter.")
     
     parser_filter_req.add_argument("--output-file",
                                    type=str,
                                    required=True,
                                    help="Name of output file")

     parser_filter_req.add_argument("--ref-fasta",
                                    type=str,
                                    required=True,
                                    help="string indicating the path to a fasta file containing the "
                                    +"sequences you used for mapping")
     
     parser_filter_opt = parser_filter.add_argument_group("optional inputs")
     parser_filter_opt.add_argument("--min-mapq",
                                    type=int,
                                    default=30,
                                    help="Minimum MAPQ for reads to be included.")

     parser_filter_opt.add_argument("--min-num-ch",
                                    type=int,
                                    default=30,
                                    help="Minimum number of CH sites for mCH level filter to be applied.")

     parser_filter_opt.add_argument("--max-mch-level",
                                    type=float,
                                    default=0.7,
                                    help="Maximum mCH level for reads to be included.")

     parser_filter_opt.add_argument("--path-to-samtools",
                                type=str,
                                default="",
                                help="Path to samtools installation")
     
     parser_filter_opt.add_argument("--buffer-line-number",
                                type=int,
                                default=100000,
                                help="size of buffer for reads to be written on hard drive.")

def add_call_mc_subparser(subparsers):
     call_mc = subparsers.add_parser("call-methylation-state",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     help="Call cytosine methylation state from BAM file")

     call_mc_req = call_mc.add_argument_group("required inputs")     
     call_mc_req.add_argument("--input-file",
                              type=str,
                              help="bam file that contains mapped bisulfite sequencing reads.")
     
     call_mc_req.add_argument("--sample",
                              type=str,
                              required=True,
                              help="String indicating the name of the sample you are processing. "
                              + "It will be included in the output files.")

     call_mc_req.add_argument("--ref-fasta",
                              type=str,
                              required=True,
                              help="string indicating the path to a fasta file containing the "
                              + "sequences you used for mapping")

     call_mc_req.add_argument("--paired-end",
                              type=str2bool,
                              required=True,
                              default=False,
                              help="Boolean indicating whether the input BAM file is from paired-end "
                              +"data.")
     
     call_mc_opt = call_mc.add_argument_group("optional inputs")

     call_mc_opt.add_argument("--path-to-output",
                              type=str,
                              default="",
                              help="Path to a directory where you would like the output to be stored. "
                              + "The default is the same directory as the input fastqs.")

     call_mc_opt.add_argument("--num-procs",
                              type=int,
                              default=1,
                              help="Number of processors you wish to use to parallelize this function")     

     call_mc_opt.add_argument("--num-upstream-bases",
                              type=int,
                              default=0,
                              help="Number of base(s) upstream of each cytosine that you wish to include "
                              + "in output file. Recommend value 1 for NOMe-seq processing since the "
                              + "upstream base is required to tell apart cytosine at GC context.")

     call_mc_opt.add_argument("--num-downstream-bases",
                              type=int,
                              default=2,
                              help="Number of base(s) downstream of each cytosine that you wish to include "
                              + "in output file. Recommend value to be at least 1 to separate cytosines at "
                              + "different sequence contexts.")

     call_mc_opt.add_argument("--generate-allc-file",
                              type=str2bool,
                              default=True,
                              help="Boolean indicating whether to generate the final output file that "
                              +" contains the methylation state of each cytosine. If set to be false, "
                              +"only alignment file (in BAM format) will be generated.")

     call_mc_opt.add_argument("--generate-mpileup-file",
                              type=str2bool,
                              default=True,
                              help="Boolean indicating whether to generate intermediate mpileup file to save "
                              +"space. However, skipping mpileup step may cause problem due to the nature of "
                              +"python. Not skipping this step is recommended.")
     
     call_mc_opt.add_argument("--compress-output",
                              type=str2bool,
                              default=True,
                              help="Boolean indicating whether to compress (by gzip) the final output "
                              + "(allc file(s)).")
     
     call_mc_opt.add_argument("--bgzip",
                              type=str2bool,
                              default=False,
                              help="Boolean indicating whether to bgzip compressed allc files and tabix index.")
     
     call_mc_opt.add_argument("--path-to-bgzip",
                              type=str,
                              default="",
                              help="Path to bgzip installation")

     call_mc_opt.add_argument("--path-to-tabix",
                              type=str,
                              default="",
                              help="Path to tabix installation")

     call_mc_opt.add_argument("--path-to-samtools",
                              type=str,
                              default="",
                              help="Path to samtools installation")

     call_mc_opt.add_argument("--remove-chr-prefix",
                              type=str2bool,
                              default=True,
                              help="Boolean indicates whether to remove in the final output the \"chr\" prefix "
                              +"in the chromosome name")

     call_mc_opt.add_argument("--add-snp-info",
                              type=str2bool,
                              default=False,
                              help="Boolean indicates whether to add extra two columns in the output (allc) file "
                              +"regarding the genotype information of each site. The first (second) column contain "
                              +"the number of basecalls that support the reference gentype (variant) for nucleotides "
                              "in the sequence context.")

     call_mc_opt.add_argument("--unmethylated-control",
                              type=str,
                              default=None,
                              help="name of the chromosome/region that you want to use to estimate "
                              + "the non-conversion rate of your sample, or the non-conversion rate "
                              + "you would like to use. Consequently, control is either a string, or "
                              + "a decimal. If control is a string then it should be in the following "
                              + "format: \"chrom:start-end\". If you would like to specify an entire "
                              + "chromosome simply use \"chrom:\"")
     
     call_mc_opt.add_argument("--binom-test",
                              type=str2bool,
                              default=False,
                              help="Indicates that you would like to perform a binomial test on each cytosine "
                              +"to delineate cytosines that are significantly methylated than noise due to "
                              +"the failure of bisulfite conversion.")
     
     call_mc_opt.add_argument("--sig-cutoff",
                              type=float,
                              default=.01,
                              help="float indicating the adjusted p-value cutoff you wish to use for "
                              + "determining whether or not a site is methylated")

     call_mc_opt.add_argument("--min-mapq",
                              type=int,
                              default=30,
                              help="Minimum MAPQ for reads to be included.")
     
     call_mc_opt.add_argument("--min-cov",
                              type=int,
                              default=0,
                              help="Integer indicating the minimum number of reads for a site to be tested.")
          
     call_mc_opt.add_argument("--min-base-quality",
                              type=int,
                              default=1,
                              help="Integer indicating the minimum PHRED quality score for a base to be "
                              +"included in the mpileup file (and subsequently to be considered for "
                              +"methylation calling).")

     call_mc_opt.add_argument("--keep-temp-files",
                              type=str2bool,
                              default=False,
                              help="Boolean indicating that you would like to keep the intermediate files "
                              +"generated by this function. This can be useful for debugging, but in general "
                              +"should be left False.")

def add_get_methylation_level_subparser(subparsers):
     add_mc = subparsers.add_parser("add-methylation-level",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     help="Get methylation level of genomic regions")

     add_mc_req = add_mc.add_argument_group("required inputs")
     add_mc_req.add_argument("--input-tsv-file",
                             type=str,
                             help="A tab-separate file that specifies genomic intervals. The file contains a header."
                             +"First three columns are required to be chromosome, start and end, which are "
                             +"1-based cooridates. It may contain additional column(s). ")

     add_mc_req.add_argument("--output-file",
                             type=str,
                             required=True,
                             help="Name of output file")

     add_mc_req.add_argument("--allc-files",
                             type=str,
                             nargs="+",
                             required=True,
                             help="List of allc files.")

     add_mc_opt = add_mc.add_argument_group("optional inputs")
     add_mc_opt.add_argument("--samples",
                             type=str,
                             nargs="+",
                             default=None,
                             help="List of space separated samples matching allc files. By default "
                             +"sample names will be inferred from allc filenames")

     add_mc_opt.add_argument("--mc-type",
                             type=str,
                             nargs="+",
                             default=["CGN"],
                             help="List of space separated mc nucleotide contexts for "
                             + "which you want to look for DMRs. These classifications "
                             + "may use the wildcards H (indicating anything but a G) and "
                             + "N (indicating any nucleotide).")

     add_mc_opt.add_argument("--extra-info",
                             type=str2bool,
                             default=False,
                             help="Boolean to indicate whether to generate two output extra files with "
                             +"the total basecalls and covered sites in each of the regions.")

     add_mc_opt.add_argument("--num-procs",
                             type=int,
                             default=1,
                             help="Number of processors you wish to use to parallelize this function")
     
     add_mc_opt.add_argument("--min-cov",
                             type=int,
                             default=0,
                             help="Minimum coverage for a site to be included")
     
     add_mc_opt.add_argument("--max-cov",
                             type=int,
                             default=None,
                             help="Maximum coverage for a site to be included. By default this cutoff is not applied.")
     
     add_mc_opt.add_argument("--buffer-line-number",
                             type=int,
                             default=100000,
                             help="size of buffer for reads to be written on hard drive.")

     add_mc_opt.add_argument("--input-no-header",
                             type=str2bool,
                             default=False,
                             help="Indicating whether input tsv file contains a header. If this is set to "
                             +"True, a header will be automatically generated in the output file.")

def add_allc2bw_subparser(subparsers):
     allc2bw = subparsers.add_parser("allc-to-bigwig",
                                         formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                         help="Get bigwig file from allc file")

     allc2bw_req = allc2bw.add_argument_group("required inputs")
     allc2bw_req.add_argument("--allc-file",
                              type=str,
                              help="input allc file to be converted to bigwig format")

     allc2bw_req.add_argument("--output-file",
                              type=str,
                              required=True,
                              help="Name of output file")

     allc2bw_req.add_argument("--ref-fasta",
                              type=str,
                              required=True,
                              help="string indicating the path to a fasta file containing the "
                              + "genome sequences")

     allc2bw_opt = allc2bw.add_argument_group("optional inputs")
     allc2bw_opt.add_argument("--mc-type",
                              type=str,
                              nargs="+",
                              default=["CGN"],
                              help="List of space separated mc nucleotide contexts for "
                              + "which you want to look for DMRs. These classifications "
                              + "may use the wildcards H (indicating anything but a G) and "
                              + "N (indicating any nucleotide).")

     allc2bw_opt.add_argument("--bin-size",
                              type=int,
                              default=100,
                              help="Genomic bin size for calculating methylation level")

     allc2bw_opt.add_argument("--min-bin-sites",
                              type=int,
                              default=0,
                              help="Minimum sites in a bin for it to be included.")
     
     allc2bw_opt.add_argument("--min-bin-cov",
                              type=int,
                              default=0,
                              help="Minimum total coverage of all sites in a bin for methylation level "
                              +"to be calculated.")

     allc2bw_opt.add_argument("--min-site-cov",
                              type=int,
                              default=0,
                              help="Minimum total coverage of a site for it to be included.")

     allc2bw_opt.add_argument("--max-site-cov",
                              type=int,
                              default=None,
                              help="Maximum total coverage of a site for it to be included.")
     
     allc2bw_opt.add_argument("--path-to-wigToBigWig",
                              type=str,
                              default="",
                              help="Path to wigToBigWig executable ")
     
     allc2bw_opt.add_argument("--path-to-samtools",
                              type=str,
                              default="",
                              help="Path to samtools installation")

     allc2bw_opt.add_argument("--remove-chr-prefix",
                              type=str2bool,
                              default=True,
                              help="Boolean indicates whether to remove \"chr\" in the chromosome names in "
                              +"genome sequence file to match chromosome names in input allc file.")

     allc2bw_opt.add_argument("--add-chr-prefix",
                              type=str2bool,
                              default=True,
                              help="Boolean indicates whether to add \"chr\" in the chromosome names in "
                              +"input allc file to match chromosome names in genome sequence file. This option "
                              +"overrides --remove-chr-prefix.")

def add_merge_allc_subparser(subparsers):
     merge_allc = subparsers.add_parser("merge-allc",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     help="Merge allc files")

     merge_allc_req = merge_allc.add_argument_group("required inputs")
     merge_allc_req.add_argument("--allc-files",
                             type=str,
                             nargs="+",
                             required=True,
                             help="List of allc files to merge.")

     merge_allc_req.add_argument("--output-file",
                                 type=str,
                                 required=True,
                                 help="String indicating the name of output file")
     
     merge_allc_opt = merge_allc.add_argument_group("optional inputs")
     merge_allc_opt.add_argument("--num-procs",
                                 type=int,
                                 default=1,
                                 help="Number of processors to use")

     merge_allc_opt.add_argument("--compress-output",
                                 type=str2bool,
                                 default=True,
                                 help="Boolean indicating whether to compress (by gzip) the final output")

     merge_allc_opt.add_argument("--skip-snp-info",
                                 type=str2bool,
                                 default=True,
                                 help="Boolean indicating whether to skip the merging of SNP information")
     
     merge_allc_opt.add_argument("--mini-batch",
                                 type=int,
                                 default=100,
                                 help="The maximum number of allc files to be merged at the same time. Since "
                                 +"OS or python may limit the number of files that can be open at once, value "
                                 +"larger than 200 is not recommended")

def add_index_allc_subparser(subparsers):
     index_allc = subparsers.add_parser("index-allc",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     help="Index allc files")

     index_allc_req = index_allc.add_argument_group("required inputs")
     index_allc_req.add_argument("--allc-files",
                             type=str,
                             nargs="+",
                             required=True,
                             help="List of allc files to index.")
     
     index_allc_opt = index_allc.add_argument_group("optional inputs")
     index_allc_opt.add_argument("--num-procs",
                                 type=int,
                                 default=1,
                                 help="Number of processors to use")

     index_allc_opt.add_argument("--reindex",
                                 type=str2bool,
                                 default=True,
                                 help="Boolean indicating whether to index allc files whose "
                                 +"index files already exist.")

def add_filter_allc_subparser(subparsers):
     filter_allc = subparsers.add_parser("filter-allc",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     help="Filter allc file")

     filter_allc_req = filter_allc.add_argument_group("required inputs")
     filter_allc_req.add_argument("--allc-files",
                                  type=str,
                                  required=True,
                                  nargs="+",
                                  help="allc files to filter.")

     filter_allc_req.add_argument("--output-files",
                                  type=str,
                                  required=True,
                                  nargs="+",
                                  help="Name of output files. Each output file matches each allc file.")

     filter_allc_opt = filter_allc.add_argument_group("optional inputs")

     filter_allc_opt.add_argument("--num-procs",
                                  type=int,
                                  default=1,
                                  help="Number of processors you wish to use to parallelize this function")

     filter_allc_opt.add_argument("--mc-type",
                                  type=str,
                                  nargs="+",
                                  default=None,
                                  help="List of space separated cytosine nucleotide contexts for "
                                  + "sites to be included in output file. These classifications "
                                  + "may use the wildcards H (indicating anything but a G) and "
                                  + "N (indicating any nucleotide).")

     filter_allc_opt.add_argument("--min-cov",
                                  type=int,
                                  default=0,
                                  help="Minimum number of reads that must cover a site for it to be "
                                  + "included in the output file.")

     filter_allc_opt.add_argument("--max-cov",
                                  type=int,
                                  default=None,
                                  help="Maximum number of reads that must cover a site for it to be "
                                  + "included in the output file. By default this cutoff is not applied.")

     filter_allc_opt.add_argument("--max-mismatch",
                                  type=int,
                                  nargs="+",
                                  default=None,
                                  help="Maximum numbers of mismatch basecalls allowed in each nucleotide in "
                                  +"the sequence context of a site for it to be included in output file. If "
                                  +"the sequence context has three nucleotides, an example of this option is \"0 1 2\". "
                                  +"It requires no mismatch basecall at the first nucleotide, at most one mismatch "
                                  +"basecall at the second nucleotide, and at most two at the third nucleotide for a site "
                                  +"to be reported.")

     filter_allc_opt.add_argument("--max-mismatch-frac",
                                  type=float,
                                  nargs="+",
                                  default=None,
                                  help="Maximum fraction of mismatch basecalls out of unambiguous basecalls allowed "
                                  +"in each nucleotide in the sequence context of a site for it to be included "
                                  +" in output file. If the sequence context has three nucleotides, an example "
                                  +"of this option is \"0 0 0.1\". It requires no mismatch basecall at the first "
                                  +"and second nucleotide, and at most 10%% mismatches out of unambiguous basecalls "
                                  +"at the third nucleotide for a site to be reported.")
     
     filter_allc_opt.add_argument("--compress-output",
                                  type=str2bool,
                                  default=True,
                                  help="Boolean indicating whether to compress (by gzip) the final output")

     filter_allc_opt.add_argument("--chroms",
                                  type=str,
                                  nargs="+",
                                  default=None,
                                  help="Space separated listing of chromosomes to be included in the output. "
                                  +"By default, data of all chromosomes in input allc file will be included.")

def add_test_allc_subparser(subparsers):
     test_allc = subparsers.add_parser("test-allc",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     help="Binomial test on allc file")

     test_allc_req = test_allc.add_argument_group("required inputs")
     test_allc_req.add_argument("--allc-file",
                                  type=str,
                                  required=True,
                                  help="allc file to be tested.")

     test_allc_req.add_argument("--sample",
                                  type=str,
                                  required=True,
                                  help="sample name")

     test_allc_req.add_argument("--unmethylated-control",
                                type=str,
                                help="name of the chromosome/region that you want to use to estimate "
                                + "the non-conversion rate of your sample, or the non-conversion rate "
                                + "you would like to use. Consequently, control is either a string, or "
                                + "a decimal. If control is a string then it should be in the following "
                                + "format: \"chrom:start-end\". If you would like to specify an entire "
                                + "chromosome simply use \"chrom:\"")

     test_allc_opt = test_allc.add_argument_group("optional inputs")

     test_allc_opt.add_argument("--path-to-output",
                                type=str,
                                default="",
                                help="Path to a directory where you would like the output to be stored. "
                                + "The default is the same directory as the input fastqs.")

     test_allc_opt.add_argument("--num-procs",
                                  type=int,
                                  default=1,
                                  help="Number of processors you wish to use to parallelize this function")

     test_allc_opt.add_argument("--min-cov",
                                  type=int,
                                  default=2,
                                  help="Minimum number of reads that must cover a site for it to be "
                                  + "tested.")

     test_allc_opt.add_argument("--compress-output",
                                  type=str2bool,
                                  default=True,
                                  help="Boolean indicating whether to compress (by gzip) the final output")

     test_allc_opt.add_argument("--sig-cutoff",
                                type=float,
                                default=.01,
                                help="Float indicating at what FDR you want to consider a result "
                                + "significant.")

     test_allc_opt.add_argument("--sort-mem",
                                type=str,
                                default="500M",
                                help="Parameter to pass to unix sort with -S/--buffer-size command")

     test_allc_opt.add_argument("--remove-chr-prefix",
                                type=str2bool,
                                default=True,
                                help="Boolean indicates whether to remove in the final output the \"chr\" prefix "
                                +"in the chromosome name")

def str2bool(v):
     ## adapted from the answer by Maxim at
     ## https://stackoverflow.com/questions/15008758/parsing-boolean-values-with-argparse
     if v.lower() in ('yes', 'true', 't', 'y', '1'):
          return(True)
     elif v.lower() in ('no', 'false', 'f', 'n', '0'):
          return(False)
     else:
          raise argparse.ArgumentTypeError('Boolean value expected.')

              
if __name__ == "__main__":
     parse_args()
