'''
Created on Oct 23, 2011

@author: Matt
'''
import sys
import subprocess
try:
    import unittest
except:
    sys.exit("DMRfind_test requires unittest")
try:
    from methylpy.DMRfind import DMRfind, collapse_dmr_windows,run_rms_tests
    from methylpy.DMRfind import __file__ as module_dir

except:
    sys.exit("DMRfind_test requires the methylpy.DMRfind module")
class DMRfind_tests(unittest.TestCase):
    
    def test_dmrfind_multi_processor(self):
        devnull = open('/dev/null', 'w')
        try:
            subprocess.check_call(["rm", "test_dmr_find_output_rms_results.tsv", "test_dmr_find_output_rms_results_collapsed.tsv"], stderr=devnull)
        except:
            pass
        DMRfind(["CNN"],{"chr1":[0,50000000]},["sampleA","sampleB"],"test_files/",save_result="test_dmr_find_output",num_procs=2,seed=1, dmr_max_dist=10,num_sig_tests=100,num_sims=10000)
        correct_file = open("test_files/DMRfind_correct_output_1_proc.tsv",'r')
        unit_test_file = open("test_dmr_find_output_rms_results_collapsed.tsv",'r')
        line1 = correct_file.readline()
        line2 = unit_test_file.readline()
        while line1 or line2:
            self.assertFalse(line1 == "" and line2,"Too many entries in output file")
            self.assertFalse(line2 == "" and line1,"Too few entries in output file")
            self.assertEqual(line1,line2,"Line from output file:\n"+line2+" is different from the line from the correct file:\n"+line1) 
            line1 = correct_file.readline()
            line2 = unit_test_file.readline()   
        correct_file.close()
        unit_test_file.close()
        subprocess.check_call(["rm", "test_dmr_find_output_rms_results.tsv", "test_dmr_find_output_rms_results_collapsed.tsv"], stderr=devnull)
        devnull.close()
    
    def test_dmrfind_single_processor(self):
        devnull = open('/dev/null', 'w')
        try:
            subprocess.check_call(["rm", "test_dmr_find_output_rms_results.tsv", "test_dmr_find_output_rms_results_collapsed.tsv"], stderr=devnull)
        except:
            pass
        DMRfind(["CNN"],{"chr1":[0,50000000]},["sampleA","sampleB"],"test_files/", save_result="test_dmr_find_output",num_procs=1,seed=1, dmr_max_dist=10,num_sig_tests=100,num_sims=10000)
        correct_file = open("test_files/DMRfind_correct_output_1_proc.tsv",'r')
        unit_test_file = open("test_dmr_find_output_rms_results_collapsed.tsv",'r')
        line1 = correct_file.readline()
        line2 = unit_test_file.readline()
        while line1 or line2:
            self.assertFalse(line1 == "" and line2,"Too many entries in output file")
            self.assertFalse(line2 == "" and line1,"Too few entries in output file")
            self.assertEqual(line1,line2,"Line from output file:\n"+line2+" is different from the line from the correct file:\n"+line1) 
            line1 = correct_file.readline()
            line2 = unit_test_file.readline()   
        correct_file.close()
        unit_test_file.close()
        subprocess.check_call(["rm", "test_dmr_find_output_rms_results.tsv", "test_dmr_find_output_rms_results_collapsed.tsv"], stderr=devnull)
        devnull.close()
    
    def test_dmrfind_sample_category_1(self):
        devnull = open('/dev/null', 'w')
        try:
            subprocess.check_call(["rm", "test_dmr_find_output_rms_results.tsv", "test_dmr_find_output_rms_results_collapsed.tsv"], stderr=devnull)
        except:
            pass
        DMRfind(["CNN"],{"chr1":[0,50000000]},["sampleA","sampleB"],"test_files/",save_result="test_dmr_find_output",num_procs=1,seed=1,collapse_samples=["sampleA","sampleB"], sample_category=[0,1], min_cluster=2, dmr_max_dist=10, sig_cutoff=0.99,num_sig_tests=100,num_sims=10000)
        unit_test_file = open("test_dmr_find_output_rms_results_collapsed.tsv",'r')
        correct_file = open("test_files/DMRfind_correct_output_sample_category_1.tsv",'r')
        line1 = correct_file.readline()
        line2 = unit_test_file.readline()
        while line1 or line2:
            self.assertFalse(line1 == "" and line2,"Too many entries in output file")
            self.assertFalse(line2 == "" and line1,"Too few entries in output file")
            self.assertEqual(line1,line2,"Line from output file:\n"+line2+" is different from the line from the correct file:\n"+line1) 
            line1 = correct_file.readline()
            line2 = unit_test_file.readline()   
        unit_test_file.close()
        correct_file.close()
        subprocess.check_call(["rm", "test_dmr_find_output_rms_results.tsv", "test_dmr_find_output_rms_results_collapsed.tsv"], stderr=devnull)
        devnull.close()
    
    def test_dmrfind_sample_category_2(self):
        devnull = open('/dev/null', 'w')
        try:
            subprocess.check_call(["rm", "test_dmr_find_output_rms_results.tsv", "test_dmr_find_output_rms_results_collapsed.tsv"], stderr=devnull)
        except:
            pass
        DMRfind(["CNN"],{"chr1":[0,50000000]},["sampleA","sampleB"],"test_files/",save_result="test_dmr_find_output",num_procs=1,seed=1,collapse_samples=["sampleA","sampleB"], sample_category=[0,1], min_cluster=2, dmr_max_dist=10, sig_cutoff=0.99,num_sig_tests=100,num_sims=10000)        
        unit_test_file = open("test_dmr_find_output_rms_results_collapsed.tsv",'r')
        correct_file = open("test_files/DMRfind_correct_output_sample_category_2.tsv",'r')
        line1 = correct_file.readline()
        line2 = unit_test_file.readline()
        while line1 or line2:
            self.assertFalse(line1 == "" and line2,"Too many entries in output file")
            self.assertFalse(line2 == "" and line1,"Too few entries in output file")
            self.assertEqual(line1,line2,"Line from output file:\n"+line2+" is different from the line from the correct file:\n"+line1) 
            line1 = correct_file.readline()
            line2 = unit_test_file.readline()   
        unit_test_file.close()
        correct_file.close()
        subprocess.check_call(["rm", "test_dmr_find_output_rms_results.tsv", "test_dmr_find_output_rms_results_collapsed.tsv"], stderr=devnull)
        devnull.close()
    
    def test_run_rms_tests(self):
        devnull = open('/dev/null', 'w')
        subprocess.check_call(["cp", module_dir[:module_dir.rfind("/")]+"/test/test_files/allc_sampleA_1.tsv", module_dir[:module_dir.rfind("/")]+"/test/test_files/temp_rms_test1"])
        subprocess.check_call(["cp", module_dir[:module_dir.rfind("/")]+"/test/test_files/allc_sampleB_1.tsv", module_dir[:module_dir.rfind("/")]+"/test/test_files/temp_rms_test2"])
        
        files = [module_dir[:module_dir.rfind("/")]+"/test/test_files/temp_rms_test1", module_dir[:module_dir.rfind("/")]+"/test/test_files/temp_rms_test2"]
        samples = ["sampleA","sampleB"]
        output = "run_rms_tests_test_output.tsv"
        run_rms_tests(files,output,samples,min_cov = 0,num_sims="10000",num_sig_tests="100",seed=1)
        correct_file = open("test_files/run_rms_tests_correct_output.tsv",'r')
        unit_test_file = open("run_rms_tests_test_output.tsv",'r')
        line1 = correct_file.readline()
        line2 = unit_test_file.readline()
        while line1 or line2:
            self.assertFalse(line1 == "" and line2,"Too many entries in output file")
            self.assertFalse(line2 == "" and line1,"Too few entries in output file")
            self.assertEqual(line1,line2,"Line from output file:\n"+line2+" is different from the line from the correct file:\n"+line1) 
            line1 = correct_file.readline()
            line2 = unit_test_file.readline()   
        correct_file.close()
        unit_test_file.close()
        subprocess.check_call(["rm", module_dir[:module_dir.rfind("/")]+"/test/run_rms_tests_test_output.tsv"], stderr = devnull)
        devnull.close()
    #def test_collapse_dmr_windows(self):
    
if __name__ == "__main__":
    unittest.main()
