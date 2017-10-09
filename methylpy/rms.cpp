#include <numeric>
#include <vector>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string>
#include <fstream>
#include <map>
#include <algorithm>
#include <functional>
#include <cctype>
#include <locale>
#include <sstream>
#include <time.h>
#include <iomanip>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
using namespace std;

int rms_test(vector< vector<int> > obs,float* stat, vector< vector<float> >* residuals,vector< vector<float> >* exp_prob);
float calculate_rms_pvalue(float stat,int num_sims,int total,vector< vector<float> > exp_prob, int num_sig_tests, int* num_positive_tests, int* num_sims_run,int seed);
static inline std::string &rtrim(std::string &s);
int run_rms_tests(vector<string> files, string output, vector<string> samples, int min_cov, int num_sims, int num_sig_tests,int seed);
std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems);
std::vector<std::string> split(const std::string &s, char delim);

static inline std::string &rtrim(std::string &s) {
		/* This function was taken from Evan Teran's answer at:
		 * http://stackoverflow.com/questions/216823/whats-the-best-way-to-trim-stdstring
		 */
        s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
        return s;
}

std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while(std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}


std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    return split(s, delim, elems);
}

int rms_test(vector< vector<int> > obs,float* stat, vector< vector<float> >* residuals,vector< vector<float> >* exp_prob){


	int num_sims;
	int nrow = obs.size();
	int ncol = obs[0].size();
	vector< vector<float> > exp(nrow,vector<float>(2,0));
	vector< vector<float> > resid_denom(nrow,vector<float>(2,0));
	int size = nrow*ncol;
	float total = 0.0;
	vector<double> rowsum(nrow);
	vector<double> colsum(nrow);
	float row_std_resid;
	for(int row_ind=0;row_ind<nrow;row_ind++){
		for(int col_ind=0;col_ind<ncol;col_ind++){
			rowsum[row_ind]+=obs[row_ind][col_ind];
			colsum[col_ind]+=obs[row_ind][col_ind];
			total+=static_cast<float>(obs[row_ind][col_ind]);
		}
	}
	//total = static_cast<float>(total);
	for(int row_ind=0;row_ind<nrow;row_ind++){
		row_std_resid = (1.0-(rowsum[row_ind]/total));
		for(int col_ind=0;col_ind<ncol;col_ind++){
			exp[row_ind][col_ind] = (rowsum[row_ind]*colsum[col_ind]) / total;
			(*exp_prob)[row_ind][col_ind] = exp[row_ind][col_ind] / total;
	        resid_denom[row_ind][col_ind] = pow((exp[row_ind][col_ind] * row_std_resid * (1-(colsum[col_ind]/total))),0.5);
	        *stat += pow((obs[row_ind][col_ind] - exp[row_ind][col_ind]),2);
	        if(resid_denom[row_ind][col_ind] == 0){
	        	(*residuals)[row_ind][col_ind] = 0;
	        }else{
	        	(*residuals)[row_ind][col_ind] = (obs[row_ind][col_ind] - exp[row_ind][col_ind]) / resid_denom[row_ind][col_ind];
	        }
		}
	}
    *stat = pow((*stat / (ncol*nrow)),0.5);

}
float calculate_rms_pvalue(float stat,int num_sims,int total,vector< vector<float> > exp_prob, int num_sig_tests, int* num_positive_tests, int* num_sims_run,int seed=-1){
	//calculate p-value

	int num_samples = exp_prob.size();
	int num_boxes = num_samples * 2;

	//set randomization seed. Need this control for reproducibility
	const gsl_rng_type *T = gsl_rng_default;
	gsl_rng *r = gsl_rng_alloc (T);
	if(seed == -1){
		gsl_rng_set(r,time(NULL));
		//srand ( time(NULL) );
	} else {
		gsl_rng_set(r,seed);
		//srand( seed );
	}

	double distribution[num_boxes];
	for(int i=0;i<num_boxes;i++){
		distribution[i] = exp_prob[i/2][i%2];

	}

	gsl_ran_discrete_t *sample_struct = gsl_ran_discrete_preproc (num_boxes, distribution);

	float pvalue;
	//Perform at least min_tests before letting the adaptive pvalue cutoff kick in
	vector< vector<int> > obs(num_samples,vector<int>(2,0));
	vector< vector<float> > residuals(num_samples,vector<float>(2));
	vector< vector<float> > new_exp_prob(num_samples,vector<float>(2));
	float new_stat;




	for(*num_sims_run=1;(*num_sims_run)<=num_sims;(*num_sims_run)++){
		new_stat = 0;
		for(int j=0;j<total;j++){
			int rand_num =gsl_ran_discrete(r,sample_struct);
			obs[rand_num/2][rand_num%2]++;
			/*
			float rand_num = (rand() % 100) / 100.0;
			for(int row_ind = 0; row_ind < num_samples;row_ind++){
				if(rand_num < exp_prob[row_ind][0]){
					obs[row_ind][0]++;
					break;
				}
				rand_num -= exp_prob[row_ind][0];
				if(rand_num < exp_prob[row_ind][1]){
					obs[row_ind][1]++;
					break;
				}
				rand_num -= exp_prob[row_ind][1];
			}
			*/
		}

		//cout << obs[0][0] <<"\t" << obs[0][1] << endl;
		//cout << obs[1][0] << "\t" << obs[1][1] << endl;
		rms_test(obs,&new_stat,&residuals,&new_exp_prob);
		if(new_stat >= stat){
			*num_positive_tests+= 1;
		}
		for(int row_ind = 0; row_ind < num_samples;row_ind++){
			obs[row_ind][0]=0;
			obs[row_ind][1]=0;
		}
		pvalue = static_cast<float>(*num_positive_tests) / *num_sims_run;
		if(*num_positive_tests >= num_sig_tests){
			return(pvalue);
		}
	}
	//Definition of pvalue from: http://www.stat.iastate.edu/preprint/articles/2009-04.pdf
	//If max number of simulations are done, then:
	//pvalue = # of randomized results more significant + 1 / max simulations
	//this prevents a 0 pvalue
	*num_positive_tests+=1;
	//since the for loop will end at one over num_sims
	*num_sims_run-=1;
	pvalue = static_cast<float>(*num_positive_tests) / num_sims;
	return(pvalue);
}
int run_rms_tests(vector<string> files, string output, vector<string> samples, int min_cov, int num_sims, int num_sig_tests,int seed=-1){
    string buffer_string = "";
    ofstream results;
    results.open(output.c_str());
    int buffer_count = 0;
    int total_tests = 0;
    int index;
    int position;
    string chrom;
    string strand;
    string mc_class;
    char methylated;
    int sample_index;
    int num_samples = samples.size();
    const char* delim = "\t";
    vector<string> lines(num_samples);
    vector< vector<string> > fields(num_samples);
    vector<ifstream*> file_handles(num_samples);
    vector<int> current_positions(num_samples);
    int precision = ceil(log10(num_sims));

    for(int sample_index = 0; sample_index<num_samples;sample_index++){
    	file_handles[sample_index] = new ifstream(files[sample_index].c_str());
    }

    for(int sample_index = 0; sample_index<num_samples;sample_index++){
    	getline(*file_handles[sample_index],lines[sample_index]);
    	rtrim(lines[sample_index]);
    	fields[sample_index] = split(lines[sample_index],*delim);
		if(fields[sample_index].size()>0){
			current_positions[sample_index]=atoi(fields[sample_index][1].c_str());
		}
		else{
			current_positions[sample_index]=1000000000;
		}
    }

    while(1){
    	int lines_length = 0;
    	for(int sample_index = 0; sample_index<num_samples;sample_index++){
    		lines_length+=lines[sample_index].length();
    	}
    	if(lines_length == 0){
    		results.close();
    	    for(int sample_index = 0; sample_index<num_samples;sample_index++){
    	    	(*file_handles[sample_index]).close();
    	    	delete file_handles[sample_index];
    	    }
    	    cout << "Ran " << total_tests << " root mean square tests. Results in " << output << "\n";
    		return total_tests;
    	}

    	int min_pos = *min_element(current_positions.begin(),current_positions.end());
    	vector< vector<int> > counts;
    	vector<string> current_samples;
    	vector<int> mc_list;
    	vector<int> h_list;
    	vector<float> frac_list;
    	int total = 0;

    	//these variables are used later to see if we even need to bother with a test
		//if the max and min frac are the same, all the fracs are the same so what's the
		//point in testing for a difference
		float min_frac = 1.0;
		float max_frac = 0.0;

    	for(int sample_index = 0; sample_index<num_samples;sample_index++){
    		if(lines[sample_index].length()!=0){
    			position = atoi(fields[sample_index][1].c_str());
    			if(position == min_pos && atoi(fields[sample_index][5].c_str()) >= min_cov){
                    chrom = fields[sample_index][0];
                    strand = fields[sample_index][2];
                    mc_class = fields[sample_index][3];
                    int mc;
                    int h;
                    float frac;

                    string methylated = fields[sample_index][6];
                    vector<int> current_counts(2);
                    
		    mc = atoi(fields[sample_index][4].c_str());
		    h = atoi(fields[sample_index][5].c_str());
		    current_counts[0] = mc;
		    current_counts[1] = h-mc;
		    frac = mc / static_cast<float>(h);

                    if(frac > max_frac){
                    	max_frac = frac;
                    }

                    if(frac < min_frac){
                    	min_frac = frac;
                    }
                    mc_list.push_back(mc);
                    h_list.push_back(h);
                    frac_list.push_back(frac);
                    counts.push_back(current_counts);
                    current_samples.push_back(samples[sample_index]);
                    total+=h;
    			} else {
    				mc_list.push_back(-1);
    				h_list.push_back(-1);
    				frac_list.push_back(-1);
    			}
    		} else {
				mc_list.push_back(-1);
				h_list.push_back(-1);
				frac_list.push_back(-1);
    		}
    	}
    	int num_current_samples = current_samples.size();
    	if(num_current_samples > 1 && min_frac != max_frac){
    		vector< vector<float> > residuals(num_current_samples,vector<float>(2));
    		vector< vector<float> > exp_prob(num_current_samples,vector<float>(2));
    		float stat = 0;
    		int num_sims_run = 1;
    		int num_positive_tests = 0;
    		rms_test(counts,&stat,&residuals,&exp_prob);
    		float pvalue = calculate_rms_pvalue(stat,num_sims,total,exp_prob,num_sig_tests,&num_positive_tests,&num_sims_run,seed);

    		if(pvalue == 0){
    			results << chrom << "\t" << min_pos << "\t" << strand << "\t" << mc_class << "\t<" << scientific << 1 / static_cast<float>(num_sims) << fixed;
    		} else{
    			results << chrom << "\t" << min_pos << "\t" << strand << "\t" << mc_class << fixed << setprecision(precision) << "\t"<< pvalue;
    		}
    		for(vector<int>::iterator it=mc_list.begin(); it < mc_list.end(); it++){
    			results << "\t" << *it;
    		}
    		for(vector<int>::iterator it=h_list.begin(); it < h_list.end(); it++){
    			results << "\t" << *it;
    		}
    		for(vector<float>::iterator it=frac_list.begin(); it < frac_list.end(); it++){
    			results << "\t" << *it;
    		}
    		int current_samples_index = 0;
    		for(int sample_index = 0;sample_index < num_samples;sample_index++){
     			if(current_samples_index < current_samples.size() && samples[sample_index] == current_samples[current_samples_index]){
	    			results << "\t" << residuals[current_samples_index][0];
					current_samples_index++;
				}
				else{
					results << "\t-1";
				}
    		}
    		current_samples_index = 0;
    		for(int sample_index = 0;sample_index < num_samples;sample_index++){
     			if(current_samples_index < current_samples.size() && samples[sample_index] == current_samples[current_samples_index]){
	    			results << "\t" << residuals[current_samples_index][1];
					current_samples_index++;
				}
				else{
					results << "\t-1";
				}
    		}
    		results << "\t" << num_positive_tests << "\t" << num_sims_run;
    		results << "\n";
    		total_tests+=1;
    	}
		for(int sample_index = 0; sample_index<num_samples;sample_index++){
			if(lines[sample_index].length()>0 && current_positions[sample_index] == min_pos){
				getline(*file_handles[sample_index],lines[sample_index]);
				rtrim(lines[sample_index]);
			}
			fields[sample_index] = split(lines[sample_index],*delim);
			if(fields[sample_index].size()>0){
				current_positions[sample_index]=atoi(fields[sample_index][1].c_str());
			}
			else{
				current_positions[sample_index]=1000000000;
			}
		}
    }
}
int main(int argc, char* argv[]){
	vector<string> files;
	vector<string> samples;
	const char* delim = ",";
	if(argc == 7){
		files=split(argv[1],*delim);
		samples=split(argv[3],*delim);
		run_rms_tests(files,argv[2],samples,atoi(argv[4]),atoi(argv[5]),atoi(argv[6]));
	}
	else if(argc == 8){
		files=split(argv[1],*delim);
		samples=split(argv[3],*delim);
		run_rms_tests(files,argv[2],samples,atoi(argv[4]),atoi(argv[5]),atoi(argv[6]),atof(argv[7]));
	}
	else{
		cout << "Usage: ./rms.out <chunk files> <output file> <samples> <min_cov> <num_sims> <num_sig_tests> <seed>\n";
	}
}
