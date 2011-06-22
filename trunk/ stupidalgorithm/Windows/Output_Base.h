#ifndef STUPIDALGO_PSO_OUTPUT_INTERFACE
#define STUPIDALGO_PSO_OUTPUT_INTERFACE

#include <iostream>
#include <boost/shared_ptr.hpp>
#include "Algo_DataStruct.h"
#include "com_alg.h"

class output_base
{
public:
	void print_run_title(std::ostream &os);
	inline void print_run_times(std::ostream &os,int m_cur_run) { os<<"\n"<<"run="<<m_cur_run; }
	void print_sep_line(std::ostream &os);
	void print_gen_stat(std::ostream &os,int cur_gen,const alg_stat &m_alg_stat);
	void print_run_stat_title(std::ostream &os);
	void print_run_stat(std::ostream &os,const alg_stat &m_alg_stat,int max_run);
	void print_gen_avg_val(std::ostream &os,const d_array &gen_vals);
	void print_gen_avg_vec(std::ostream &os,const d_mat &gen_vec);
	void print_best_x(std::ostream &os,const individual &bst_ind);
	void print_all_gbest_val(std::ostream &os,const d_array &gbest_vec);
	void print_avg_time(std::ostream &os,double sec);
protected:
	void write_initial_pop(std::ostream &os,const population &pop);
	inline void print_gen_val(std::ostream &os,int cur_gen,double value) { os<<cur_gen<<"\t"<<value<<"\n"; }
	void print_gen_vec(std::ostream &os,int cur_gen,const d_array &vec);
};// end class output_base declaration

void print_algo_stat_title(std::ostream &os,bool max_gen_set);
void print_para_stat_title(std::ostream &os,std::string para_var_list,std::string para_num_str,bool max_gen_set);
void print_common_filelist_title(std::ostream &os);
inline void print_file_name( std::ostream &os,int func_type,int algo_type,std::string file_name ) { os<<func_type<<"\t"<<algo_type<<"\t"<<file_name<<"\n"; }
#endif