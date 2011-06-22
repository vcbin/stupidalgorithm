#ifndef STUPIDALGO_ALL_ALGORITHM_SCHEDULER
#define STUPIDALGO_ALL_ALGORITHM_SCHEDULER

#include <iostream>
#include <strstream>

#include "FuncDef.h"
#include "com_alg.h"
#include "para_base.h"

class scheduler
{
public:
	// ctor
	scheduler(std::string conf_path):
		m_batch_conf_path(conf_path) {}

	  int batch_run(int max_thread);
private:
	void algo_spec_output(int func_type,int algo_type,const boost::shared_ptr<alg_base> &pAlg,std::string str_common,std::string str_seq_desc);
	// functions for output description text
	std::string gen_func_eval_name(const boost::shared_ptr<func_base> &pFunc);
	std::string gen_out_prefix_common(int algo_type,const boost::shared_ptr<func_base> &pFunc,const boost::shared_ptr<para_base> &ppara);
	std::string gen_out_prefix_para(const boost::shared_ptr<para_base> &ppara,std::string delimiter);
	std::string gen_out_para_val_str(const boost::shared_ptr<para_base> &ppara);
	std::string gen_seq_var_name_str(const boost::shared_ptr<para_base> &ppara);
	std::string gen_seq_val_str(const boost::shared_ptr<para_base> &ppara);
	std::string gen_seq_var_num_str(const boost::shared_ptr<para_base> &ppara);
	void output_running_prompt(int counter,const boost::shared_ptr<para_base> &ppara);
	// functions for output stat file names

	void gen_common_file_path(common_path &com_out_path,std::string str_common,std::string str_seq_desc);
	std::string gen_stat_file_name(std::string str_common,std::string str_seq_desc,std::string str_append);
	std::string gen_all_bst_val_file_name(std::string str_common,std::string str_seq_desc);
	std::string gen_stat_file_name(std::string str_common,std::string str_seq_desc);
	std::string gen_avg_bst_file_name(std::string str_common,std::string str_seq_desc);
	std::string gen_avg_div_file_name(std::string str_common,std::string str_seq_desc);
	std::string gen_avg_rad_file_name(std::string str_common,std::string str_seq_desc);
	std::string gen_para_stat_file_name(std::string str_common);

	std::string gen_avg_vel_div_file_name(std::string str_common,std::string str_seq_desc);
	std::string gen_avg_attr_file_name(std::string str_common,std::string str_seq_desc);
	std::string gen_avg_f_file_name(std::string str_common,std::string str_seq_desc);
	std::string gen_avg_pr_file_name(std::string str_common,std::string str_seq_desc);
	std::string gen_avg_eta_file_name(std::string str_common,std::string str_seq_desc);
	std::string gen_avg_bi_norm_stat_file_name(std::string str_common,std::string str_seq_desc);
	void print_common_filelist(std::ostream &os,int func_type,int algo_type,const common_path &com_out_path);

	std::string m_batch_conf_path;
	std::string m_algo_name;
	std::string m_func_eval_name;
	std::string m_para_out_prefix;
	std::string m_para_val_prefix;

	void print_stat_file(const std::vector<boost::shared_ptr<std::stringstream> >& vec_line);

	struct locked_para_stat_os
	{
		boost::shared_ptr<std::ostream> pos;
		boost::shared_ptr<boost::mutex> pmut;
	};
};

#endif