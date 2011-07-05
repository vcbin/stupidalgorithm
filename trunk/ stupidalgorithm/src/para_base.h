#ifndef STUPIDALGO_BASE_ALL_PARAMETERS
#define STUPIDALGO_BASE_ALL_PARAMETERS

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/tokenizer.hpp>
#include <boost/token_functions.hpp>
#include <boost/algorithm/string.hpp> // split
#include <boost/any.hpp>
#include <map>
#include <sstream> // operator<< definition
#include <iostream>

#include "algo_datastruct.h"
#include "cmdline_para_type.h"

namespace po = boost::program_options;

typedef std::map<std::string,std::vector<boost::any> > num_var_map;// sequential numeric variables mapping

template <typename T>struct type_aux
{
  typedef T type;
};

class para_base
{
public:
	para_base(std::string conf_path)
		:m_conf_path(conf_path),
		general("General options")
	{
		general.add_options()
			("help", "produce a help message")
			("version", "output the version number")

			("func_type",po::value<int>(),"test problem type")
			("value_to_reach",po::value<std::string>(&vtr)->default_value("0"),"Value To Reach(precision of solution quality)")
			("pop_size", po::value<int>(),"population size" )
			("stop_type",po::value<int>(),"stopping criterion type")
			("stop_threshold",po::value<double>(),"stopping criterion associated threshold value")
			("ini_type",po::value<int>(&ini_type)->default_value(1),"initialization type")
			("ini_pop_file",po::value<std::string>(&ini_pop_path)->default_value(""),"associated initialization population file path")
			("ini_range_file",po::value<std::string>(&ini_rng_path)->default_value(""),"associated initialization value ranges file path")
			("bnd_value_type",po::value<int>(&bnd_val_type)->default_value(bndvaltype_all),"boundary value type")
			("bnd_range_file",po::value<std::string>(&bnd_rng_path)->default_value(""),"associated boundary value ranges file path")
			("run",po::value<int>(),"run/trial number")
			("dimension",po::value<std::string>(),"number of real variable/ dimensionality of test problem")
			("val_lower_bound",po::value<double>(),"lower bound of ALL variables")
			("val_upper_bound",po::value<double>(),"upper bound of ALL variables")
			("ini_lower_bound",po::value<double>(),"initialization lower bound of ALL variables")
			("ini_upper_bound",po::value<double>(),"initialization upper bound of ALL variables")
			("out_interval",po::value<int>(&out_interval)->default_value(1),"output interval")
			("q_number",po::value<std::string>(&q_number)->default_value("10"),"orthogonal initialization granularity")

			("num_obj",po::value<int>(&num_obj)->default_value(1),"number of objectives")
			("trunc_type",po::value<std::string>(&trunc_type)->default_value("2"),"truncation type")
			("plot",po::value<bool>(&plot)->default_value(true),"plotting flag")
			("plot_script",po::value<std::string>(&plot_script)->default_value("plot_all_point_2d.p"),"plotting script path")
			;
	}
	virtual ~para_base(){}
	virtual int read_para(std::string conf_path)=0;

	inline	int get_obj_func() {return func_type;}
	inline double get_vtr() { return get_seq_val<double>("vtr"); }
	inline	int get_max_run() {return run;}
	inline	int get_pop_size() {return pop_size;}
	inline	int get_stop_type() {return stop_type;}
	inline	double get_stop_val() {return stop_val;}
	inline	int get_ini_type() {return ini_type;}
	inline int get_bnd_val_type() {return bnd_val_type;}
	inline	std::string get_ini_pop_file_path() {return ini_pop_path;}
	inline	int get_max_gen() {return max_gen;}
	inline	int get_dim() {return get_seq_val<int>("dimension");}
	inline	val_range& get_val_bnd() {return val_bounds;}
	inline	val_range& get_ini_rng() {return ini_bounds;}
	inline int get_out_interval() {return out_interval;}
	inline int get_q_number() {return get_seq_val<int>("q_number");}

	inline int get_obj_num() {return num_obj;}
	inline int get_trunc_type() {return get_seq_val<int>("trunc_type");}
	inline bool get_plot_flag() {return plot;}
	inline std::string get_plot_script() {return plot_script;}
	// parameter value assignment operator
	inline para_base& operator=(const para_base& rhs)
	{
		// boost::options_descriptions is NON-copyable
		m_conf_path=rhs.m_conf_path;
		func_type=rhs.func_type;
		vtr=rhs.vtr;
		run=rhs.run;
		pop_size=rhs.pop_size;
		stop_type=rhs.stop_type;
		stop_val=rhs.stop_val;
		ini_type=rhs.ini_type;
		ini_pop_path=rhs.ini_pop_path;
		ini_rng_path=rhs.ini_rng_path;
		bnd_val_type=rhs.bnd_val_type;
		bnd_rng_path=rhs.bnd_rng_path;
		max_gen=rhs.max_gen;
		dims=rhs.dims;
		val_bounds=rhs.val_bounds;
		ini_bounds=rhs.ini_bounds;
		q_number=rhs.q_number;

		num_obj=rhs.num_obj;
		plot=rhs.plot;
		plot_script=rhs.plot_script;

		m_seq_var_map=rhs.m_seq_var_map;
		m_var_val_idx=rhs.m_var_val_idx;
		m_seq_dim_size=rhs.m_seq_dim_size;
		m_seq_var_name=rhs.m_seq_var_name;

		return *this;
	} // end operator=

	// functions for sequential parameters(of numeric type int/double wrapped by boost::any)
	inline bool has_seq_value() {return !m_seq_var_name.empty();}
	inline const std::vector<std::string>&get_seq_var_name_vec() {return m_seq_var_name;}
	inline int get_cur_seq_val_idx(const std::string& para) { return m_var_val_idx[para]; }
	inline boost::any get_cur_seq_val(const std::string& para)
	{
		int cur_idx=get_cur_seq_val_idx(para);
		return m_seq_var_map[para][cur_idx];
	}

	inline void set_cur_seq_val_idx(const std::string& para,int idx) { m_var_val_idx[para]=idx; }
	inline void set_cur_seq_val(const std::string& para,double val)
	{
		int cur_idx=get_cur_seq_val_idx(para);
		m_seq_var_map[para][cur_idx]=val;
	}
	inline const std::vector<int>& get_seq_dim_size_vec() { return m_seq_dim_size; }
	/*inline void set_cur_seq_val(const std::string& para,int cur_idx,double val)
	{
		m_seq_var_map[para][cur_idx]=val;
		set_cur_seq_val_idx(para,cur_idx);
	}*/
	template <typename T>
	inline T get_seq_val(const std::string& para)
	{
		return get_val(para,type_aux<T>());
	}

protected:
	template <typename T>
	inline T get_val(const std::string& para,type_aux<T> type)
	{
		T ret_val;
		boost::any any_val=get_cur_seq_val(para);
		if ( any_val.type()==typeid(T) )
			ret_val=boost::any_cast<T>(any_val);
		else // target type and parsed type doesn't match,try to correct it
		{
			// can't convert to target type,throw bad_cast error
			std::stringstream str_err;
			str_err<<"\n"
				<<"Bad Cast In para_base!"
				<<"required type is"
				<<" "
				<<typeid(T).name()
				<<","
				<<"while exist type is"
				<<" "
				<<any_val.type().name()
				<<"."
				<<"\n";
			std::cout<<str_err.str();
			throw std::bad_cast();
		}
		return ret_val;
	}// end function get_seq_val

	inline int get_val(const std::string& para,type_aux<int>)
	{
		int ret_val;
		boost::any any_val=get_cur_seq_val(para);
		if ( any_val.type()==typeid(int) )
			ret_val=boost::any_cast<int>(any_val);
		else // target type and parsed type doesn't match,try to correct it
		{
			try
			{
				ret_val=static_cast<int>(boost::any_cast<double>(get_cur_seq_val(para)));
			}
			catch (std::bad_cast &bc)
			{
				std::stringstream str_err;
				str_err<<"\n"
					<<"Bad Cast In para_base!"
					<<"required type is"
					<<" "
					<<typeid(int).name()
					<<","
					<<"while exist type is"
					<<" "
					<<any_val.type().name()
					<<"."
					<<"\n";
				std::cout<<str_err.str();
				throw bc;
			}
		}
		return ret_val;
	}// end function get_seq_val<int>

	inline double get_val(const std::string& para,type_aux<double>)
	{
		double ret_val;
		boost::any any_val=get_cur_seq_val(para);
		if ( any_val.type()==typeid(double) )
			ret_val=boost::any_cast<double>(any_val);
		else // target type and parsed type doesn't match,try to correct it
		{
			try
			{
				ret_val=static_cast<double>(boost::any_cast<int>(get_cur_seq_val(para)));
			}
			catch (std::bad_cast &bc)
			{
				std::stringstream str_err;
				str_err<<"\n"
					<<"Bad Cast In para_base!"
					<<"required type is"
					<<" "
					<<typeid(double).name()
					<<","
					<<"while exist type is"
					<<" "
					<<any_val.type().name()
					<<"."
					<<"\n";
				std::cout<<str_err.str();
				throw bc;
			}
		}
		return ret_val;
	}// end function get_seq_val<double>

	void load_range_file(std::string rng_path);
	void load_ini_range() { return load_range_file(ini_rng_path); }
	void load_bnd_val_range() { return load_range_file(bnd_rng_path); }

	int parse_seq_var(po::variables_map &m_ext_vm,std::string opt_name,std::string literal_va);

	std::string m_conf_path;

	// int alg_type;// pso algorithm type
	int func_type;// objective function type
	std::string vtr;
	int run;// run/trial number
	int pop_size;// population size
	int stop_type;// stop criterion type
	double stop_val;// stop value:max_generation or max_evaluation or max_stagnation
	int ini_type;// initialization type
	std::string ini_pop_path;// OPTIONAL,associated initialization population file path
	std::string ini_rng_path;// OPTIONAL,associated initialization value ranges file path
	int bnd_val_type;// boundary value type
	std::string bnd_rng_path;// OPTIONAL,associated boundary value ranges file path
	int max_gen;// max generation number
	std::string dims;// number of real variable
	val_range val_bounds;// value boundry of dimensions of every variable,for EXTENSION purpose
	val_range ini_bounds;// value boundry of dimensions of every variable,for asymmetric initialization
	int out_interval;// output interval
	std::string q_number;// orthogonal initialization granularity

	int num_obj;// number of test problem objects
	std::string trunc_type;// truncation type
	bool plot;// plot flag
	std::string plot_script;// plotting script path

	num_var_map m_seq_var_map;// PARSED sequential parameter valueS of numerical type
	std::map<std::string,int> m_var_val_idx;// current value indicator of sequential value
	std::vector<int> m_seq_dim_size;// size of every sequential value
	std::vector<std::string> m_seq_var_name;// auxiliary name array of sequential value
	po::options_description general;

	po::variables_map m_ext_vm;
};// end class para_base

#endif
