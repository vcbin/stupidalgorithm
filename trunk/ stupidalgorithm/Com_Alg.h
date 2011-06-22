#ifndef JERRYLIU_PSO_ALG_BASE
#define JERRYLIU_PSO_ALG_BASE

#include <fstream>
#include <boost/shared_ptr.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/timer.hpp>

#include "FuncDef.h"
#include "para_base.h"
#include "initializer.h"
#include "problem_base.h"
#include "output_base.h"
#include "alg_base.h"
#include "CmdLine_Para_Type.h"

using boost::lock_guard;

// extern boost::mutex run_stat_file_mutex;
extern boost::mutex para_stat_file_mutex;

template <typename real_para>
class com_alg:public alg_base,public problem_base,public output_base
{
public:
	// ctor
	com_alg(std::string conf_path):
	  m_ppara(new real_para(conf_path)) {}
	  com_alg(boost::shared_ptr<real_para> ppara):// set algorithm parameters pointer directly
	  m_ppara(ppara) {}
	  // dtor
	  virtual ~com_alg() {}

	  void initialize()
	  {
		  set_diag_len();

		  // initialize stat variable
		  int stop_type=m_ppara->get_stop_type();
		  if ( stoptype_stag!=stop_type )
			  m_alg_stat.stop_stag=false;
		  else
			  m_alg_stat.stop_stag=true;
		  int max_run=m_ppara->get_max_run();
		  int pop_size=m_ppara->get_pop_size();
		  int num_dims=m_ppara->get_dim();
		  m_alg_stat.alloc_space(max_run,pop_size,num_dims);
		  m_alg_stat.initialize(max_run);
	  }// end function initialize

	  void set_para_val_desc(std::string str_para_val) { m_para_val_desc=str_para_val; }
	  void load_parameter(const boost::shared_ptr<para_base>& ppara)
	  { 
		  try
		  {
			  m_ppara=boost::dynamic_pointer_cast<real_para>(ppara);// downcast base parameters type to algorithm specific parameters type
		  }
		  catch (std::bad_cast &bc)
		  {
			  std::cout<<"bad_cast caught in Load_Parameters:"
				  <<bc.what();
			  throw;
		  }
	  }
	  void set_com_file_path(const common_path &com_out_path) { m_com_out_path=com_out_path; }

	  void print_para_stat(std::ostream &os,bool max_gen_set,boost::shared_ptr<boost::mutex> pmut)
	  {
		  lock_guard<boost::mutex> lock(*pmut);
		  os<<"\n"
			  <<m_para_val_desc;
		  os<<"\t"
			  <<m_alg_stat.bst_ind.obj[0]
		  <<"\t"
			  <<m_alg_stat.wst_ind.obj[0]
		  <<"\t"
			  <<m_alg_stat.run_avg_bst
			  <<"\t"
			  <<m_alg_stat.run_bst_std
			  <<"\t"
			  <<m_alg_stat.conv_count
			  <<"\t"
			  <<m_alg_stat.conv_ratio
			  <<"\t"
			  <<m_alg_stat.avg_conv_eval_num
			  <<"\t"
			  <<m_alg_stat.run_avg_time
			  <<"\t"
			  <<m_alg_stat.all_eval_num
			  <<"\t"
			  <<m_alg_stat.all_ob_num
			  <<"\t"
			  <<m_alg_stat.run_avg_stag_gen;
		  if ( max_gen_set )
			  os<<"\t"
			  <<m_alg_stat.run_avg_gen;
	  }// end function print_para_stat

	  void print_algo_stat( std::ostream &os,int func_type,int algo_type,bool max_gen_set )
	  {
		  // lock_guard<boost::mutex> lock(run_stat_file_mutex);
		  // std::stringstream os;
		  os<<"\n"
			  <<func_type
			  <<"\t"
			  <<algo_type
			  <<"\t"
			  <<m_alg_stat.bst_ind.obj[0]
		  <<"\t"
			  <<m_alg_stat.wst_ind.obj[0]
		  <<"\t"
			  <<m_alg_stat.run_avg_bst
			  <<"\t"
			  <<m_alg_stat.run_bst_std
			  <<"\t"
			  <<m_alg_stat.conv_count
			  <<"\t"
			  <<m_alg_stat.conv_ratio
			  <<"\t"
			  <<m_alg_stat.avg_conv_eval_num
			  <<"\t"
			  <<m_alg_stat.run_avg_time
			  <<"\t"
			  <<m_alg_stat.all_eval_num
			  <<"\t"
			  <<m_alg_stat.all_ob_num
			  <<"\t"
			  <<m_alg_stat.run_avg_stag_gen;
		  if ( max_gen_set )
			  os<<"\t"
			  <<m_alg_stat.run_avg_gen;
		  // os_dst=os.str();
	  }// end function print_algo_stat
protected:
	void set_diag_len()
	{
		int num_dims=m_ppara->get_dim();
		val_range& val_bound=m_ppara->get_val_bnd();
		int i;
		m_diag_len=0.0;
		double diff;
		for ( i=0;i<num_dims;i++ )
		{
			diff=val_bound[i].max_val-val_bound[i].min_val;
			m_diag_len += diff*diff;
		}
		m_diag_len=sqrt(m_diag_len);
	}

	void update_diversity(population &pop)
	{
		int pop_size=m_ppara->get_pop_size();
		int num_dims=m_ppara->get_dim();

		int i;
		// update centroid individual
		m_alg_stat.cen_ind.assign(num_dims,0.0);
		for ( i=0;i<pop_size;i++ )
		{
			m_alg_stat.cen_ind += pop[i].x;
		}// for every particle
		m_alg_stat.cen_ind /=  pop_size;

		m_alg_stat.pos_diver=0.0;
		for ( i=0;i<pop_size;i++ )
		{
			m_alg_stat.pos_diver += vec_distance(pop[i].x,m_alg_stat.cen_ind);
		}// for every particle
		m_alg_stat.pos_diver /=  pop_size;
		m_alg_stat.pos_diver /= (m_diag_len*1.0);
	}// end function update_diversity

	void calc_gen_avg_bst()
	{
		// find MAX gen
		int max_run=m_ppara->get_max_run();
		int max_gen;
		int i;
		int cur_gen_num;
		max_gen=0;
		for ( i=0;i<max_run;i++ )
		{
			cur_gen_num=m_alg_stat.gen_bst_val[i].size();
			if ( cur_gen_num > max_gen )
				max_gen=cur_gen_num;
		}
		// fill the blank
		int cur_max_gen;
		int j;
		for ( i=0;i<max_run;i++ )
		{
			cur_max_gen=m_alg_stat.gen_bst_val[i].size();
			if ( cur_max_gen < max_gen )
				m_alg_stat.gen_bst_val[i].resize(max_gen,m_alg_stat.run_gbest_val[i]);
		}
		// calculate gen avg best
		m_alg_stat.gen_avg_bst_val.assign(max_gen,0.0);
		for ( j=0;j<max_gen;j++ ) // for every generation
			for ( i=0;i<max_run;i++ ) // for every run
				m_alg_stat.gen_avg_bst_val[j] += m_alg_stat.gen_bst_val[i][j];
		m_alg_stat.gen_avg_bst_val /= max_run;
	}// end function calc_gen_avg_bst

	virtual void calc_gen_avg_vals()
	{
		// find MIN gen
		int max_run=m_ppara->get_max_run();
		int min_gen;
		int i;
		int cur_gen_num;
		min_gen=m_alg_stat.gen_div_val[0].size();
		for ( i=1;i<max_run;i++ )
		{
			cur_gen_num=m_alg_stat.gen_div_val[i].size();
			if ( cur_gen_num < min_gen )
				min_gen=cur_gen_num;
		}

		m_alg_stat.gen_avg_div_val.assign(min_gen,0.0);
		m_alg_stat.gen_avg_rad_val.assign(min_gen,0.0);
		int j;
		for ( j=0;j<min_gen;j++ ) // for every generation
			for ( i=0;i<max_run;i++ ) // for every run
			{
				m_alg_stat.gen_avg_div_val[j] += m_alg_stat.gen_div_val[i][j];// calculate gen avg diversity
				m_alg_stat.gen_avg_rad_val[j] += m_alg_stat.gen_rad_val[i][j];// calculate gen avg radius
			}
			m_alg_stat.gen_avg_div_val /= max_run;
			m_alg_stat.gen_avg_rad_val /= max_run;
	}// end function calc_gen_avg_div

	void calc_avg_gen()
	{
		int max_run=m_ppara->get_max_run();
		int i;
		m_alg_stat.run_avg_gen=0.0;
		for ( i=0;i<max_run;i++ )
			m_alg_stat.run_avg_gen += m_alg_stat.run_gen_val[i];
		m_alg_stat.run_avg_gen /= max_run;
	}

	void calc_avg_stag_gen()
	{
		int max_run=m_ppara->get_max_run();
		m_alg_stat.run_avg_stag_gen /= max_run;
	}

	void calc_avg_conv_eval_num()
	{
		if ( m_alg_stat.conv_count!=0 )
			m_alg_stat.avg_conv_eval_num /= m_alg_stat.conv_count;
		else
			m_alg_stat.avg_conv_eval_num = 0.0;
	}

	void stat_run(population &pop,int cur_run)
	{
		m_alg_stat.run_gbest_val.push_back(m_alg_stat.gbest.obj[0]);
		m_alg_stat.run_gen_val.push_back(m_cur_gen);
		m_alg_stat.run_avg_stag_gen+=m_alg_stat.num_stag;
		int max_run=m_ppara->get_max_run();
		if ( 0==cur_run )// first run
		{
			m_alg_stat.bst_ind=m_alg_stat.gbest;
			m_alg_stat.wst_ind=m_alg_stat.gbest;
			m_alg_stat.run_avg_bst=m_alg_stat.gbest.obj[0];
			m_alg_stat.run_bst_std=0.0;
			// convergence ratio
			m_alg_stat.conv_ratio = m_alg_stat.conv_count/(max_run*1.0);
			if ( 1==max_run ) // run once
			{
				// calculate generational values
				calc_avg_stag_gen();
				calc_gen_avg_bst();
				calc_gen_avg_vals();
				calc_avg_gen();
				calc_avg_conv_eval_num();
			}
			return;
		}
		// update best individual of all trials
		if ( m_alg_stat.gbest < m_alg_stat.bst_ind )
			m_alg_stat.bst_ind=m_alg_stat.gbest;
		// update worst individual of gbest of all trials
		if ( m_alg_stat.gbest > m_alg_stat.wst_ind )
			m_alg_stat.wst_ind=m_alg_stat.gbest;

		m_alg_stat.run_avg_bst += m_alg_stat.gbest.obj[0];

		if ( cur_run==(max_run-1) )
		{
			// average fitness value of gbest
			m_alg_stat.run_avg_bst /= (max_run*1.0);
			// calculate std of gbest
			int i;
			m_alg_stat.run_bst_std=0.0;
			double diff;
			for ( i=0;i<max_run;i++ )
			{
				diff=m_alg_stat.run_gbest_val[i]-m_alg_stat.run_avg_bst;
				m_alg_stat.run_bst_std += diff*diff;
			}
			m_alg_stat.run_bst_std /= (max_run*1.0);
			m_alg_stat.run_bst_std = sqrt(m_alg_stat.run_bst_std);
			// convergence ratio
			m_alg_stat.conv_ratio = m_alg_stat.conv_count/(max_run*1.0);
			// calculate generational average values of all runs
			calc_avg_stag_gen();
			calc_gen_avg_bst();
			calc_gen_avg_vals();
			calc_avg_gen();
			calc_avg_conv_eval_num();
		}
	}// end function stat_run

	virtual int load_ini_pop(population &pop,std::string ini_pop_path)
	{
		int pop_size=m_ppara->get_pop_size();
		std::ifstream ini_pop_data(ini_pop_path.c_str());
		int num_dims=m_ppara->get_dim();
		int i=0;// entry read number
		int j;
		std::string inbuf;
		d_array x_tmp(num_dims);
		bool inconsist_flag=false;
		if ( ini_pop_data.fail() )
		{
			std::stringstream str_err;
			str_err<<"\nopen file "
				<<ini_pop_path
				<<" ERROR!\n"
				<<"please check file path!\n";
			throw std::logic_error(str_err.str());
		}
		while ( !ini_pop_data.eof() )
		{
			std::getline(ini_pop_data,inbuf);
			if ( false==is_dataLine(inbuf) )
				continue;
			if ( i==pop_size ) 
				inconsist_flag=true;
			if ( !inconsist_flag )
			{
				std::stringstream line(inbuf);
				for ( j=0;j<num_dims;j++ )
				{
					line>>x_tmp[j];
					bound_check(x_tmp[j],j);

					pop[i].x[j]=x_tmp[j];
				}
				if ( ini_pop_data.fail() ) // read error
				{
					std::stringstream str_err("ERROR:Read data from file");
					str_err<<" "
						<<ini_pop_path.c_str()
						<<" "
						<<"failed."
						<<"Check your data file for data validity."
						<<"\n";
					throw std::logic_error(str_err.str());
				}
			}// if consistent
			i++;
		}// while !eof
		if ( inconsist_flag ) // entry count is inconsistent with pop_size in config file,data file is invalid.
		{
			std::stringstream str_err("ERROR:Invalid Initialization File data!");
			str_err<<"Initial population file "
				<<ini_pop_path.c_str()
				<<"'s"
				<<" population size="
				<<i
				<<","
				<<"while pop_size in algorithm configuration file is"
				<<" "
				<<pop_size
				<<".";
			throw std::logic_error(str_err.str());
		}
		return 0;
	}// end function load_orig_pop

	void set_orig_pop(population &pop)
	{
		int ini_type=m_ppara->get_ini_type();
		std::string ini_pop_path=m_ppara->get_ini_pop_file_path();
		if ( initype_file==ini_type )
		{
			if ( ini_pop_path.empty() )
			{
				std::stringstream str_err("ERROR:ini_pop_file entry in config file is empty!");
				throw std::logic_error(str_err.str());
			}
			load_ini_pop(pop,ini_pop_path);
		}
		else
		{
			val_range& ini_range=m_ppara->get_ini_rng();
			initializer initializer;// original population initializer
			// initialize original pop
			if ( initype_ortho==ini_type )
				initializer.ini_orthogonal(pop,m_ppara->get_q_number(),ini_range);
			else
				initializer.ini_real_uni(pop,ini_range);
		}

		// evaluate original pop and update x_tmp range
		eval_ini_pop(pop,*m_pfunc,m_alg_stat);

		std::ifstream read_test_file;
		read_test_file.open(ini_pop_path.c_str(), std::ifstream::in);
		read_test_file.close();
		// write to file only if ini_pop file dos not exist
		if( read_test_file.fail() )
		{
			std::ofstream ini_pop_file(ini_pop_path);
			write_initial_pop(ini_pop_file,pop);
			read_test_file.clear(std::ios::failbit);// file already exist,reset ios flag
		}

		stat_ini_pop(pop,m_alg_stat);
	}// end function set_orig_pop

	void update_search_radius()
	{
		int num_dims=m_ppara->get_dim();
		int pop_size=m_ppara->get_pop_size();
		int i,j;
		m_alg_stat.radius.assign(pop_size,0.0);
		m_alg_stat.avg_radius=0.0;
		for ( i=0;i<pop_size;i++ )
		{
			for ( j=0;j<num_dims;j++ )
			{
				double diff=m_alg_stat.delta_x[i][j];
				m_alg_stat.radius[i] += diff*diff;
			}

			m_alg_stat.radius[i]=sqrt(m_alg_stat.radius[i]);
			m_alg_stat.avg_radius += m_alg_stat.radius[i];
		}// for every particle
		m_alg_stat.avg_radius /= pop_size;
	}// end function update_search_radius

	// default boundaries check
	bool bound_check(double &x,int dim)
	{
		const val_range& val_bounds=m_ppara->get_val_bnd();
		bool high_bnd,low_bnd;
		low_bnd=x < val_bounds[dim].min_val;
		if ( low_bnd )
		{
			x=val_bounds[dim].min_val;
			m_alg_stat.all_ob_num++;
			return true;
		}
		high_bnd=x > val_bounds[dim].max_val;
		if ( high_bnd )
		{
			x=val_bounds[dim].max_val;
			m_alg_stat.all_ob_num++;
			return true;
		}
		return false;
	}

	// record average gbest of all runs in every generation 
	void write_avg_gbest_per_gen()
	{
		std::ofstream avg_bst_file(m_com_out_path.avg_bst_path);
		print_gen_avg_val(avg_bst_file,m_alg_stat.gen_avg_bst_val);
	}

	// record average diversity of all runs in every generation
	void write_avg_diver_per_gen()
	{
		std::ofstream avg_div_file(m_com_out_path.avg_div_path);
		print_gen_avg_val(avg_div_file,m_alg_stat.gen_avg_div_val);
	}

	// record average search radius of all runs in every generation
	void write_avg_rad_per_gen()
	{
		std::ofstream avg_rad_file(m_com_out_path.avg_rad_path);
		print_gen_avg_val(avg_rad_file,m_alg_stat.gen_avg_rad_val);
	}

	void write_all_bst_val()
	{
		std::ofstream all_bst_val_file(m_com_out_path.all_bst_val_path);
		print_all_gbest_val(all_bst_val_file,m_alg_stat.run_gbest_val);
	}

	// record stat values of run
	void write_stat_vals()
	{
		write_all_bst_val();
		write_avg_gbest_per_gen();
		write_avg_diver_per_gen();
		write_avg_rad_per_gen();
	}

	void reset_run_stat()
	{
		m_alg_stat.reset_run_stat();
		m_convergent=false;
		m_alg_stat.avg_radius=0.0;
	}

	void update_conv_stat(double vtr)
	{
		if ( !m_convergent && m_pfunc->is_convergent( m_alg_stat.gbest.obj[0],vtr) )// reach global optimum BEFORE stop condition is met
		{
			m_alg_stat.conv_count++;
			m_alg_stat.avg_conv_eval_num += m_alg_stat.eval_num;
			m_convergent=true;
		}
	}

	std::string m_out_prefix;
	std::string m_para_val_desc;
	std::string m_outfile_common;
	std::string m_outfile_para;
	common_path m_com_out_path;
	double m_diag_len; // length of longest diagonal of search space

	int m_cur_run;// run times counter
	int m_cur_gen;// generation counter

	boost::shared_ptr<real_para> m_ppara;
	alg_stat m_alg_stat;// algorithm statistics

	class stop_cond // virtual base class for stop criterion
	{
	public:
		stop_cond(com_alg *pself)
			:self(pself) {}
		virtual operator bool()=0;
	protected:
		com_alg *self;// pimp idiom:self pointer to container class for nested class
	};

	class stop_gen:public stop_cond
	{
	public:
		stop_gen(com_alg *pself)
			:stop_cond(pself) {}
		operator bool() { return self->m_cur_gen >= (self->m_ppara->get_max_gen()); }
	};

	class stop_eval:public stop_cond
	{
	public:
		stop_eval(com_alg *pself)
			:stop_cond(pself) {}
		operator bool() { return self->m_alg_stat.eval_num >= (self->m_ppara->get_stop_val()); }
	};

	class stop_stag:public stop_cond
	{
	public:
		stop_stag(com_alg *pself)
			:stop_cond(pself) {}
		operator bool() { return self->m_alg_stat.num_stag >= (self->m_ppara->get_stop_val()); }
	};
	boost::shared_ptr<stop_cond> m_pstop_cond;

	void alloc_stop_cond()
	{
		int stop_type=m_ppara->get_stop_type();
		switch (stop_type)
		{
		case stoptype_eval:
			m_pstop_cond=shared_ptr<stop_cond>(new stop_eval(this));
			break;
		case stoptype_stag:
			m_pstop_cond=shared_ptr<stop_cond>(new stop_stag(this));
			break;
		default:
			m_pstop_cond=shared_ptr<stop_cond>(new stop_gen(this));
		};
	}

	d_array radius;// search radius of every individual,stat variable
	double avg_radius;// mean search radius of population

	boost::shared_ptr<func_base> m_pfunc;
	void set_eval_func(boost::shared_ptr<func_base> &pfunc_Eval) { m_pfunc=pfunc_Eval; }

	bool m_convergent;

	void print_avg_gen(std::ostream &os,double avg_gen)
	{
		int stoptype=m_ppara->get_stop_type();
		if ( stoptype_stag==stoptype || stoptype_delta==stoptype )
			os<<"\n"
			<<"average_gen_till_stag:"
			<<avg_gen;
	}

	bool is_output_gen(int gen_idx,int out_interval) {return 0==((gen_idx+1)%out_interval);}

	bool is_learn_gen(int cur_gen,int learn_p) { return ( 0==(cur_gen+1)%learn_p ); }

	bool is_final_run(int cur_run,int max_run) { return (max_run-1)==cur_run;}

	//// OBSELETE function
	//void alloc_prog_indicator(boost::shared_ptr<boost::progress_display> &ppd)
	//{
	//	int max_gen=m_ppara->get_max_gen();
	//	int max_run=m_ppara->get_max_run();// run/trial number
	//	if ( 1==max_run )
	//		ppd=shared_ptr<progress_display>(new progress_display(max_gen-1));
	//	else
	//		ppd=shared_ptr<progress_display>(new progress_display(max_run));
	//}
};// end class com_alg declaration

#endif
