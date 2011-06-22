#ifndef STUPIDALGO_BASIC_DIFFERENTIAL_EVOLUTION_ALGORITHM
#define STUPIDALGO_BASIC_DIFFERENTIAL_EVOLUTION_ALGORITHM

#include <boost/progress.hpp>
#include <boost/timer.hpp>
#include <boost/shared_ptr.hpp>

#include "com_alg.h"
#include "de_para.h"
#include "Algo_DataStruct.h"

namespace de
{
	class de_alg:public com_alg<de_para>
	{
	public:
		// ctor
		de_alg(std::string conf_path):
		  com_alg(conf_path) {}
		  de_alg(boost::shared_ptr<de_para> ppara):// set algorithm parameters pointer directly
		  com_alg(ppara) {}
		  // dtor
		  virtual ~de_alg() { }

		  void initialize();
		  int run();
		  inline void  set_gen_avg_f_file_name( std::string str_filename ){ m_avg_f_path=str_filename; }
		  inline void  set_gen_avg_pr_file_name( std::string str_filename ){ m_avg_pr_path=str_filename; }
	protected:
		static const int stra_num=2;
		enum stra_name{pr,f};

		void update_pop(population &pop,const population &trial_pop);

		void bound_check(double &x,int dim);// SPECIAL boundaries check

		void calc_de_para_stat(const population& pop);
		struct de_stat
		{
			de_stat():
				f_mean(0.0),f_std(0.0),
				pr_mean(0.0),pr_std(0.0) {}

		void initialize(int max_run);
		void reset() { f_mean=f_std=pr_mean=pr_std=0.0; }
		double f_mean;
		double f_std;
		double pr_mean;
		double pr_std;
		d_mat gen_f;
		d_array gen_avg_f;
		d_mat gen_pr;
		d_array gen_avg_pr;
		void alloc(int max_run);
		};
		de_stat m_de_stat;
		//d_mat m_f;// 2 dimensional array for extensibility,maximum strategy parameters granularity: DIMENSION
		//d_mat m_pr;
		//d_mat m_trial_f;// array to record generated trial value of strategy parameters
		//d_mat m_trial_pr;
		void calc_gen_avg_de_para();
		void record_avg_f(int cur_run);
		void record_avg_pr(int cur_run);
		void record_de_para_stat(int cur_run);
		void write_avg_de_para_per_gen();
		void write_stat_vals();

		std::string m_avg_f_path;
		std::string m_avg_pr_path;
		void stat_run(population &pop,int cur_run);

		void print_run_title(std::ostream &os);
		void print_gen_stat(std::ostream &os,int cur_gen,const alg_stat &alg_stat);
	};// end class de_alg declaration
}// end namespace de

#endif
