#ifndef STUPIDALGO_NORMAL_SELFADAPTIVE_DIFFERENTIAL_EVOLUTION
#define STUPIDALGO_NORMAL_SELFADAPTIVE_DIFFERENTIAL_EVOLUTION

#include "bbde_alg.h"

namespace de
{
	namespace nsde
	{
		class nsde_alg:public bbde::bbde_alg
		{
		public:
			nsde_alg(std::string conf_path):
			  bbde_alg(conf_path) {}
			  inline void set_gen_avg_bi_norm_stat_file_name(std::string name){m_gen_avg_bi_norm_stat_path=name;}
			  void initialize();
			  int run();
		protected:
			std::string m_gen_avg_bi_norm_stat_path;
			bi_norm_var m_bi_norm_var;
			struct prob_stat
			{
				int gen;
				bi_norm_var var;
				prob_stat() {}
				prob_stat(int r_gen,const bi_norm_var& r_var)
					:gen(r_gen),var(r_var){}
			};
			inline void record_gen_bi_norm_stat(bi_norm_var &norm_stat,int cur_gen,int cur_run) { m_gen_bi_norm_stat[cur_run].push_back(prob_stat(cur_gen,norm_stat)); }
			void calc_gen_avg_bi_norm_stat();
			void stat_run(population &pop,int cur_run);
			void write_stat_vals();
			//void bound_check(double &tri_x,int dim);// 
			//void record_gen_vals(alg_stat &alg_stat,int cur_run);

			void update_pop(population &pop,const population &trial_pop);

			enum ini_f_type{uni=1,norm=2};

			std::vector<prob_stat> m_avg_bi_norm_var;// mean stat over all run
			d_array m_succ_f;
			d_array m_succ_cr;
			std::vector<std::vector<prob_stat> > m_gen_bi_norm_stat;
		};//end class nsde_alg

	}// end namespace nsde
}// end namespace de

#endif