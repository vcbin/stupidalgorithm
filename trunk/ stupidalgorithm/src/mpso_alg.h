#ifndef STUPIDALGO_MPSO_ALGORITHM
#define STUPIDALGO_MPSO_ALGGORITHM

#include "pso_alg.h"

namespace pso
{
	namespace mpso
	{
		class mpso_alg:public pso_alg
		{
		public:
			// ctor
			mpso_alg(std::string conf_path):
			  pso_alg(conf_path) {}

			  void initialize();
			  int run();

			  void print_run_title(std::ostream &os);
			  void print_gen_stat(std::ostream &os,int cur_gen,const alg_stat &alg_stat);
			  // functions for batch run overall stat file writing
			  inline void  set_gen_avg_attra_file_name( std::string str_filename ){ m_avg_attr_path=str_filename; }
		private:
			void print_rep_stat(std::ostream &os);
			void print_rep_run_stat(std::ostream &os);
			void update_speed(population &pop);
			inline void set_diver_bound() { m_d_l=m_dl_cof*m_diag_len;m_d_h=m_dh_cof*m_diag_len; }
			void update_accle(population &pop,int cur_gen);

			inline void update_rep_size() { m_mpso_stat.rep_size++; }
			void update_rep_stat(int cur_gen);
			void calc_run_rep_stat();
			void stat_rep_run();
			void calc_run_avg_rep_stat();
			void calc_gen_avg_vals();

			std::string m_avg_attr_path;
			void write_avg_attr_per_gen();
			void write_stat_vals();

			d_array accle;// acceleration a1 of every dimenson for every particle
			double m_d_l;// attractive distance in genotype
			double m_d_h;// repulsive distance in genotype
			double m_dl_cof;// d_l coefficient for attractive distance calculation
			double m_dh_cof;// d_h coefficient for repulsive distance calculation

			// stat variables
			struct Mpso_stat
			{
				Mpso_stat():
			rep_size(0),
				rep_perc(0.0),
				rep_gen_num(0),
				rep_gen_perc(0.0),
				avg_rep_size(0.0),
				avg_rep_perc(0.0),
				first_rep_gen(-1),
				first_rep_flag(true),
				rep_run(0),
				avg_rep_gen_num(0.0),
				avg_rep_gen_perc(0.0),
				avg_avg_rep_size(0.0),
				avg_avg_rep_size_perc(0.0) {}
			inline void alloc(int max_run) { gen_attr_val.resize(max_run); }
			int rep_size;
			double rep_perc;
			int rep_gen_num;
			double rep_gen_perc;
			double avg_rep_size;
			double avg_rep_perc;
			int rep_gen;
			int first_rep_gen;
			bool first_rep_flag;// auxiliary variable
			int rep_run;// repulsion run/trial number
			double avg_rep_gen_num;
			double avg_rep_gen_perc;
			double avg_avg_rep_size;
			double avg_avg_rep_size_perc;
			d_mat gen_attr_val;
			d_array gen_avg_attr_val;
			inline void reset_gen_stat() {rep_size=0;rep_perc=0.0;}
			inline void reset_run_stat() {rep_gen_num=0;rep_gen_perc=0.0;first_rep_gen=-1;first_rep_flag=true;avg_rep_size=0.0;avg_rep_perc=0.0;}
			} m_mpso_stat;

			void record_attr_perc(Mpso_stat &mpso_stat,int cur_run);
			void record_gen_vals(alg_stat &alg_stat,int cur_run);
		};// end calss mpso_alg
	}// namespace mpso

}// namespace pso


#endif
