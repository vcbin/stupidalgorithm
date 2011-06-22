#ifndef STUPIDALGO_ARPSO_ALGORITHM
#define STUPIDALGO_ARPSO_ALGORITHM

#include "pso_alg.h"

namespace pso
{
	namespace arpso
	{
		class arpso_alg:public pso_alg
		{
		public:
			//ctor
			arpso_alg(std::string conf_path):
					pso_alg(conf_path),
					m_dir_accl(1.0) {}

			 void initialize();
			 int run();
		protected:
			void print_rep_stat(std::ostream &os);
			void print_rep_run_stat(std::ostream &os);

			struct arpso_stat
			{
				arpso_stat():
					rep_gen_num(0),
					first_rep_gen(-1),
					rep_gen_perc(0.0),
					first_rep_flag(true),
					rep_run(0),
					avg_rep_gen_num(0.0),
					avg_rep_gen_perc(0.0) {}
				int rep_gen_num;// repulsion generation number
				int first_rep_gen;// first repulsion generation
				double rep_gen_perc;
				bool first_rep_flag;// auxiliary variable
				int rep_run;// repulsion run/trial number
				double avg_rep_gen_num;
				double avg_rep_gen_perc;
				void reset() {rep_gen_num=0;first_rep_gen=-1;rep_gen_perc=0.0;first_rep_flag=true;}
			} m_arpso_stat;
			void update_speed(population &pop);
			void update_dir(population &pop);
			void calc_rep_gen_perc();
			void calc_run_avg_rep_stat();
			void stat_rep_gen(int cur_gen);
			void stat_rep_run();
			double m_dir_accl;// direction of speed
			double m_d_l;// attractive distance in genotype
			double m_d_h;// repulsive distance in genotype
		};// end class arpso_alg
	}// namespace arpso
}// namespace pso

#endif