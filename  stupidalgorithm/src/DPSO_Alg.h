#ifndef STUPIDALGO_DPSO_ALGORITHM
#define STUPIDALGO_DPSO_ALGORITHM

#include "PSO_Alg.h"
#include "Algo_DataStruct.h"
#include "Output_Base.h"

namespace pso
{
	namespace dpso
	{

		static const unsigned pop_count=2;

		struct dpso_stat
		{
			dpso_stat():all_exch_num(0),delta_avg(0.0) {}
			std::vector<individual> gbest_inds;
			d_mat cen_inds;
			std::vector<int> gbest_idx;
			// d_array divers;
			d_array stds;
			d_array avgs;
			d_array pos_divers;
			d_mat v_dim_divers;
			d_array vel_divers;
			int all_exch_num;
			double delta_avg;// difference of two subpopulation's average fitness
			void alloc_space(size_t num_dims);
			void reset();
		};

		typedef std::vector<population> SubPops;

		class dpso_alg:public pso_alg
		{
		public:
			// ctor
			dpso_alg(std::string conf_path):
				pso_alg(conf_path) {}

			  void initialize();
			  int run();
			  // overridden output function
			  void print_gen_stat(std::ostream &os,int cur_gen,const alg_stat &alg_stat);
			  void print_run_title(std::ostream &os);
		protected:

			// dpso algo specific
			void update_speed(population &pop);
			void update_diversity(population &pop);
			virtual void update_momentum(population &pop);
			void update_temp();
			void update_pool(population &pop);
			void diffuse(population &pop,func_base &func_eval);

			// overridden function to handle additional algo stat
			void stat_ini_pop(population &pop,alg_stat &alg_stat);
			void stat_pop(population &pop,alg_stat &alg_stat);

			void print_Exch_stat(std::ostream &os);

			d_array m_moment;// momentum of all particle
			double m_t;// temperature of all particle
			std::vector<std::vector<individual> > m_diff_pool;
			dpso_stat m_dpso_stat;
		};
	}// namespace dpso
}// namespace pso



#endif
