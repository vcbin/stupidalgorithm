#ifndef STUPIDALGO_PSOBC_ALGORITHM
#define STUPIDALGO_PSOBC_ALGORITHM

#include "arpso_alg.h"

namespace pso
{
	namespace psobc
	{
		struct psobc_stat
		{
			std::vector<individual> pworst;
			individual gworst;
			d_array pre_pop_val;
			void alloc_spaces(size_t pop_size,int num_dims);
		};

		bool operator>(individual &lhs,double val);
		int copy_obj_val_from_pop(d_array &pop_val,const population &pop);

		class psobc_alg:public arpso::arpso_alg
		{
		public:
			// ctor
			psobc_alg(std::string conf_path):
			  arpso_alg(conf_path) {}

			  void initialize();
			  int run();

			  // overridden function to adapt to psobc algo specific procedure
			  void stat_ini_pop( population &pop,
								 alg_stat &alg_stat
								);
			  //void eval_pop(population &pop,func_base &func_eval,bool first_time);
			  void stat_pop(population &pop,alg_stat &alg_stat);
			  
			  void update_speed(population &pop);
		private:
			psobc_stat m_psobc_stat;
			std::vector<std::vector<bool> > m_ob_handled;
			d_mat m_prev_dir;// for O.B condition handling auxiliary variable
		};

	}// namespace psobc
}// namespace pso

#endif