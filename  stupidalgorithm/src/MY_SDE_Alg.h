#ifndef STUPIDALGO_MY_SELF_ADAPTIVE_DIFFERENTIAL_EVOLUTION_ALGORITHM
#define STUPIDALGO_MY_SELF_ADAPTIVE_DIFFERENTIAL_EVOLUTION_ALGORITHM

#include "bbde_alg.h"

namespace de
{
	namespace sde
	{
		class my_sde_alg:public bbde::bbde_alg
		{
		public:
			// ctor
			my_sde_alg(std::string conf_path):
			  bbde_alg(conf_path) {}
			  // set algorithm parameters pointer directly
			my_sde_alg(boost::shared_ptr<de_para> ppara): 
			  bbde_alg(ppara) {}
			  void initialize();

			  int run();
		protected:
			enum ini_f_type{uni=1,norm=2};
			enum pr_type{stra_norm=1,stra_learn};
			// overridden
			inline bool is_learn_gen(int cur_gen,int learn_p);
			void update_pop(population &pop,const population &trial_pop);
			float m_ini_f;// initial F value
			d_mat m_succ_f;
			d_array m_succ_pr;
		};// end class sde_alg declaration

	}// end namespace sde
}// end namespace de
#endif