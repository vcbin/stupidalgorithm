#ifndef STUPIDALGO_SP_DIFFERENTIAL_EVOLUTION_ALGORITHM 
#define STUPIDALGO_SP_DIFFERENTIAL_EVOLUTION_ALGORITHM

#include "BBDE_Alg.h"

namespace de
{
	namespace spde
	{
		class spde_alg:public bbde::bbde_alg
		{
		public:
			// ctor
			spde_alg(std::string conf_path):
			  bbde_alg(conf_path) {}
			  // set algorithm parameters pointer directly
			  spde_alg(boost::shared_ptr<de_para> ppara): 
			  bbde_alg(ppara) {}

			  int run();
		protected:
			float m_ini_f;// initial F value
		};// end class spde_alg declaration

	}// end namespace spde
}// end namespace de
#endif
