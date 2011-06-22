#ifndef JERRYLIU_SELF_ADAPTIVE_DIFFERENTIAL_EVOLUTION_ALGORITHM
#define JERRYLIU_SELF_ADAPTIVE_DIFFERENTIAL_EVOLUTION_ALGORITHM

#include "BBDE_Alg.h"

namespace de
{
	namespace sde
	{
		class sde_alg:public bbde::bbde_alg
		{
		public:
			// ctor
			sde_alg(std::string conf_path):
					bbde_alg(conf_path) {}
			 // set algorithm parameters pointer directly
			 sde_alg(boost::shared_ptr<de_para> ppara): 
					bbde_alg(ppara) {}

			  int run();
		};// end class sde_alg declaration

	}// end namespace sde
}// end namespace de
#endif