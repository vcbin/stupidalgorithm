#ifndef JERRYLIU_J_SELF_ADAPTIVE_DIFFERENTIAL_EVOLUTION_ALGORITHM
#define JERRYLIU_J_SELF_ADAPTIVE_DIFFERENTIAL_EVOLUTION_ALGORITHM

#include "BBDE_Alg.h"

namespace de
{
	namespace jde
	{
		class jde_alg:public bbde::bbde_alg
		{
		public:
			// ctor
			jde_alg(std::string conf_path):
					bbde_alg(conf_path) {}
			 // set algorithm parameters pointer directly
			 jde_alg(boost::shared_ptr<de_para> ppara): 
					bbde_alg(ppara) {}

			  int run();
		};// end class jde_alg declaration

	}// end namespace jde
}// end namespace de


#endif