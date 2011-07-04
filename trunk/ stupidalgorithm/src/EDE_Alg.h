#ifndef STUPIDALGO_ENHANCED_DIFFERENTIAL_EVOLUTION_ALGORITHM
#define STUPIDALGO_ENHANCED_DIFFERENTIAL_EVOLUTION_ALGORITHM

#include "DE_Alg.h"

namespace de
{
	namespace ede
	{
		class ede_alg:public de_alg
		{
		public:
			ede_alg(std::string conf_path):
		  de_alg(conf_path) {}
		protected:
			int run();
		};// end class ede_alg

	}// namespace ede
}// namespace de


#endif