#ifndef STUPIDALGO_EDA_DIRECTED_DIFFERENTIAL_EVOLUTION_ALGORITHM
#define STUPIDALGO_EDA_DIRECTED_DIFFERENTIAL_EVOLUTION_ALGORITHM

#include "de-eda_alg.h"

namespace de
{
	namespace dmde
	{
		class dmde_alg:public de_eda::de_eda_alg
		{
		public:
			dmde_alg(std::string conf_path):
			  de_eda_alg(conf_path) {}
		protected:
			int run();
		};// end class ede_alg

	}// namespace dmde
}// namespace de


#endif