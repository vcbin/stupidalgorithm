#ifndef STUPIDALGO_IMPROVED_FAST_EVOLUTIONARY_PROGRAMMING
#define STUPIDALGO_IMPROVED_FAST_EVOLUTIONARY_PROGRAMMING

#include "com_alg.h"
#include "FEP_Alg.h"
#include "ep_para.h"

namespace ep
{
	namespace ifep
	{
		class ifep_alg:public fep::fep_alg
		{
		public:
			// ctor
			ifep_alg(std::string conf_path):
					fep_alg(conf_path) {}
			ifep_alg(boost::shared_ptr<ep_para> ppara):// set algorithm parameters pointer directly
					fep_alg(ppara) {}

			  int run();
		protected:
			int load_ini_pop(population &pop,std::string ini_pop_path);
			void generate_offspring(population &parent,population &offspring);
		};// end class ifep_algo declaration

	}// end namespace ifep
}// end namespace ep


#endif