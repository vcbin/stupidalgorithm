#ifndef JERRYLIU_DPSO_MY_ALGORITHM
#define JERRYLIU_DPSO_MY_ALGORITHM

#include "DPSO_Alg.h"

namespace pso
{
	namespace dpso_m
	{
		class dpso_m_alg:public dpso::dpso_alg
		{
		public:
			// ctor
			dpso_m_alg(std::string conf_path):
				dpso_alg(conf_path) {}
		private:
			void update_momentum(population &pop);
		};

	}// namespace dpso_m
}// namespace pso


#endif