#ifndef STUPIDALGO_HYBRID_DIFFERENTIAL_EVOLUTION_ALGORITHM
#define STUPIDALGO_HYBRID_DIFFERENTIAL_EVOLUTION_ALGORITHM

#include "de_alg.h"

namespace de
{
	namespace de_eda
	{
		class de_eda_alg:public de_alg
		{
		public:
			de_eda_alg(std::string conf_path):
			  de_alg(conf_path) {}
		protected:
			void initialize();
			int run();

			void stat_eda_dist(population &pop);
			double gen_eda_x(int dim);
		private:
			d_array m_x_mean;
			d_array m_x_std;
		};// end class de_ede_alg

	}// namespace de_eda
}// namespace de
#endif