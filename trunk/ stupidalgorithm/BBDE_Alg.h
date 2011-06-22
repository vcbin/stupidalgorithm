#ifndef STUPIDALGO_DIFFERENTIAL_EVOLUTION_ALGORITHM
#define STUPIDALGO_DIFFERENTIAL_EVOLUTION_ALGORITHM

#include <boost/progress.hpp>
#include <boost/timer.hpp>
#include <boost/shared_ptr.hpp>

#include "DE_Alg.h"
#include "Algo_DataStruct.h"

namespace de
{
	namespace bbde
	{
		class bbde_alg:public de_alg
		{
		public:
			// ctor
			bbde_alg(std::string conf_path):
			  de_alg(conf_path) {}
			  bbde_alg(boost::shared_ptr<de_para> ppara):// set algorithm parameters pointer directly
			  de_alg(ppara) {}
			  // dtor
			  virtual ~bbde_alg() { }

			  int run();
		protected:
			double generate_rnd_pr(double mean=0.0,double sigma=1.0);
		};// end class bbde_alg declaration
	}// end namespace bbde
}// end namespace de

#endif
