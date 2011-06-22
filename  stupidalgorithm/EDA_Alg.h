#ifndef JERRYLIU_BASIC_ESTIMATION_OF_DISTRIBUTION_ALGORITHM
#define JERRYLIU_BASIC_ESTIMATION_OF_DISTRIBUTION_ALGORITHM

#include <boost/progress.hpp>
#include <boost/timer.hpp>
#include <boost/shared_ptr.hpp>

#include "com_alg.h"
#include "eda_para.h"
#include "Algo_DataStruct.h"

namespace eda
{
	class eda_alg:public com_alg<eda_para>
	{
	public:
		// ctor
		eda_alg(std::string conf_path):
		  com_alg(conf_path) {}
		  eda_alg(boost::shared_ptr<eda_para> ppara):// set algorithm parameters pointer directly
		  com_alg(ppara) {}
		  // dtor
		  virtual ~eda_alg() { }

		  void initialize();
		  int run();

	protected:
		void update_pop(population &pop,const population &trial_pop);
		void bound_check(double &x,int dim);// SPECIAL boundaries check

		void stat_eda_dist(population &pop);
		double gen_eda_x(int dim);
		d_array m_x_mean;
		d_array m_x_std;
	};// end class eda_alg declaration

}// end namespace eda

#endif
