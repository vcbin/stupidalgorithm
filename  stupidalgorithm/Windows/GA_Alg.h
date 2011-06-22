#ifndef DIVERSITY_GUIDED_EVOLUTIONARY_ALGORITHM
#define DIVERSITY_GUIDED_EVOLUTIONARY_ALGORITHM

#include <boost/progress.hpp>
#include <boost/timer.hpp>
#include <boost/shared_ptr.hpp>

#include "com_alg.h"
#include "ga_para.h"
#include "de_para.h"
#include "Algo_DataStruct.h"

namespace ga
{
	namespace dgea
	{
		class dgea_alg:public com_alg<ga_para>
		{
		public:
			// ctor
			dgea_alg(std::string conf_path)
				:com_alg(conf_path) {}
			dgea_alg(boost::shared_ptr<ga_para> ppara)
				:com_alg(ppara) {}
			// dtor
			virtual ~dgea_alg() {}

			int run();
		protected:
			void crossover(const individual& ind1,const individual& ind2, individual& ind3, individual& ind4);
			void update_mode(int& mode,double d_l,double d_h);
			void gen_child(int mode,const population& pop,population& child_pop);
			void tour_breed(const population& pop, population& new_pop);
			void mutation(population& pop);
			void mutation_ind(individual& ind);
			void select(population& pop, population& child_pop);
			const individual& tournament(const individual& ind1,const individual& ind2);
			enum ga_mode{explore=1, exploit};

		};// end class dgea_alg
	}// end namespace dgea
}// end namespace ga

#endif
