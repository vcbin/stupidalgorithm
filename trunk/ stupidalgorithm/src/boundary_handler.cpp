#include "rand_val.h"

#include "boundary_handler.h"

using boost::uniform_01;
using boost::mt19937;
using boost::variate_generator;

extern boost::mt19937 gen;

namespace boundary_condition
{
	void bh_damp::handle(double &velocity) 
	{
		uniform_01<> dist;
		variate_generator<mt19937&, uniform_01<> > rnd_num(gen, dist);
		velocity *= (-1.0*rnd_num());
	}

}// end namespace boundary_condition
