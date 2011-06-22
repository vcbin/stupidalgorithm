#include "Rand_Val.h"

using boost::normal_distribution;
using boost::mt19937;
using boost::variate_generator;

extern boost::mt19937 gen;

double gen_rnd_norm(double mean,double sigma) 
{
	normal_distribution<> dist_norm(mean,sigma);
	variate_generator<mt19937&, normal_distribution<> > rnd_norm(gen, dist_norm);
	return rnd_norm(); 
}