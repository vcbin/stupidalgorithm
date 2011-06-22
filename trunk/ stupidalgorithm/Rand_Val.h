#ifndef JERRYLIU_BOOST_RAND_NUMBER_GENERATOR_HEADER
#define JERRYLIU_BOOST_RAND_NUMBER_GENERATOR_HEADER

// boost random floating-point uniform [0,1) distribution header
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/cauchy_distribution.hpp>
#include <boost/random/variate_generator.hpp>

double gen_rnd_norm(double mean,double sigma);

#endif