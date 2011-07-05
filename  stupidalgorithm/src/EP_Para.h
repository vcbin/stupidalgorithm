#ifndef STUPIDALGO_FEP_ALL_PARAMETERS
#define STUPIDALGO_FEP_ALL_PARAMETERS

#include <fstream>
#include <iostream>

#include "para_base.h"

namespace ep
{
	class ep_para:public para_base
	{
	public:
		// ctor
		ep_para(std::string conf_path="")
			:para_base(conf_path) 
		{
			general.add_options()
				("tour_size",po::value<std::string>(),"tournament size of tournament selection")
				("mut_type",po::value<std::string>(&mut_type)->default_value("1"),"random number type of mean step size randomization")
				("ini_eta",po::value<std::string>(&ini_eta)->default_value("3.0"),"initial value of scaling factor eta");
		}
		virtual ~ep_para(){}

		int read_para(std::string conf_path="");

		int get_tour_size() {return get_seq_val<int>("tour_size");}
		int get_mut_type() {return get_seq_val<int>("mut_type");}
		double get_ini_eta() {return get_seq_val<double>("ini_eta");}
	private:
		std::string tour_size;// tournament size of tournament selection
		std::string mut_type;// random number type of mean step size randomization
		std::string ini_eta;// 
	};// end class declaration
}// namespace ep
#endif
