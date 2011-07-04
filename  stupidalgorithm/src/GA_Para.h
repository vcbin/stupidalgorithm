#ifndef DGEA_ALL_PARAMETERS
#define DGEA_ALL_PARAMETERS

#include <fstream>
#include <iostream>

#include "Para_Base.h"

namespace ga
{
	class ga_para:public para_base
	{
	public:
		// ctor
		ga_para(std::string conf_path="")
			:para_base(conf_path)
		{
			general.add_options()
				("pr", po::value<std::string>(), "probability of performing crossover")
				("pm", po::value<std::string>(), "probability of mutation")
				("d_l",po::value<std::string>(&d_l)->default_value("5.0E-6"),"default value for lower bound of attration")
				("d_h",po::value<std::string>(&d_h)->default_value("0.25"),"default value for upper bound of attration");
		}
		virtual ~ga_para(){}
		int read_para(std::string conf_path="");

		double get_pr() {return get_seq_val<double>("pr");}
		double get_pm() {return get_seq_val<double>("pm");}
		double get_dl() {return get_seq_val<double>("d_l");}
		double get_dh() {return get_seq_val<double>("d_h");}
	private:
		std::string pr;	// probability of crossover
		std::string pm;	// probability of mutation
		std::string d_l;// OPTIONAL,default value for lower bound of attration
		std::string d_h;// OPTIONAL,default value for upper bound of attration
	};// end class ga_para
}//end namespace ga

#endif
