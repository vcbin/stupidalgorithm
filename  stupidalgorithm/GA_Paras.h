#ifndef DGEA_ALL_PARAMETERS
#define DGEA_ALL_PARAMETERS

#include <fstream>
#include <iostream>

#include "Para_Base.h"

namespace ga
{
	class GA_Paras:public Para_Base
	{
	public:
		// ctor
		GA_Paras(std::string conf_path="")
			:Para_Base(conf_path)
		{
			general.add_options()
				("pr", po::value<std::string>(), "probability of performing crossover")
				("pm", po::value<std::string>(), "probability of mutation")
				("d_l",po::value<std::string>(&d_l)->default_value("5.0E-6"),"default value for lower bound of attration")
				("d_h",po::value<std::string>(&d_h)->default_value("0.25"),"default value for upper bound of attration");
		}
		int ReadPara(std::string conf_path="");

		double GetPr() {return GetSeqVal<double>("pr");}
		double GetPm() {return GetSeqVal<double>("pm");}
		double GetD_l() {return GetSeqVal<double>("d_l");}
		double GetD_h() {return GetSeqVal<double>("d_h");}
	private:
		std::string pr;	// probability of crossover
		std::string pm;	// probability of mutation
		std::string d_l;// OPTIONAL,default value for lower bound of attration
		std::string d_h;// OPTIONAL,default value for upper bound of attration
	};// end class GA_Paras
}//end namespace ga

#endif
