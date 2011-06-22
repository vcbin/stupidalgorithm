#ifndef JERRYLIU_FEP_ALL_PARAMETERS
#define JERRYLIU_FEP_ALL_PARAMETERS

#include <fstream>
#include <iostream>

#include "Para_Base.h"

namespace ep
{
	class EP_Paras:public Para_Base
	{
	public:
		// ctor
		EP_Paras(std::string conf_path="")
			:Para_Base(conf_path) 
		{
			general.add_options()
				("tour_size",po::value<std::string>(),"tournament size of tournament selection")
				("mut_type",po::value<std::string>(&mut_type)->default_value("1"),"random number type of mean step size randomization")
				("ini_eta",po::value<std::string>(&ini_eta)->default_value("3.0"),"initial value of scaling factor eta");
		}
		int ReadPara(std::string conf_path="");

		int GetTourSize() {return GetSeqVal<int>("tour_size");}
		int GetMutType() {return GetSeqVal<int>("mut_type");}
		double GetIniEta() {return GetSeqVal<double>("ini_eta");}
	private:
		std::string tour_size;// tournament size of tournament selection
		std::string mut_type;// random number type of mean step size randomization
		std::string ini_eta;// 
	};// end class declaration
}// namespace ep
#endif
