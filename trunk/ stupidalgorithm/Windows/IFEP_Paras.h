#ifndef JERRYLIU_IMPROVED_FEP_ALL_PARAMETERS
#define JERRYLIU_IMPROVED_FEP_ALL_PARAMETERS

#include <fstream>
#include <iostream>

#include "Para_Base.h"

namespace ep
{
	class IFEP_Paras:public Para_Base
	{
	public:
		// ctor
		IFEP_Paras(std::string conf_path="")
			:Para_Base(conf_path) 
		{
			general.add_options()
				("tour_size",po::value<int>(),"tournament size of tournament selection");
		}
		int ReadPara(std::string conf_path="");

		int GetTourSize() {return tour_size;}
	private:
		int tour_size;// tournament size of tournament selection
	};// end class declaration
}// namespace ep
#endif

