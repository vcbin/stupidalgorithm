#ifndef JERRYLIU_EDA_ALL_PARAMETERS
#define JERRYLIU_EDA_ALL_PARAMETERS

#include <fstream>
#include <iostream>

#include "Para_Base.h"

namespace eda
{
	class EDA_Paras:public Para_Base
	{
	public:
		// ctor
		EDA_Paras(std::string conf_path="")
			:Para_Base(conf_path) 
		{
			general.add_options()
				("m_ratio",po::value<std::string>(&m_ratio)->default_value("0.5"),"ratio of parents population");
		}
		int ReadPara(std::string conf_path="");

		double GetMRatio(){return GetSeqVal<double>("m_ratio");}
	private:
		std::string m_ratio;// ratio of parents population
	};// end class EDA_Paras

}// namespace eda
#endif