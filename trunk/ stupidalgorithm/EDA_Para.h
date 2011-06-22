#ifndef JERRYLIU_EDA_ALL_PARAMETERS
#define JERRYLIU_EDA_ALL_PARAMETERS

#include <fstream>
#include <iostream>

#include "para_base.h"

namespace eda
{
	class eda_para:public para_base
	{
	public:
		// ctor
		eda_para(std::string conf_path="")
			:para_base(conf_path) 
		{
			general.add_options()
				("m_ratio",po::value<std::string>(&m_ratio)->default_value("0.5"),"ratio of parents population");
		}
		int read_para(std::string conf_path="");

		double get_m_ratio(){return get_seq_val<double>("m_ratio");}
	private:
		std::string m_ratio;// ratio of parents population
	};// end class eda_para

}// namespace eda
#endif