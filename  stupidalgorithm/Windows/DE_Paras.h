#ifndef JERRYLIU_DE_ALL_PARAMETERS
#define JERRYLIU_DE_ALL_PARAMETERS

#include <fstream>
#include <iostream>

#include "Para_Base.h"

namespace de
{
	class DE_Paras:public Para_Base
	{
	public:
		// ctor
		DE_Paras(std::string conf_path="")
			:Para_Base(conf_path) 
		{
			general.add_options()
				("pr_val",po::value<std::string>(&pr_val)->default_value("0.9"),"crossover probability fixed value")
				("f_val",po::value<std::string>(&f_val)->default_value("0.5"),"amplification factor F fixed value")
				("pr_mean",po::value<std::string>(&pr_mean)->default_value("0.5"),"mean of crossover probability normal distribution")
				("pr_sigma",po::value<std::string>(&pr_sigma)->default_value("0.15"),"sigma of crossover probability normal distribution")
				("f_mean",po::value<std::string>(&f_mean)->default_value("0.7"),"mean of control parameter F normal distribution")
				("f_sigma",po::value<std::string>(&f_sigma)->default_value("0.3"),"sigma of control parameter F normal distribution")
				("ini_f_type",po::value<std::string>(&ini_f_type)->default_value("1"),"control parameter F initialization type")
				("ini_f_uni_low_bound",po::value<std::string>(&ini_f_uni_low_bnd)->default_value("0.0"),"lower bound of initial uniform distribution of control parameter F")
				("ini_f_uni_up_bound",po::value<std::string>(&ini_f_uni_up_bnd)->default_value("1.0"),"upper bound of initial uniform distribution of control parameter F")
				("ini_f_mean",po::value<std::string>(&ini_f_mean)->default_value("0.5"),"mean of initial normal distribution of control parameter F")
				("ini_f_sigma",po::value<std::string>(&ini_f_sigma)->default_value("0.15"),"sigma of initial normal distribution of control parameter F")
                ("learn_p",po::value<std::string>(&learn_p)->default_value("50"),"learning period coefficient of control parameter")
				("f_per_dim",po::value<std::string>(&f_per_dim)->default_value("0"),"control parameter F granularity")
				("f_lower_bound",po::value<std::string>(&f_low_bnd)->default_value("0.1"),"lower bound of control parameter F")
				("f_upper_bound",po::value<std::string>(&f_up_bnd)->default_value("0.9"),"upper bound of control parameter F")
				("tau_1",po::value<std::string>(&tau_1)->default_value("0.1"),"control parameter F update probability")
				("tau_2",po::value<std::string>(&tau_2)->default_value("0.1"),"control parameter Pr update probability")
				("pr_stra",po::value<std::string>(&pr_stra)->default_value("1"),"control parameter Pr update strategy")
				("delta",po::value<std::string>(&delta)->default_value("0.5"),"probability of DE offspring generation scheme")
				("m_ratio",po::value<std::string>(&m_ratio)->default_value("0.5"),"ratio of parents population")
				;
		}
		int ReadPara(std::string conf_path="");

		double GetPr() {return GetSeqVal<double>("pr_val");}
		double GetF() {return GetSeqVal<double>("f_val");}
		double GetPrMean() {return GetSeqVal<double>("pr_mean");}
		double GetPrSigma() {return GetSeqVal<double>("pr_sigma"); }
		double GetFMean() {return GetSeqVal<double>("f_mean");}
		double GetFSigma() {return GetSeqVal<double>("f_sigma");}
		int GetIniFType() {return GetSeqVal<int>("ini_f_type");}
		double GetIniFUniLowBound() {return GetSeqVal<double>("ini_f_uni_low_bound");}
		double GetIniFUniUpBound() {return GetSeqVal<double>("ini_f_uni_up_bound");}
		double GetIniFMean() {return GetSeqVal<double>("ini_f_mean");}
		double GetIniFSigma() {return GetSeqVal<double>("ini_f_sigma");}
        int GetLearnPeriod() {return GetSeqVal<int>("learn_p");}
		int GetFPerDim() {return GetSeqVal<int>("f_per_dim");}
		double GetFLowerBound() {return GetSeqVal<double>("f_lower_bound");}
		double GetFUpperBound() {return GetSeqVal<double>("f_upper_bound");}
		double GetTau1() {return GetSeqVal<double>("tau_1");}
		double GetTau2() {return GetSeqVal<double>("tau_2");}
		int GetPrStra() {return GetSeqVal<int>("pr_stra");}
		double GetDelta() {return GetSeqVal<double>("delta");}
		double GetMRatio(){return GetSeqVal<double>("m_ratio");}
	private:
		std::string pr_val;// basic DE SPECIFIC,crossover probability fixed value
		std::string f_val;// basic DE SPECIFIC,amplification factor F fixed value
		std::string pr_mean;// crossover probability normal distribution mean
		std::string pr_sigma;// crossover probability normal distribution standard deviation:sigma
		std::string ini_f_type;// MY sde SPECIFIC,control parameter F initialization type
		std::string ini_f_uni_low_bnd;// MY sde SPECIFIC,lower bound of initial uniform distribution of control parameter F
		std::string ini_f_uni_up_bnd;// MY sde SPECIFIC,upper bound of initial uniform distribution of control parameter F
		std::string ini_f_mean;// sde/MY sde SPECIFIC,initial value of control parameter F normal distribution mean
		std::string ini_f_sigma;// sde/MY sde SPECIFIC,initial value of control parameter F normal distribution sigma
		std::string f_mean; // control parameter F normal distribution mean
		std::string f_sigma;// control parameter F normal distribution standard deviation:sigma
        std::string learn_p;// MY sde SPECIFIC,learning frequency/span
		std::string f_per_dim;// MY sde SPECIFIC,control parameter F granularity
		std::string f_low_bnd;// jde SPECIFIC,lower bound of control parameter F
		std::string f_up_bnd;// jde SPECIFIC,upper bound of control parameter F
		std::string tau_1;// jde SPECIFIC,control parameter F update probability
		std::string tau_2;// jde SPECIFIC,control parameter Pr update probability
		std::string pr_stra;// MY sde SPECIFIC,control parameter Pr update strategy
		std::string delta;// DE/EDA SPECIFIC,probability of DE offspring generation scheme
		std::string m_ratio;// DE/EDA SPECIFIC,ratio of parents population
	};
}// namespace de
#endif
