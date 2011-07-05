#ifndef STUPIDALGO_DE_ALL_PARAMETERS
#define STUPIDALGO_DE_ALL_PARAMETERS

#include <fstream>
#include <iostream>

#include "para_base.h"

namespace de
{
	class de_para:public para_base
	{
	public:
		// ctor
		de_para(std::string conf_path="")
			:para_base(conf_path) 
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

				("max_archive",po::value<std::string>(&max_archive)->default_value("50"),"maximum size of external elitist archive")
				;
		}
		virtual ~de_para(){}

		int read_para(std::string conf_path="");

		double get_pr() {return get_seq_val<double>("pr_val");}
		double get_f() {return get_seq_val<double>("f_val");}
		double get_pr_mean() {return get_seq_val<double>("pr_mean");}
		double get_pr_sigma() {return get_seq_val<double>("pr_sigma"); }
		double get_f_mean() {return get_seq_val<double>("f_mean");}
		double get_f_sigma() {return get_seq_val<double>("f_sigma");}
		int get_ini_f_type() {return get_seq_val<int>("ini_f_type");}
		double get_ini_f_uni_low_bnd() {return get_seq_val<double>("ini_f_uni_low_bound");}
		double get_ini_f_uni_up_bnd() {return get_seq_val<double>("ini_f_uni_up_bound");}
		double get_ini_f_mean() {return get_seq_val<double>("ini_f_mean");}
		double get_ini_f_sigma() {return get_seq_val<double>("ini_f_sigma");}
        int get_learn_period() {return get_seq_val<int>("learn_p");}
		int get_f_per_dim() {return get_seq_val<int>("f_per_dim");}
		double get_f_low_bnd() {return get_seq_val<double>("f_lower_bound");}
		double get_f_up_bnd() {return get_seq_val<double>("f_upper_bound");}
		double get_tau_1() {return get_seq_val<double>("tau_1");}
		double get_tau_2() {return get_seq_val<double>("tau_2");}
		int get_pr_stra() {return get_seq_val<int>("pr_stra");}
		double get_delta() {return get_seq_val<double>("delta");}
		double get_m_ratio(){return get_seq_val<double>("m_ratio");}

		int get_max_archive() {return get_seq_val<int>("max_archive");}
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

		std::string max_archive;// MOSADE SPECIFIC,maximum size of external elitist archive
	};
}// namespace de
#endif
