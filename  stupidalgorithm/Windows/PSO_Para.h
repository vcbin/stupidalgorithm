#ifndef STUPIDALGO_PSO_ALL_PARAMETERS
#define STUPIDALGO_PSO_ALL_PARAMETERS

#include "para_base.h"

namespace pso
{
	class pso_para:public para_base
	{
	public:
		// ctor
		pso_para(std::string conf_path="")
			:para_base(conf_path)
		{
			general.add_options()
				("Vmax_type",po::value<std::string>(),"type of V_max")
				("Vmax_cof",po::value<std::string>(),"V_max related coefficient")
				("omega_max",po::value<std::string>(),"maximum value of the inertia weight coefficient omega")
				("omega_min",po::value<std::string>(),"minimum value of the inertia weight coefficient omega")
				("cr1",po::value<std::string>(),"value of speed update formula coefficient cr1")
				("cr2",po::value<std::string>(),"value of speed update formula coefficient cr2")
				("ob_type",po::value<std::string>(),"boundary condition type")
				("d_l",po::value<std::string>(&d_l)->default_value("5.0E-6"),"default value for lower bound of attration")
				("d_h",po::value<std::string>(&d_h)->default_value("0.25"),"default value for upper bound of attration");
		}
		int read_para(std::string conf_path="");

		inline double get_omega_max() {return get_seq_val<double>("omega_max");}
		inline double get_omega_min() {return get_seq_val<double>("omega_min");}
		inline double get_cr1() { return get_seq_val<double>("cr1"); }
		inline double get_cr2() {return get_seq_val<double>("cr2");}
		inline int get_vmax_type() {return get_seq_val<int>("Vmax_type");}
		inline double get_vmax_cof() {return get_seq_val<double>("Vmax_cof");}
		inline int get_ob_type() {return get_seq_val<int>("ob_type");}
		inline double get_dl() {return get_seq_val<double>("d_l");}
		inline double get_dh() {return get_seq_val<double>("d_h");}
	private:
		// literal value of PSO specific parameter
		std::string w_max;// omega_start/omega_max
		std::string w_min;// omega_end/omega_min
		std::string cr1;// speed coefficient 1
		std::string cr2;// speed coefficient 2
		std::string vmax_type;// type of V_max
		std::string vmax_cof;// V_max related coefficient
		std::string ob_type;// Out Bound handler type

		std::string d_l;// OPTIONAL,default value for lower bound of attration
		std::string d_h;// OPTIONAL,default value for upper bound of attration

	};// end class para

}// namespace pso

#endif