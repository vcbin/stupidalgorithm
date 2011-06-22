#ifndef JERRYLIU_PSO_OUTPUT_INTERFACE
#define JERRYLIU_PSO_OUTPUT_INTERFACE

#include <iostream>
#include <boost/shared_ptr.hpp>
#include "Algo_DataStruct.h"
#include "Alg_Interface.h"

class Output_Interface
{
public:
	void Print_Run_Title(std::ostream &os);
	inline void Print_Run_Times(std::ostream &os,int m_cur_run) { os<<"\n"<<"run="<<m_cur_run; }
	void Print_Sep_Line(std::ostream &os);
	void Print_Gen_Stat(std::ostream &os,int cur_gen,const Alg_Stat &m_alg_stat);
	void Print_Run_Stat_Title(std::ostream &os);
	void Print_Run_Stat(std::ostream &os,const Alg_Stat &m_alg_stat,int max_run);
	void Print_Gen_Avg_Vals(std::ostream &os,const std::vector<double> &gen_vals);
	void Print_Gen_Avg_Vec(std::ostream &os,const std::vector<std::vector<double> > &gen_vec);
	void Print_Best_X(std::ostream &os,const Individual &bst_ind);
	void Print_All_Gbest_Vals(std::ostream &os,const std::vector<double> &gbest_vec);
	void Print_Avg_Time(std::ostream &os,double sec);
protected:
	void write_initial_pop(std::ostream &os,const Population &pop);
	inline void print_gen_val(std::ostream &os,int cur_gen,double value) { os<<cur_gen<<"\t"<<value<<"\n"; }
	void print_gen_vec(std::ostream &os,int cur_gen,const std::vector<double> &vec);
};// end class Output_Interface declaration

void Print_Algo_Stat_Title(std::ostream &os,bool has_max_gen);
void Print_Para_Stat_Title(std::ostream &os,std::string para_var_list,std::string para_num_str,bool has_max_gen);
void Print_Common_Filelist_Title(std::ostream &os);
inline void Print_FileName( std::ostream &os,int func_type,int algo_type,std::string file_name ) { os<<func_type<<"\t"<<algo_type<<"\t"<<file_name<<"\n"; }
#endif