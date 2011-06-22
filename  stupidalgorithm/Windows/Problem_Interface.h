#ifndef JERRYLIU_PSO_PROBLEM_INTERFACE
#define JERRYLIU_PSO_PROBLEM_INTERFACE

#include "Algo_DataStruct.h"
#include "FuncDef.h"

class Problem_Interface
{
public:
	// evaluation function
	void eval_ind( Individual &ind,
				  Func_Interface &func_eval,
				   Alg_Stat &alg_stat
				 );
	void eval_pop( Population &pop,
				  Func_Interface &func_eval,
				   Alg_Stat &alg_stat
				 );
	void eval_ini_pop( Population &pop,
					  Func_Interface &func_eval,
					   Alg_Stat &alg_stat
					 );
	virtual void stat_ini_pop(Population &pop,Alg_Stat &alg_stat);
	virtual void stat_pop(const Population &pop,Alg_Stat &alg_stat);

	void record_bst_so_far(Alg_Stat &alg_stat,int cur_run); 
	void record_diver(Alg_Stat &alg_stat,int cur_run);
	void record_radius(Alg_Stat &alg_stat,int cur_run);
	virtual void record_gen_vals(Alg_Stat &alg_stat,int cur_run);
};



#endif