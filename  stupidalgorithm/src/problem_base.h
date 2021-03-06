#ifndef STUPIDALGO_SOP_PROBLEM_INTERFACE
#define STUPIDALGO_SOP_PROBLEM_INTERFACE

#include "algo_datastruct.h"
#include "funcdef.h"

using namespace benchmark;

class problem_base
{
public:
	virtual ~problem_base() {}
	// evaluation function
	void eval_ind( individual &ind,
				  func_base &func_eval,
				   alg_stat &alg_stat
				 );
	void eval_pop( population &pop,
				  func_base &func_eval,
				   alg_stat &alg_stat
				 );
	void eval_ini_pop( population &pop,
					  func_base &func_eval,
					   alg_stat &alg_stat
					 );
	virtual void stat_ini_pop(population &pop,alg_stat &alg_stat);
	virtual void stat_pop(const population &pop,alg_stat &alg_stat);

	void record_bst_so_far(alg_stat &alg_stat,int cur_run); 
	void record_diver(alg_stat &alg_stat,int cur_run);
	void record_radius(alg_stat &alg_stat,int cur_run);
	virtual void record_gen_vals(alg_stat &alg_stat,int cur_run);
};// end class problem_base



#endif
