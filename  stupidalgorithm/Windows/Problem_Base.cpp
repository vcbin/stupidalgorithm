#include "problem_base.h"

// #include <algorithm> // for std::min_element

void problem_base::eval_ind( individual &ind,
										func_base &func_eval,
										 alg_stat &alg_stat
									   )
{
		func_eval(ind);// polymorphic invocation
		alg_stat.eval_num++;
		alg_stat.all_eval_num++;
}

void problem_base::stat_ini_pop(population &pop,alg_stat &alg_stat)
{
	alg_stat.gbest.obj[0]=INF;
	size_t pop_size=pop.size();
	size_t num_dims=pop[0].x.size();

	alg_stat.avg=0.0;
	bool stag_flag=true;
	size_t i;
	for ( i=0;i<pop_size;i++ )
	{
		// stat the pop step by step
		alg_stat.avg += pop[i].obj[0];

		if ( pop[i] < alg_stat.gbest )
		{
			alg_stat.gbest=pop[i];
			stag_flag=false;
			alg_stat.num_stag=0;
		}
	}// for every particle
	alg_stat.avg /=  pop_size;// average fitness value

	if (stag_flag)
		alg_stat.num_stag++;// increment stagnation count
	// calclate standard deviation
	alg_stat.std=0.0;
	for ( i=0;i<pop_size;i++ )
	{
		double diff=pop[i].obj[0]-alg_stat.avg;
		alg_stat.std += diff*diff;
	}
	alg_stat.std /=  pop_size;
	alg_stat.std=sqrt(alg_stat.std);

	alg_stat.pbest=pop;
}

void problem_base::eval_ini_pop( population &pop,
									 func_base &func_eval,
									  alg_stat &alg_stat
									 )
{
	size_t pop_size=pop.size();
	unsigned i;

	for ( i=0;i<pop_size;i++ )
	{
		func_eval(pop[i]);// polymorphic invocation
		alg_stat.eval_num++;
		alg_stat.all_eval_num++;
	}//  end evaluation
}// end function eval_ini_pop


void problem_base::eval_pop( population &pop,
								 func_base &func_eval,
								  alg_stat &alg_stat
								 )
{
	size_t pop_size=pop.size();

	unsigned i;
	for ( i=0;i<pop_size;i++ )
	{
		func_eval(pop[i]);// polymorphic invocation
		alg_stat.eval_num++;
		alg_stat.all_eval_num++;
	}//  end evaluation
}

void problem_base::stat_pop(const population &pop,alg_stat &alg_stat)
{
	size_t pop_size=pop.size();

	alg_stat.avg=0.0;
	bool stag_flag=true;
	unsigned i;
	for ( i=0;i<pop_size;i++ )
	{
		// stat the pop step by step
		alg_stat.avg += pop[i].obj[0];
		
		if ( pop[i] < alg_stat.pbest[i] )
			alg_stat.pbest[i]=pop[i];// operator=
		if ( pop[i] < alg_stat.gbest )
		{
			alg_stat.gbest=pop[i];
			stag_flag=false;
			alg_stat.num_stag=0;
		}
	}// for every stating particle

	if (stag_flag)
		alg_stat.num_stag++;// increment stagnation count
	alg_stat.avg /=  pop_size;// average fitness value

	// calclate standard deviation
	alg_stat.std=0.0;
	for ( i=0;i<pop_size;i++ )
	{
		double diff=pop[i].obj[0]-alg_stat.avg;
		alg_stat.std += diff*diff;
	}
	alg_stat.std /=  pop_size;
	alg_stat.std=sqrt(alg_stat.std);
}

void problem_base::record_bst_so_far(alg_stat &alg_stat,int cur_run)
{
	alg_stat.gen_bst_val[cur_run].push_back(alg_stat.gbest.obj[0]);
}

void problem_base::record_diver(alg_stat &alg_stat,int cur_run)
{
	alg_stat.gen_div_val[cur_run].push_back(alg_stat.pos_diver);
}

void problem_base::record_radius(alg_stat &alg_stat,int cur_run)
{
	alg_stat.gen_rad_val[cur_run].push_back(alg_stat.avg_radius);
}

void problem_base::record_gen_vals(alg_stat &alg_stat,int cur_run)
{
	record_bst_so_far(alg_stat,cur_run);
	record_diver(alg_stat,cur_run);
	record_radius(alg_stat,cur_run);
}
