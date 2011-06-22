#include "Problem_Interface.h"

// #include <algorithm> // for std::min_element

void Problem_Interface::eval_ind( Individual &ind,
										Func_Interface &func_eval,
										 Alg_Stat &alg_stat
									   )
{
		func_eval(ind);// polymorphic invocation
		alg_stat.eval_num++;
		alg_stat.all_eval_num++;
}

void Problem_Interface::stat_ini_pop(Population &pop,Alg_Stat &alg_stat)
{
	alg_stat.gbest.obj[0]=std::numeric_limits<double>::max();
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

void Problem_Interface::eval_ini_pop( Population &pop,
									 Func_Interface &func_eval,
									  Alg_Stat &alg_stat
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


void Problem_Interface::eval_pop( Population &pop,
								 Func_Interface &func_eval,
								  Alg_Stat &alg_stat
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

void Problem_Interface::stat_pop(const Population &pop,Alg_Stat &alg_stat)
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

void Problem_Interface::record_bst_so_far(Alg_Stat &alg_stat,int cur_run)
{
	alg_stat.gen_bst_val[cur_run].push_back(alg_stat.gbest.obj[0]);
}

void Problem_Interface::record_diver(Alg_Stat &alg_stat,int cur_run)
{
	alg_stat.gen_div_val[cur_run].push_back(alg_stat.pos_diver);
}

void Problem_Interface::record_radius(Alg_Stat &alg_stat,int cur_run)
{
	alg_stat.gen_rad_val[cur_run].push_back(alg_stat.avg_radius);
}

void Problem_Interface::record_gen_vals(Alg_Stat &alg_stat,int cur_run)
{
	record_bst_so_far(alg_stat,cur_run);
	record_diver(alg_stat,cur_run);
	record_radius(alg_stat,cur_run);
}
