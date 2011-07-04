#include "EDA_Alg.h"
#include "Initializer.h"
#include "Rand_Val.h"

using std::vector;
using std::cout;
using std::string;
using std::ofstream;
using std::ostream;
using std::swap;
using std::endl;
using std::ostream;
using std::ios;

using boost::shared_ptr;
// using boost::progress_display;
using boost::timer;

using boost::uniform_01;
using boost::uniform_int;
using boost::normal_distribution;
using boost::mt19937;
using boost::variate_generator;

extern boost::mt19937 gen;

namespace eda
{
	void eda_alg::initialize()
	{
		com_alg::initialize();

		//
		int num_dims=m_ppara->get_dim();
		m_x_mean.resize(num_dims,0.0);
		m_x_std.resize(num_dims,0.0);
	}

	void eda_alg::update_pop(population &pop,const population &trial_pop)
	{
		int pop_size=pop.size();
		int num_dims=m_ppara->get_dim();
		int i;
		for ( i=0;i<pop_size;i++ )
		{
			if ( trial_pop[i]<pop[i] )
			{
				m_alg_stat.delta_x[i]=trial_pop[i].x-pop[i].x;// update step size delta_x
				pop[i]=trial_pop[i];
			}
			else
				m_alg_stat.delta_x[i].assign(num_dims,0.0);// step size==0
		}// for every individual
	}// end function update_pop


	// SPECIAL boundaries check
	void eda_alg::bound_check(double &x,int dim)
	{
		const val_range& val_bounds=m_ppara->get_val_bnd();
		bool out_up_bnd,out_low_bnd;
		double low_bnd,high_bnd;
		low_bnd=val_bounds[dim].min_val;
		high_bnd=val_bounds[dim].max_val;
		out_low_bnd=x < low_bnd;
		out_up_bnd=x > high_bnd;
		if ( out_up_bnd || out_low_bnd )
		{
			m_alg_stat.all_ob_num++;
			/*normal_distribution<> norm_dist;
			variate_generator<mt19937&, normal_distribution<> > rnd_norm(gen, norm_dist);
			x=m_alg_stat.cen_ind[dim]+rnd_norm();*/

			// uniformly randomize x within [lower_bound,upper_bound]
			uniform_01<> dist;
			variate_generator<mt19937&, uniform_01<> > rnd_num(gen, dist);
			x=low_bnd+rnd_num()*(high_bnd-low_bnd);
			bound_check(x,dim);// RECURSIVE invocation
		}
	}// end function bound_check

	void eda_alg::stat_eda_dist(population &pop)
	{
		int pop_size=m_ppara->get_pop_size();
		int num_dims=m_ppara->get_dim();
		double m_ratio=m_ppara->get_m_ratio();
		int parent_size=static_cast<int>(pop_size*m_ratio);
		typedef population::iterator pop_itr;
		pop_itr itr_middle=pop.begin()+parent_size;
		pop_itr itr_begin=pop.begin();
		partial_sort(itr_begin,itr_middle,pop.end());// best parent_size in second parent_size
		// calc mean of best parent_size
		m_x_mean.assign(num_dims,0.0);
		pop_itr itr_cur;
		for ( itr_cur=itr_begin;itr_cur!=itr_middle;++itr_cur )
			m_x_mean += itr_cur->x;
		m_x_mean /= parent_size;
		// calc std of best parent_size
		int j;
		m_x_std.assign(num_dims,0.0);
		for ( itr_cur=itr_begin;itr_cur!=itr_middle;++itr_cur )
		{
			for ( j=0;j<num_dims;j++ )
			{
				double diff=(itr_cur->x[j]-m_x_mean[j]);
				m_x_std[j] += diff*diff;
			}
		}
		m_x_std /= parent_size;
		for ( j=0;j<num_dims;j++ )
			m_x_std[j]=sqrt(m_x_std[j]);
	}// end function stat_eda_dist

	double eda_alg::gen_eda_x(int dim)
	{
		normal_distribution<> dist_norm(m_x_mean[dim],m_x_std[dim]);
		variate_generator<mt19937&, normal_distribution<> > rnd_x(gen, dist_norm);
		double x=rnd_x();
		bound_check(x,dim);
		return x;
	}// end function gen_eda_ind

	int eda_alg::run()
	{
		if ( !m_ppara )
			return -1;

		timer elapsed_t;
		// retrieve algorithm parameters
		int pop_size=m_ppara->get_pop_size();
		int num_dims=m_ppara->get_dim();
		double vtr=m_ppara->get_vtr();

		int m_cur_run;
		int max_run=m_ppara->get_max_run();// run/trial number
		// shared_ptr<progress_display> pprog_dis;// algorithm progress indicator from boost
		// alloc_prog_indicator(pprog_dis);

		// allocate original pop and trial pop
		population pop(pop_size);
		allocate_pop(pop,num_dims);
		population trial_pop(pop_size);
		allocate_pop(trial_pop,num_dims);

		// generate algorithm statistics output file name
		ofstream stat_file(m_com_out_path.stat_path.c_str());
		// allocate stop condition object dynamically
		alloc_stop_cond();


		// random U(0,1) generator
		uniform_01<> dist_01;
		variate_generator<mt19937&, uniform_01<> > rnd_01(gen, dist_01);

		// iteration start
		for ( m_cur_run=0;m_cur_run<max_run;m_cur_run++ )
		{
			reset_run_stat();

			set_orig_pop(pop);
			update_diversity(pop);
			stat_eda_dist(pop);// eda SPECIFIC
			record_gen_vals(m_alg_stat,m_cur_run);

			print_run_times(stat_file,m_cur_run+1);
			print_run_title(stat_file);
			// output original population statistics
			print_gen_stat(stat_file,1,m_alg_stat);

			m_cur_gen=1;
			while ( false==(*m_pstop_cond) ) // for every iteration
			{
				int i,j;

				trial_pop=pop;// operator =
				for ( i=0;i<pop_size;i++ )
				{
					for ( j=0;j<num_dims;j++ )
					{
						trial_pop[i].x[j]=gen_eda_x(j);
						// boundaries check
						bound_check(trial_pop[i].x[j],j);
					}// for every dimension
				}// for every particle
				// evaluate pop
				eval_pop(trial_pop,*m_pfunc,m_alg_stat);
				update_pop(pop,trial_pop);
				stat_eda_dist(pop);// eda SPECIFIC
				stat_pop(pop,m_alg_stat);
				update_search_radius();
				update_diversity(pop);

				record_gen_vals(m_alg_stat,m_cur_run);

				print_gen_stat(stat_file,m_cur_gen+1,m_alg_stat);
				update_conv_stat(vtr);

				/*if ( run_once )
				++(*pprog_dis);*/

				m_cur_gen++;
			}// while single run termination criterion is not met

			// single run end
			stat_run(pop,m_cur_run);// stat single run for algorithm analysis
			if ( is_final_run(m_cur_run,max_run) )
				print_run_stat(stat_file,m_alg_stat,max_run);
			/*if ( !run_once )
			++(*pprog_dis);*/
		}// for every run

		print_avg_gen(stat_file,m_alg_stat.run_avg_gen);
		// stat and output average time per run by second
		m_alg_stat.run_avg_time=elapsed_t.elapsed();
		m_alg_stat.run_avg_time /= (max_run*1.0);
		print_avg_time(stat_file,m_alg_stat.run_avg_time);

		print_best_x(stat_file,m_alg_stat.bst_ind);
		write_stat_vals();

		cout<<endl;// flush cout output
		return 0;
	}// end function Run

}// end namespace eda
