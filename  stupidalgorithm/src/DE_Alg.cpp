#include "de_alg.h"
#include "initializer.h"
#include "rand_val.h"

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

namespace de
{
	void de_alg::de_stat::alloc(int max_run)
	{
		gen_f.resize(max_run);
		gen_pr.resize(max_run);
	}

	void de_alg::de_stat::initialize(int max_run)
	{
		int i;
		for ( i=0;i<max_run;i++ )
		{
			gen_f[i].clear();
			gen_pr[i].clear();
		}
	}

	void de_alg::initialize()
	{
		com_alg::initialize();

		// initialize de_stat
		int max_run=m_ppara->get_max_run();
		m_de_stat.alloc(max_run);
		m_de_stat.initialize(max_run);
	}

	void de_alg::calc_gen_avg_de_para()
	{
		// find minimum gen in all runs
		int max_run=m_ppara->get_max_run();
		int min_gen;
		int i;
		int cur_gen_num;
		min_gen=m_de_stat.gen_f[0].size();
		for ( i=1;i<max_run;i++ )
		{
			cur_gen_num=m_de_stat.gen_f[i].size();
			if ( cur_gen_num < min_gen )
				min_gen=cur_gen_num;
		}

		// calculate average generational  F and Pr value
		m_de_stat.gen_avg_f.resize(min_gen,0.0);
		m_de_stat.gen_avg_pr.resize(min_gen,0.0);
		int j;
		for ( j=0;j<min_gen;j++ )
		{
			for ( i=0;i<max_run;i++ )
			{
				m_de_stat.gen_avg_f[j] += m_de_stat.gen_f[i][j];
				m_de_stat.gen_avg_pr[j] += m_de_stat.gen_pr[i][j];
			}// for every run
		}// for every generation
		m_de_stat.gen_avg_f /= max_run;
		m_de_stat.gen_avg_pr /= max_run;
	}// end function calc_gen_avg_de_para


	void de_alg::update_pop(population &pop,const population &trial_pop)
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

	void de_alg::calc_de_para_stat(const population& pop)
	{
		int pop_size=m_ppara->get_pop_size();
		int num_dims=m_ppara->get_dim();
		int f_per_dim=m_ppara->get_f_per_dim();
		int f_dims=(f_per_dim ? num_dims : 1);
		int i,j;
		int f_total_num=(f_per_dim ? pop_size*num_dims : pop_size);
		// calculate mean of f,pr
		m_de_stat.f_mean=0.0;
		m_de_stat.pr_mean=0.0;
		for ( i=0;i<pop_size;i++ )
		{
			m_de_stat.pr_mean += pop[i].stra[pr][0];//  pop[i].stra[pr][0];
			for ( j=0;j<f_dims;j++ )
				m_de_stat.f_mean +=  pop[i].stra[f][j]; // pop[i].stra[f][j];
		}
		m_de_stat.f_mean /= f_total_num;
		m_de_stat.pr_mean /= pop_size;
		// calculate standard deviation
		m_de_stat.f_std=m_de_stat.pr_std=0.0;
		for ( i=0;i<pop_size;i++ )
		{
			double diff;
			diff=pop[i].stra[pr][0]-m_de_stat.pr_mean;
			m_de_stat.pr_std += diff*diff;
			if ( false==f_per_dim )
			{
				diff=pop[i].stra[f][0]-m_de_stat.f_mean;
				m_de_stat.f_std += diff*diff;
			}
			else
			{
				for ( j=0;j<num_dims;j++ )
				{
					diff=pop[i].stra[f][j]-m_de_stat.f_mean;
					m_de_stat.f_std += diff*diff;
				}
			}
		}
		m_de_stat.pr_std /= pop_size;
		m_de_stat.pr_std=sqrt(m_de_stat.pr_std);
		m_de_stat.f_std /= f_total_num;
		m_de_stat.f_std=sqrt(m_de_stat.f_std);
	}// end function calc_de_para_stat

	inline void de_alg::record_avg_f(int cur_run)
	{
		m_de_stat.gen_f[cur_run].push_back(m_de_stat.f_mean);
	}
	inline void de_alg::record_avg_pr(int cur_run)
	{
		m_de_stat.gen_pr[cur_run].push_back(m_de_stat.pr_mean);
	}

	void de_alg::record_de_para_stat(int cur_run)
	{
		record_avg_f(cur_run);
		record_avg_pr(cur_run);
	}

	// record average gbest of every generation in all runs
	void de_alg::write_avg_de_para_per_gen()
	{
		ofstream avg_f_file(m_avg_f_path.c_str());
		ofstream avg_pr_file(m_avg_pr_path.c_str());
		print_gen_avg_val(avg_f_file,m_de_stat.gen_avg_f);
		print_gen_avg_val(avg_pr_file,m_de_stat.gen_avg_pr);
	}

	void de_alg::stat_run(population &pop,int cur_run)
	{
		com_alg::stat_run(pop,cur_run);
		int max_run=m_ppara->get_max_run();
		if ( cur_run==max_run-1 )
			calc_gen_avg_de_para();
	}

	// SPECIAL boundaries check
	void de_alg::bound_check(double &x,int dim)
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

	void de_alg::print_run_title(ostream &os)
	{
		output_base::print_run_title(os);
		os<<"\t"
			<<"\t"
			<<"avg_radius"
			<<"\t"
			<<"mean_of_F"
			<<"\t"
			<<"std_of_F"
			<<"\t"
			<<"mean_of_Pr"
			<<"\t"
			<<"std_of_Pr";
	}// end function print_run_title

	void de_alg::print_gen_stat(ostream &os,int cur_gen,const alg_stat &alg_stat)
	{
		output_base::print_gen_stat(os,cur_gen,alg_stat);
		// os.setf(ios::fixed,ios::floatfield);
		os.precision(8);
		os<<"\t"
			<<"\t"
			<<m_alg_stat.avg_radius;
		os.precision(4);
		os<<"\t"
			<<m_de_stat.f_mean
			<<"\t"
			<<m_de_stat.f_std
			<<"\t"
			<<m_de_stat.pr_mean
			<<"\t"
			<<m_de_stat.pr_std;

	}// end function print_gen_stat

	void de_alg::write_stat_vals()
	{
		com_alg::write_stat_vals();
		write_avg_de_para_per_gen();
	}

	int de_alg::run()
	{
		if ( !m_ppara )
			return -1;

		timer elapsed_t;
		// retrieve algorithm parameters
		int pop_size=m_ppara->get_pop_size();
		int num_dims=m_ppara->get_dim();
		double vtr=m_ppara->get_vtr();
		double pr_val,f_val;
		pr_val=m_ppara->get_pr();
		f_val=m_ppara->get_f();

		int m_cur_run;
		int max_run=m_ppara->get_max_run();// run/trial number
		// shared_ptr<progress_display> pprog_dis;// algorithm progress indicator from boost
		// alloc_prog_indicator(pprog_dis);

		// allocate original pop and trial pop
		population pop(pop_size);
		allocate_pop(pop,num_dims,stra_num);
		population trial_pop(pop_size);
		allocate_pop(trial_pop,num_dims,stra_num);

		// generate algorithm statistics output file name
		ofstream stat_file(m_com_out_path.stat_path.c_str());
		// allocate stop condition object dynamically
		alloc_stop_cond();

		vector<int> vec_idx(pop_size-1);
		int shuffle_size=pop_size-1;

		// random U(0,1) generator
		uniform_01<> dist_01;
		variate_generator<mt19937&, uniform_01<> > rnd_01(gen, dist_01);

		// generator for random DIMENSION index
		uniform_int<> dist_dim(0,num_dims-1);
		variate_generator<mt19937&, uniform_int<> > rnd_dim_idx(gen, dist_dim);

		// iteration start
		for ( m_cur_run=0;m_cur_run<max_run;m_cur_run++ )
		{
			reset_run_stat();
			m_de_stat.reset();
			int z;
			for ( z=0;z<pop_size;z++ )
			{
				pop[z].stra[f].assign(1,f_val);
				pop[z].stra[pr].assign(1,pr_val);
			}
			set_orig_pop(pop);
			update_diversity(pop);

			record_gen_vals(m_alg_stat,m_cur_run);
			calc_de_para_stat(pop);
			record_de_para_stat(m_cur_run);

			print_run_times(stat_file,m_cur_run+1);
			print_run_title(stat_file);
			// output original population statistics
			print_gen_stat(stat_file,1,m_alg_stat);

			m_cur_gen=1;
			while ( false==(*m_pstop_cond) ) // for every iteration
			{
				m_de_stat.reset();
				int rnd_dim;
				double dim_mut_chance;
				int i,j,k;

				trial_pop=pop;// operator =
				for ( i=0;i<pop_size;i++ )
				{
					// generating three mutually different individual index other than i using random shuffle
					// initialize index vector
					for ( k=0;k<shuffle_size;k++ )
					{
						if ( k<i )
							vec_idx[k]=k;
						else
							vec_idx[k]=(k+1)%pop_size;// EXCLUDE i
					}
					// random shuffle
					for ( k=0;k<shuffle_size;k++ )
					{
						// generator for random SHUFFLE VECTOR index
						uniform_int<> dist_uni_shuf(k,shuffle_size-1);
						variate_generator<mt19937&, uniform_int<> > rnd_shuf_idx(gen, dist_uni_shuf);
						int idx_tmp=rnd_shuf_idx();
						swap(vec_idx[k],vec_idx[idx_tmp]);
					}
					int i1,i2,i3;// i!=i1!=i2!=i3
					i1=vec_idx[0];
					i2=vec_idx[1];
					i3=vec_idx[2];

					// generating F and Pr from specified normal distribution and stat its value
					rnd_dim=rnd_dim_idx();
					
					for ( j=0;j<num_dims;j++ )
					{
						dim_mut_chance=rnd_01();
						if ( rnd_dim==j || dim_mut_chance<=pr_val )
						{
							trial_pop[i].x[j]=trial_pop[i1].x[j]+f_val*(trial_pop[i2].x[j]-trial_pop[i3].x[j]);
							// boundaries check
							bound_check(trial_pop[i].x[j],j);
						}
					}// for every dimension
				}// for every particle
				// evaluate pop
				eval_pop(trial_pop,*m_pfunc,m_alg_stat);
				update_pop(pop,trial_pop);
				stat_pop(pop,m_alg_stat);
				update_search_radius();
				update_diversity(pop);

				calc_de_para_stat(pop);
				record_de_para_stat(m_cur_run);
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

}// end namespace de
