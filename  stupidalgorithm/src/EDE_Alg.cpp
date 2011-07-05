#include "ede_alg.h"

using std::partial_sort;

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
	namespace ede
	{
		int ede_alg::run()
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

			int shuffle_size;
			vector<int> vec_idx;

			// random U(0,1) generator
			uniform_01<> dist_01;
			variate_generator<mt19937&, uniform_01<> > rnd_01(gen, dist_01);

			// generator for random DIMENSION index
			uniform_int<> dist_dim(0,num_dims-1);
			variate_generator<mt19937&, uniform_int<> > rnd_dim_idx(gen, dist_dim);

			// generator for random INDIVIDUAL index
			uniform_int<> dist_ind(0,pop_size-1);
			variate_generator<mt19937&, uniform_int<> > rnd_ind_idx(gen, dist_ind);

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
					int d,i,j,k;

					trial_pop=pop;// operator =
					for ( i=0;i<pop_size;i++ )
					{
						do
						{
							d=rnd_ind_idx();
						}
						while ( pop[d]>pop[i] );
						// generating three mutually different individual index other than i using random shuffle
						// initialize index vector
						int equ;
						if ( i!=d )
						{
							shuffle_size=pop_size-2;
							equ=false;
						}
						else
						{
							shuffle_size=pop_size-1;
							equ=true;
						}
						vec_idx.resize(shuffle_size);
						for ( k=0;k<shuffle_size;k++ )
						{
							int less,more;
							if ( !equ )
							{
								less=(i<d?i:d);
								more=(i>d?i:d);
							}
							else
								less=i;

							if ( k<less )
								vec_idx[k]=k;
							else
							{
								vec_idx[k]=(k+1)%pop_size;// EXCLUDE i
								if ( !equ )
								{
									if ( vec_idx[k]>=more )
										vec_idx[k]=(vec_idx[k]+1)%pop_size;// EXCLUDE d
								}
							}
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
						int i2,i3;// i!=d!=i2!=i3
						i2=vec_idx[1];
						i3=vec_idx[2];

						// generating F and Pr from specified normal distribution and stat its value
						rnd_dim=rnd_dim_idx();

						for ( j=0;j<num_dims;j++ )
						{
							dim_mut_chance=rnd_01();
							if ( rnd_dim==j || dim_mut_chance<=pr_val )
							{
								trial_pop[i].x[j]=0.5*(trial_pop[d].x[j]+trial_pop[i].x[j]) +
									f_val*( trial_pop[d].x[j]-trial_pop[i].x[j]+trial_pop[i2].x[j]-trial_pop[i3].x[j] );
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
		}
	}// namespace ede
}// namespace de
