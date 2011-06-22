#include "JDE.h"
#include "Rand_Val.h"

using std::vector;
using std::cout;
using std::string;
using std::ofstream;
using std::ostream;
using std::swap;
using std::endl;

using boost::shared_ptr;
using boost::progress_display;
using boost::timer;

using boost::uniform_01;
using boost::uniform_int;
using boost::normal_distribution;
using boost::mt19937;
using boost::variate_generator;

extern boost::mt19937 gen;

namespace de
{
	namespace jde
	{
		int jde_alg::run()
		{
			if ( !m_ppara )
				return -1;

			timer elapsed_t;
			// retrieve algorithm parameters
			size_t pop_size=m_ppara->get_pop_size();
			size_t num_dims=m_ppara->get_dim();
			int max_gen=m_ppara->get_max_gen();
			double vtr=m_ppara->get_vtr();
			double f_low_bnd,f_up_bnd;
			f_low_bnd=m_ppara->get_f_low_bnd();
			f_up_bnd=m_ppara->get_f_up_bnd();
			double tau_1,tau_2;
			tau_1=m_ppara->get_tau_1();
			tau_2=m_ppara->get_tau_2();
			
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
			ofstream stat_file(m_com_out_path.stat_path);
			// allocate stop condition object dynamically
			bool run_once=(1==max_run);
			alloc_stop_cond();

			size_t shuffle_size=pop_size-1;
			vector<int> vec_idx1(shuffle_size);

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
				// jde SPECIFIC,initialize F value vector EVERY single run
				size_t i;
				for ( i=0;i<pop_size;i++ )
				{
					pop[i].stra[pr].assign(1,0.0);
					pop[i].stra[f].assign(1,0.0);
				}
				set_orig_pop(pop);
				update_diversity(pop);

				calc_de_para_stat(pop);
				record_de_para_stat(m_cur_run);
				record_gen_vals(m_alg_stat,m_cur_run);

				print_run_times(stat_file,m_cur_run+1);
				print_run_title(stat_file);
				// output original population statistics
				print_gen_stat(stat_file,1,m_alg_stat);

				m_cur_gen=1;
				while ( false==(*m_pstop_cond) ) // for every iteration
				{
					m_de_stat.reset();
					int rnd_dim;
					double f_i;
					double pr_i;
					double dim_mut_chance;
					size_t i,j,k;

					trial_pop=pop;// operator =
					for ( i=0;i<pop_size;i++ )
					{
						// generating three mutually different individual index other than i using random shuffle
						// initialize index vector
						for ( k=0;k<shuffle_size;k++ )
						{
							if ( k<i )
								vec_idx1[k]=k;
							else
								// EXCLUDE i
								vec_idx1[k]=(k+1)%pop_size;
						}
						// random shuffle
						for ( k=0;k<shuffle_size;k++ )
						{
							// generator for random SHUFFLE VECTOR index
							uniform_int<> dist_uni_shuf(k,shuffle_size-1);
							variate_generator<mt19937&, uniform_int<> > rnd_shuf_idx(gen, dist_uni_shuf);
							int idx_tmp=rnd_shuf_idx();
							swap(vec_idx1[k],vec_idx1[idx_tmp]);
						}

						int i1,i2,i3;// i!=i1!=i2!=i3
						i1=vec_idx1[0];
						i2=vec_idx1[1];
						i3=vec_idx1[2];

						rnd_dim=rnd_dim_idx();

						double pr_chance=rnd_01();
						if ( pr_chance<=tau_2 )
						{
							pr_i=rnd_01();
							trial_pop[i].stra[pr][0]=pr_i;
						}
						else
							pr_i=trial_pop[i].stra[pr][0];
						
						// scaling factor F self-adaptive update equation
						double f_chance=rnd_01();
						if ( f_chance<=tau_1 )
						{
							f_i=f_low_bnd+rnd_01()*f_up_bnd;
							trial_pop[i].stra[f][0]=f_i;
						}
						else
							f_i=trial_pop[i].stra[f][0];// keep unchanged at thsi iteration

						for ( j=0;j<num_dims;j++ )
						{
							dim_mut_chance=rnd_01();
							if ( rnd_dim==j || dim_mut_chance<=pr_i )
							{
								trial_pop[i].x[j]=trial_pop[i1].x[j]+f_i*(trial_pop[i2].x[j]-trial_pop[i3].x[j]);
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
				}// while single run stop_condition is false

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

	}// end namespace jde
}// end namespace de
