#include "NSDE_Alg.h"

#include "SDE_Alg.h"
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
using boost::uniform_real;
using boost::normal_distribution;
using boost::mt19937;
using boost::variate_generator;

extern boost::mt19937 gen;

namespace de
{
	namespace nsde
	{
		void nsde_alg::initialize()
		{
			bbde_alg::initialize();

			int pop_size=m_ppara->get_pop_size();
			m_succ_cr.resize(pop_size);
			int num_dims=m_ppara->get_dim();
			int f_per_dim=m_ppara->get_f_per_dim();
			int f_dims=(f_per_dim ? num_dims : 1);
			m_succ_f.resize(f_dims);
			int max_run=m_ppara->get_max_run();
			m_gen_bi_norm_stat.resize(max_run);
		}// end function initialize

		void nsde_alg::update_pop(population &pop,const population &trial_pop)
		{
			int pop_size=pop.size();
			int num_dims=m_ppara->get_dim();
			int learn_p=m_ppara->get_learn_period();
			int f_per_dim=m_ppara->get_f_per_dim();
			int f_dims=(f_per_dim ? num_dims : 1);
			int i;
			for ( i=0;i<pop_size;i++ )
			{
				if ( trial_pop[i]<pop[i] )
				{
					m_alg_stat.delta_x[i]=trial_pop[i].x-pop[i].x;// update step size delta_x
					pop[i]=trial_pop[i];
					m_succ_f.push_back(pop[i].stra[f][0]);// record succeeded F values
					m_succ_cr.push_back(pop[i].stra[pr][0]);// record succeeded CR values
				}
				else
					m_alg_stat.delta_x[i].assign(num_dims,0.0);// trial failed,step size delta_x == 0
			}// for every individual
		}// end function update_pop

		void nsde_alg::calc_gen_avg_bi_norm_stat()
		{
			// find minimum checkpoint count in all runs
			int max_run=m_ppara->get_max_run();
			int min_ch;
			int i;
			int cur_ch_num;
			min_ch=m_gen_bi_norm_stat[0].size();
			for ( i=1;i<max_run;i++ )
			{
				cur_ch_num=m_gen_bi_norm_stat[i].size();
				if ( cur_ch_num < min_ch )
					min_ch=cur_ch_num;
			}

			// calculate average generational  F and Pr value
			bi_norm_var var_tmp={0,0,0,0,0};
			m_avg_bi_norm_var.resize(min_ch,prob_stat(0,var_tmp));
			int j;
			for ( j=0;j<min_ch;j++ )
			{
				for ( i=0;i<max_run;i++ )
				{
					m_avg_bi_norm_var[j].gen=m_gen_bi_norm_stat[i][j].gen;
					m_avg_bi_norm_var[j].var += m_gen_bi_norm_stat[i][j].var;
				}// for every run
				m_avg_bi_norm_var[j].var /= max_run;
			}// for every generation
		}// end function calc_gen_avg_bi_norm_stat

		void nsde_alg::write_stat_vals()
		{
			bbde_alg::write_stat_vals();
			int i;
			int max_check=m_avg_bi_norm_var.size();
			ofstream out_file(m_gen_avg_bi_norm_stat_path);
			for ( i=0;i<max_check;i++ )
			{
				out_file<<"\n"<<m_avg_bi_norm_var[i].gen+1<<"\t"<<m_avg_bi_norm_var[i].var;
			}
		}// end function write_stat_vals

		void nsde_alg::stat_run(population &pop,int cur_run)
		{
			bbde_alg::stat_run(pop,cur_run);
			int max_run=m_ppara->get_max_run();
			if ( cur_run==max_run-1 )
				calc_gen_avg_bi_norm_stat();
		}

		int nsde_alg::run()
		{
			if ( !m_ppara )
				return -1;

			timer elapsed_t;
			// retrieve algorithm parameters
			int pop_size=m_ppara->get_pop_size();
			int num_dims=m_ppara->get_dim();
			double vtr=m_ppara->get_vtr();
			int max_gen=m_ppara->get_max_gen();
			int f_per_dim=m_ppara->get_f_per_dim();
			int f_dims=(f_per_dim ? num_dims : 1);
			int ini_f_type=m_ppara->get_ini_f_type();
			double ini_f_uni_low_bnd=m_ppara->get_ini_f_uni_low_bnd();
			double ini_f_uni_up_bnd=m_ppara->get_ini_f_uni_up_bnd();
			double ini_f_mean,ini_f_sigma;
			double pr_mean,pr_sigma;
			pr_mean=m_ppara->get_pr_mean();
			pr_sigma=m_ppara->get_pr_sigma();
			ini_f_mean=m_ppara->get_ini_f_mean();
			ini_f_sigma=m_ppara->get_ini_f_sigma();
			int out_interval=m_ppara->get_out_interval();
			int learn_p=m_ppara->get_learn_period();
			int pr_stra=m_ppara->get_pr_stra();

			int m_cur_run;
			int max_run=m_ppara->get_max_run();// run/trial number
			// shared_ptr<progress_display> pprog_dis;// algorithm progress indicator from boost
			// alloc_prog_indicator(pprog_dis);

			// allocate original pop and trial pop
			population pop(pop_size);
			allocate_pop(pop,num_dims,stra_num);
			population trial_pop(pop_size);
			allocate_pop(trial_pop,num_dims,stra_num);
			reallocate_stra(pop,trial_pop,f,f_dims);

			// generate algorithm statistics output file name
			ofstream stat_file(m_com_out_path.stat_path);
			// allocate stop condition object dynamically
			bool run_once=(1==max_run);
			alloc_stop_cond();

			int shuffle_size=pop_size-1;
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
				m_succ_f.clear();
				m_succ_cr.clear();
				// nsde SPECIFIC,initialize F value vector EVERY single run
				int i;
				for ( i=0;i<pop_size;i++ )
				{
					m_succ_f.clear();
					m_succ_cr.clear();

					uniform_01<> dist_01;
					variate_generator<mt19937&, uniform_01<> > rnd_01(gen, dist_01);
					pop[i].stra[pr].assign(1,rnd_01());
					// pop[i].stra[pr].assign(1,0.0);// Pr initial value: 0.0

					m_bi_norm_var.mu_x=ini_f_mean;
					m_bi_norm_var.sig_x=ini_f_sigma;
					if ( norm==ini_f_type )
						pop[i].stra[f][0]=::gen_rnd_norm(ini_f_mean,ini_f_sigma);
					else if ( uni==ini_f_type )
					{
						uniform_real<> dist_uni(ini_f_uni_low_bnd,ini_f_uni_up_bnd);
						variate_generator<mt19937&, uniform_real<> > rnd_uni(gen, dist_uni);
						pop[i].stra[f][0]=rnd_uni();
					}
				}
				set_orig_pop(pop);
				update_diversity(pop);

				calc_de_para_stat(pop);
				record_de_para_stat(m_cur_run);
				record_gen_vals(m_alg_stat,m_cur_run);
				// nsde SPECIFIC
				m_bi_norm_var.mu_y=m_de_stat.pr_mean;
				m_bi_norm_var.sig_y=m_de_stat.pr_std;
				m_bi_norm_var.rho=0.0;// initial correlation coefficient 

				print_run_times(stat_file,m_cur_run+1);
				print_run_title(stat_file);
				// output original population statistics
				print_gen_stat(stat_file,1,m_alg_stat);

				m_cur_gen=1;
				while ( false==(*m_pstop_cond) ) // for every iteration
				{
					m_de_stat.reset();
					int rnd_dim;
					double f_i,cr_i;
					double dim_mut_chance;
					int i,j,k;

					trial_pop=pop;// operator =
					//  recalculate succ_f_mean and succ_f_sigma periodically
					if ( is_learn_gen(m_cur_gen,learn_p) )
					{
						calc_bi_norm_var(m_succ_f,m_succ_cr,m_bi_norm_var);
						m_succ_f.clear();	
						m_succ_cr.clear();
					}// if ( is_learn_gen(m_cur_gen,learn_p) )
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

						gen_bi_norm_num(m_bi_norm_var,f_i,cr_i);
						if ( cr_i < 0.0 ) 
						{
							//double fric_part=ceil(f_i)-f_i;
							//f_i=1-fric_part;// assure positive
							cr_i=ceil(cr_i)-cr_i;;// assure positive
						}
						if ( cr_i > 1.0 )
							cr_i -= floor(cr_i);// truncate integral part//f_i=1.0;
						trial_pop[i].stra[pr][0]=cr_i;
						trial_pop[i].stra[f][0]=f_i;
						// scaling factor F self-adaptive update equation
						////// F->(0,1],repair F value if outbound
						//if ( f_i < 0.0 ) 
						//{
						//	//double fric_part=ceil(f_i)-f_i;
						//	//f_i=1-fric_part;// assure positive
						//	f_i=ceil(f_i)-f_i;;// assure positive
						//}
						//if ( f_i > 1.0 )
						//	f_i -= floor(f_i);// truncate integral part//f_i=1.0;
						for ( j=0;j<num_dims;j++ )
						{
							dim_mut_chance=rnd_01();
							if ( rnd_dim==j || dim_mut_chance<=cr_i )
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
					// initialize succeeded f and cr value repository
					if ( 1==m_cur_gen )
					{
						for ( i=0;i<pop_size;i++ )
						{
							m_succ_f.clear();
							m_succ_f.push_back(pop[i].stra[f][0]);
							m_succ_cr.clear();
							m_succ_cr.push_back(pop[i].stra[pr][0]);
						}
					}
					if ( is_output_gen(m_cur_gen,out_interval) )
					{
						stat_pop(pop,m_alg_stat);
						update_search_radius();
						update_diversity(pop);

						calc_de_para_stat(pop);
						record_de_para_stat(m_cur_run);
						record_gen_vals(m_alg_stat,m_cur_run);

						print_gen_stat(stat_file,m_cur_gen+1,m_alg_stat);
					}
					if ( is_learn_gen(m_cur_gen,learn_p) )
						record_gen_bi_norm_stat(m_bi_norm_var,m_cur_gen,m_cur_run);

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

	}// end namespace nsde
}// end namespace de
