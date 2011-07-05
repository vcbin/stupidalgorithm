#include "my_sde_alg.h"

#include "sde_alg.h"
#include "rand_val.h"

using std::vector;
using std::cout;
using std::string;
using std::ofstream;
using std::ostream;
using std::swap;
using std::endl;

using boost::shared_ptr;
// using boost::progress_display;
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
	namespace sde
	{
		void my_sde_alg::initialize()
		{
			bbde_alg::initialize();

			int pop_size=m_ppara->get_pop_size();
			m_succ_pr.resize(pop_size);
			int num_dims=m_ppara->get_dim();
			int f_per_dim=m_ppara->get_f_per_dim();
			int f_dims=(f_per_dim ? num_dims : 1);
			m_succ_f.resize(f_dims);
		}// end function initialize

		inline bool my_sde_alg::is_learn_gen(int cur_gen,int learn_p)
		{
			return ( 0==(cur_gen+1)%learn_p );
		}

		void my_sde_alg::update_pop(population &pop,const population &trial_pop)
		{
			int pop_size=pop.size();
			int num_dims=m_ppara->get_dim();
			int f_per_dim=m_ppara->get_f_per_dim();
			int f_dims=(f_per_dim ? num_dims : 1);
			int i,j;
			for ( i=0;i<pop_size;i++ )
			{
				if ( trial_pop[i]<pop[i] )
				{
					m_alg_stat.delta_x[i]=trial_pop[i].x-pop[i].x;// update step size delta_x
					pop[i]=trial_pop[i];
					for ( j=0;j<f_dims;j++ )
						m_succ_f[j].push_back(pop[i].stra[f][j]);// record succeeded F values
				}
				else
					m_alg_stat.delta_x[i].assign(num_dims,0.0);// trial failed,step size delta_x == 0
			}// for every individual
		}// end function update_pop

		int my_sde_alg::run()
		{
			if ( !m_ppara )
				return -1;

			timer elapsed_t;
			// retrieve algorithm parameters
			int pop_size=m_ppara->get_pop_size();
			int num_dims=m_ppara->get_dim();
			double vtr=m_ppara->get_vtr();
			int f_per_dim=m_ppara->get_f_per_dim();
			int f_dims=(f_per_dim ? num_dims : 1);
			double pr_mean,pr_sigma;
			int ini_f_type=m_ppara->get_ini_f_type();
			double ini_f_uni_low_bnd=m_ppara->get_ini_f_uni_low_bnd();
			double ini_f_uni_up_bnd=m_ppara->get_ini_f_uni_up_bnd();
			double ini_f_mean,ini_f_sigma;
			pr_mean=m_ppara->get_pr_mean();
			pr_sigma=m_ppara->get_pr_sigma();
			ini_f_mean=m_ppara->get_ini_f_mean();
			ini_f_sigma=m_ppara->get_ini_f_sigma();
			d_array succ_f_mean,succ_f_sigma;
			succ_f_mean.resize(f_dims);
			succ_f_sigma.resize(f_dims);
            int learn_p=m_ppara->get_learn_period();
			int pr_stra=m_ppara->get_pr_stra();
			bool succ_pr_empty=true;
			double succ_pr_mean,succ_pr_sigma;
			int out_interval=m_ppara->get_out_interval();

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
			ofstream stat_file(m_com_out_path.stat_path.c_str());
			// allocate stop condition object dynamically
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
				m_succ_pr.clear();
				succ_pr_empty=true;
                // sde SPECIFIC,initialize F value vector EVERY single run
               int i,j;
                for ( i=0;i<pop_size;i++ )
                {
					uniform_01<> dist_01;
					variate_generator<mt19937&, uniform_01<> > rnd_01(gen, dist_01);
					pop[i].stra[pr].assign(1,rnd_01());
					// pop[i].stra[pr].assign(1,0.0);// Pr initial value: 0.0
					for ( j=0;j<f_dims;j++ )
					{
						m_succ_f[j].clear();
						succ_f_mean[j]=ini_f_mean;
						succ_f_sigma[j]=ini_f_sigma;
						if ( norm==ini_f_type )
							pop[i].stra[f][j]=::gen_rnd_norm(ini_f_mean,ini_f_sigma);
						else if ( uni==ini_f_type )
						{
							uniform_real<> dist_uni(ini_f_uni_low_bnd,ini_f_uni_up_bnd);
							variate_generator<mt19937&, uniform_real<> > rnd_uni(gen, dist_uni);
							pop[i].stra[f][j]=rnd_uni();
						}
					}// for ( j=0;j<f_dims;j++ )
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
                    double rnd_norm_Pr;
                    double dim_mut_chance;
                    int i,j,k;

                    trial_pop=pop;// operator =
                    //  recalculate succ_f_mean and succ_f_sigma periodically
                    if ( is_learn_gen(m_cur_gen,learn_p) )
                    {
						for ( j=0;j<f_dims;j++ )
						{
							int f_size=m_succ_f[j].size();
							if ( f_size ) 
							{
								succ_f_mean[j]=0.0;
								int m;
								for ( m=0;m<f_size;m++ )
									succ_f_mean[j] += m_succ_f[j][m];
								succ_f_mean /= f_size;
								succ_f_sigma[j]=0.0;
								double diff;
								for ( m=0;m<f_size;m++ )
								{
									diff=m_succ_f[j][m]-succ_f_mean[j];
									succ_f_sigma[j] += diff*diff;
								}
								succ_f_sigma[j] /= f_size;
								succ_f_sigma[j] = sqrt(succ_f_sigma[j]);
							}
							else
							{
								succ_f_mean[j]=ini_f_mean;
								succ_f_sigma[j]=ini_f_sigma;
							}
						}
						for ( j=0;j<f_dims;j++ )
							m_succ_f[j].clear();
						for ( i=0;i<pop_size;i++ )
						{
							int pr_size=m_succ_pr.size();
							if ( pr_size ) 
							{
								succ_pr_empty=false;
								succ_pr_mean=0.0;
								int m;
								for ( m=0;m<pr_size;m++ )
									succ_pr_mean += m_succ_pr[m];
								succ_pr_mean /= pr_size;
								succ_pr_sigma=0.0;
								double diff;
								for ( m=0;m<pr_size;m++ )
								{
									diff=m_succ_pr[m]-succ_pr_mean;
									succ_pr_sigma += diff*diff;
								}
								succ_pr_sigma /= pr_size;
								succ_pr_sigma = sqrt(succ_pr_sigma);
							}
							else
								succ_pr_empty=true;
						}
						m_succ_pr.clear();
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
						if ( stra_norm==pr_stra )
						{
							rnd_norm_Pr=generate_rnd_pr(pr_mean,pr_sigma);// generate random N(mean,sigma) crossover probability
							trial_pop[i].stra[pr][0]=rnd_norm_Pr;
						}
						else
						{
							if ( !succ_pr_empty )
							{
								rnd_norm_Pr=generate_rnd_pr(succ_pr_mean,succ_pr_sigma);// generate random N(mean,sigma) crossover probability
								trial_pop[i].stra[pr][0]=rnd_norm_Pr;
							}
							else
							{
								uniform_01<> dist_01;
								variate_generator<mt19937&, uniform_01<> > rnd_01(gen, dist_01);
								trial_pop[i].stra[pr][0]=rnd_01();
							}
						}

						// scaling factor F self-adaptive update equation
						f_i=::gen_rnd_norm(succ_f_mean[0],succ_f_sigma[0]);
						////// F->(0,1],repair F value if outbound
						//if ( f_i < 0.0 ) 
						//{
						//	//double fric_part=ceil(f_i)-f_i;
						//	//f_i=1-fric_part;// assure positive
						//	f_i=ceil(f_i)-f_i;;// assure positive
						//}
						//if ( f_i > 1.0 )
						//	f_i -= floor(f_i);// truncate integral part//f_i=1.0;
						trial_pop[i].stra[f][0]=f_i;

                        for ( j=0;j<num_dims;j++ )
						{
							if ( f_per_dim && j>0 )
							{
								f_i=::gen_rnd_norm(succ_f_mean[j],succ_f_sigma[j]);
								////// F->(0,1],repair F value if outbound
								//if ( f_i < 0.0 ) 
								//{
								//	//double fric_part=ceil(f_i)-f_i;
								//	//f_i=1-fric_part;// assure positive
								//	f_i=ceil(f_i)-f_i;;// assure positive
								//}
								//if ( f_i > 1.0 )
								//	f_i -= floor(f_i);// truncate integral part//f_i=1.0;
								pop[i].stra[f][j]=f_i;
							}

							dim_mut_chance=rnd_01();
							if ( rnd_dim==j || dim_mut_chance<=trial_pop[i].stra[pr][0] )
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
							for ( j=0;j<f_dims;j++ )
								m_succ_f[j].push_back(pop[i].stra[f][j]);
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

	}// end namespace sde
}// end namespace de
