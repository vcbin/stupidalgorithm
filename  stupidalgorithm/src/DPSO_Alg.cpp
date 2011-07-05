#include <boost/random/uniform_int.hpp>

#include "dpso_alg.h"

using std::vector;
using std::cout;
using std::string;
using std::ofstream;
using std::endl;
using std::ostream;
using std::ios;

using boost::shared_ptr;
// using boost::progress_display;
using boost::timer;

using boost::uniform_01;
using boost::uniform_int;
using boost::mt19937;
using boost::variate_generator;

extern boost::mt19937 gen;

namespace pso
{
	namespace dpso
	{
		enum {a=0,b};

		void dpso_stat::alloc_space(size_t num_dims)
		{
			gbest_inds.resize(pop_count);
			cen_inds.resize(pop_count);
			allocate_ind(gbest_inds[a],num_dims);
			allocate_ind(gbest_inds[b],num_dims);
			gbest_idx.resize(pop_count);
			avgs.resize(pop_count,0.0);
			stds.resize(pop_count,0.0);
			pos_divers.resize(pop_count,0.0);
			vel_divers.resize(pop_count,0.0);
			v_dim_divers.resize(pop_count);
			size_t i;
			for ( i=0;i<pop_count;i++ )
				v_dim_divers[i].assign(num_dims,0.0);
		}

		void dpso_stat::reset()
		{
			gbest_inds[a].obj[0]=gbest_inds[b].obj[0]=INF;
			gbest_idx[a]=gbest_idx[b]=-1;
			avgs.assign(avgs.size(),0.0);
			stds.assign(stds.size(),0.0);
		}

		void dpso_alg::print_run_title(ostream &os)
		{
			output_base::print_run_title(os);
			os<<"\t"
				<<"delta_of_avg"
				<<"\t"
				<<"("
				<<"gbest_pop_a"
				<<" "
				<<"gbest_pop_b"
				<<")"
				<<"\t"
				<<"pos_diver_of_pop_a"
				<<"\t"
				<<"pos_diver_of_pop_b"
				<<"\t"
				<<"vel_diver_of_pop_a"
				<<"\t"
				<<"vel_diver_of_pop_b"
				<<"std_of_pop_a"
				<<"\t"
				<<"std_of_pop_b";
		}

		void dpso_alg::print_gen_stat(ostream &os,int cur_gen,const alg_stat &alg_stat)
		{
			pso_alg::print_gen_stat(os,cur_gen,alg_stat);
			os.precision(8);// set float-point precision for output
			os.setf(ios::scientific,ios::floatfield);
			os<<"\t"
				<<"\t"
				<<m_dpso_stat.delta_avg
				<<"\t"
				<<"("
				<<m_dpso_stat.gbest_inds[a]
				<<" "
				<<m_dpso_stat.gbest_inds[b]
				<<")";

			os.precision(4);
			os.setf(ios::fixed,ios::floatfield);
			os<<"\t"
				<<m_dpso_stat.pos_divers[a]
				<<"\t"
				<<m_dpso_stat.pos_divers[b]
				<<"\t"
				<<m_dpso_stat.vel_divers[a]
				<<"\t"
				<<m_dpso_stat.vel_divers[b];

			os.setf(ios::scientific,ios::floatfield);
			os<<"\t"
				<<m_dpso_stat.stds[a]
				<<"\t"
				<<m_dpso_stat.stds[b];
		}

		void dpso_alg::initialize()
		{
			pso_alg::initialize();
			// initialize particle velocity 2D vector
			size_t pop_size=m_ppara->get_pop_size();
			int num_dims=m_ppara->get_dim();

			m_diff_pool.resize(pop_count);
			m_moment.resize(pop_size);
			m_dpso_stat.alloc_space(num_dims);
		}

		void dpso_alg::update_speed(population &pop)
		{
			double CR1,CR2;
			CR1=m_ppara->get_cr1();
			CR2=m_ppara->get_cr2();
			size_t num_dims=m_ppara->get_dim();
			// random 01 generator
			uniform_01<> dist;
			variate_generator<mt19937&, uniform_01<> > rnd_num(gen, dist);
			double r1,r2;

			size_t pop_size=pop.size();
			size_t subpop_size=pop_size/2;
			size_t i,j;
			for ( i=0;i<pop_size;i++ )
			{
				for ( j=0;j<num_dims;j++ )
				{
					if ( !m_outbound[i][j] ) // if particle didn't outbound
					{
						r1=rnd_num();
						r2=rnd_num();
						if ( i<subpop_size )
						{
							m_alg_stat.delta_x[i][j]= m_omega*m_alg_stat.delta_x[i][j] + 
								r1*CR1*(m_alg_stat.pbest[i].x[j]-pop[i].x[j]) +
								r2*CR2*(m_dpso_stat.gbest_inds[a].x[j]-pop[i].x[j]);
						}
						else
						{
							m_alg_stat.delta_x[i][j]= m_omega*m_alg_stat.delta_x[i][j] + 
								r1*CR1*(m_alg_stat.pbest[i].x[j]-pop[i].x[j]) +
								r2*CR2*(m_dpso_stat.gbest_inds[b].x[j]-pop[i].x[j]);
						}
						// boundaries check
						if ( abs(m_alg_stat.delta_x[i][j]) > m_vmax[j] )
						{
							if ( m_alg_stat.delta_x[i][j] > 0 )
								m_alg_stat.delta_x[i][j]=m_vmax[j];
							else
								m_alg_stat.delta_x[i][j]=-1.0*m_vmax[j];
						}
					}
					else
					{
						m_pBH->handle(m_alg_stat.delta_x[i][j]);
						m_outbound[i][j]=false;
					}
				}// for every dimenson
			}// for every particle
		}// end function update_speed

		void dpso_alg::stat_ini_pop( population &pop,alg_stat &alg_stat )
		{
			int num_dims=m_ppara->get_dim();

			alg_stat.gbest.obj[0]=INF;
			// dpso algo specific
			m_dpso_stat.gbest_inds[a].obj[0]=INF;
			m_dpso_stat.gbest_inds[b].obj[0]=INF;

			m_alg_stat.avg=0.0;
			m_alg_stat.cen_ind.assign(num_dims,0.0);
			// dpso specific
			m_dpso_stat.avgs[a]=m_dpso_stat.avgs[b]=0.0;
			m_dpso_stat.cen_inds[a].assign(num_dims,0.0);
			m_dpso_stat.cen_inds[b].assign(num_dims,0.0);

			size_t pop_size=pop.size();
			size_t subpop_size=pop_size/2;
			bool stag_flag=true;
			size_t i;
			for ( i=0;i<pop_size;i++ )
			{
				// stat the pop step by step
				m_alg_stat.cen_ind += pop[i].x;
				m_alg_stat.avg += pop[i].obj[0];

				if ( pop[i] < alg_stat.gbest )
				{
					alg_stat.gbest=pop[i];
					stag_flag=false;
					alg_stat.num_stag=0;
				}

				// record position of gbest_a,gbest_b and avgerage fitness of pop_a,pop_b,
				// dpso specific
				if ( i<subpop_size )
				{
					m_dpso_stat.avgs[a] += pop[i].obj[0];
					if ( pop[i] < m_dpso_stat.gbest_inds[a] )
					{
						m_dpso_stat.gbest_inds[a]=pop[i];
						m_dpso_stat.gbest_idx[a]=i;
					}
				}
				else
				{
					if ( i==subpop_size )
					{
						m_dpso_stat.cen_inds[a]=alg_stat.cen_ind;
						m_dpso_stat.avgs[a] /= (1.0*subpop_size);
					}
					m_dpso_stat.cen_inds[b] += pop[i].x;
					m_dpso_stat.avgs[b] += pop[i].obj[0];

					if ( pop[i] < m_dpso_stat.gbest_inds[b] )
					{
						m_dpso_stat.gbest_inds[b]=pop[i];
						m_dpso_stat.gbest_idx[b]=i;
					}
				}
			}// for every particle
			// end evaluation
			m_dpso_stat.cen_inds[a] /= (subpop_size*1.0);
			m_dpso_stat.cen_inds[b] /= (subpop_size*1.0);

			m_alg_stat.cen_ind /=  pop_size;// update centroid
			m_alg_stat.avg /=  pop_size;// average fitness value

			m_dpso_stat.avgs[b] /= (subpop_size*1.0);
			m_dpso_stat.delta_avg=m_dpso_stat.avgs[a]-m_dpso_stat.avgs[b];

			if (stag_flag)
				alg_stat.num_stag++;// increment stagnation count

			// calclate standard deviation
			m_alg_stat.std=0.0;
			m_dpso_stat.stds[a]=m_dpso_stat.stds[b]=0.0;
			for ( i=0;i<pop_size;i++ )
			{
				double diff=pop[i].obj[0]-m_alg_stat.avg;
				m_alg_stat.std += diff*diff;
				if ( i<subpop_size )
				{
					diff=pop[i].obj[0]-m_dpso_stat.avgs[a];
					m_dpso_stat.stds[a] += diff*diff;
				}
				else
				{
					diff=pop[i].obj[0]-m_dpso_stat.avgs[b];
					m_dpso_stat.stds[b] += diff*diff;
				}
			}
			m_alg_stat.std/=(pop_size*1.0);
			m_alg_stat.std=sqrt(m_alg_stat.std);
			m_dpso_stat.stds[a]/=(subpop_size*1.0);
			m_dpso_stat.stds[a]=sqrt(m_dpso_stat.stds[a]);
			m_dpso_stat.stds[b]/=(subpop_size*1.0);
			m_dpso_stat.stds[b]=sqrt(m_dpso_stat.stds[b]);

			alg_stat.pbest=pop;
		}// end function stat_ini_pop

		void dpso_alg::stat_pop(population &pop,alg_stat &alg_stat)
		{
			int num_dims=m_ppara->get_dim();

			size_t pop_size=pop.size();
			size_t subpop_size=pop_size/2;

			m_alg_stat.avg=0.0;
			m_alg_stat.cen_ind.assign(m_alg_stat.cen_ind.size(),0.0);

			// dpso specific
			m_dpso_stat.avgs[a]=m_dpso_stat.avgs[b]=0.0;
			m_dpso_stat.cen_inds[a].assign(num_dims,0.0);
			m_dpso_stat.cen_inds[b].assign(num_dims,0.0);

			bool stag_flag=true;
			unsigned i;
			for ( i=0;i<pop_size;i++ )
			{
				// stat the pop step by step
				m_alg_stat.cen_ind = m_alg_stat.cen_ind+pop[i].x;
				m_alg_stat.avg += pop[i].obj[0];

				if ( pop[i] < m_alg_stat.pbest[i] )
					m_alg_stat.pbest[i]=pop[i];// operator=
				if ( pop[i] < m_alg_stat.gbest )
				{
					m_alg_stat.gbest=pop[i];
					stag_flag=false;
					m_alg_stat.num_stag=0;
				}
				// dpso algo specific
				if ( i<subpop_size )
				{
					// update gbest of pop_a
					// record position of gbest_a,gbest_b
					// dpso specific
					if ( pop[i] < m_dpso_stat.gbest_inds[a] )
					{
						m_dpso_stat.gbest_inds[a]=pop[i];
						m_dpso_stat.gbest_idx[a]=i;
					}
				}
				else
				{
					if ( i==subpop_size )
					{
						m_dpso_stat.cen_inds[a]=alg_stat.cen_ind;
						m_dpso_stat.avgs[a]=m_alg_stat.avg/(1.0*subpop_size);
					}
					m_dpso_stat.cen_inds[b] += pop[i].x;
					m_dpso_stat.avgs[b] += pop[i].obj[0];
					// update gbest of pop_b
					if ( pop[i] < m_dpso_stat.gbest_inds[b] )
					{
						m_dpso_stat.gbest_inds[b]=pop[i];
						m_dpso_stat.gbest_idx[b]=i;
					}
				}
			}// for every particle

			if (stag_flag)
				m_alg_stat.num_stag++;// increment stagnation count
			m_dpso_stat.cen_inds[a] /= (subpop_size*1.0);
			m_dpso_stat.cen_inds[b] /= (subpop_size*1.0);
			m_alg_stat.cen_ind /=  pop_size;// update centroid
			m_alg_stat.avg /=  pop_size;// average fitness value
			m_dpso_stat.avgs[b] /= (subpop_size*1.0);
			m_dpso_stat.delta_avg=m_dpso_stat.avgs[a]-m_dpso_stat.avgs[b];

			// calclate standard deviation
			m_alg_stat.std=0.0;
			m_dpso_stat.stds[a]=m_dpso_stat.stds[b]=0.0;
			for ( i=0;i<pop_size;i++ )
			{
				double diff;
				diff=pop[i].obj[0]-m_alg_stat.avg;
				m_alg_stat.std += diff*diff;
				if ( i<subpop_size )
				{
					diff=pop[i].obj[0]-m_dpso_stat.avgs[a];
					m_dpso_stat.stds[a] += diff*diff;
				}
				else
				{
					diff=pop[i].obj[0]-m_dpso_stat.avgs[b];
					m_dpso_stat.stds[b] += diff*diff;
				}
			}
			m_alg_stat.std/=(pop_size*1.0);
			m_alg_stat.std=sqrt(m_alg_stat.std);
			m_dpso_stat.stds[a]/=(subpop_size*1.0);
			m_dpso_stat.stds[a]=sqrt(m_dpso_stat.stds[a]);
			m_dpso_stat.stds[b]/=(subpop_size*1.0);
			m_dpso_stat.stds[b]=sqrt(m_dpso_stat.stds[b]);
		}// end function stat_pop

		void dpso_alg::update_diversity(population &pop)
		{
			size_t pop_size=m_ppara->get_pop_size();
			size_t subpop_size=pop_size/pop_count;
			int num_dims=m_ppara->get_dim();

			size_t i;
			int j;
			m_alg_stat.pos_diver=0.0;
			for ( i=0;i<pop_size;i++ )
			{
				m_alg_stat.pos_diver += vec_distance(pop[i].x,m_alg_stat.cen_ind);
				if ( i>=subpop_size )
				{
					if ( i==subpop_size )
						m_dpso_stat.pos_divers[a]=m_alg_stat.pos_diver;
					m_dpso_stat.pos_divers[b] += vec_distance(pop[i].x,m_dpso_stat.cen_inds[b]);
				}
			}// for every particle
			double len=m_diag_len*pop_size;
			m_alg_stat.pos_diver /= len;
			m_dpso_stat.pos_divers[a] /= len;
			m_dpso_stat.pos_divers[b] /= len;
			// calculate velocity diversity
			vector<double> vec_avg(num_dims,0.0);
			vector<double> vec_pop_a_avg(num_dims,0.0);
			vector<double> vec_pop_b_avg(num_dims,0.0);
			for ( j=0;j<num_dims;j++ )
			{
				for ( i=0;i<pop_size;i++ )
				{
					vec_avg[j] += m_alg_stat.delta_x[i][j];
					if ( i>=subpop_size )
					{
						if ( i==subpop_size )
							vec_pop_a_avg[j]=vec_avg[j];
						vec_pop_b_avg[j] += m_alg_stat.delta_x[i][j];
					}
				}// for every particle
				vec_avg[j] /=  pop_size;
				vec_pop_a_avg /= (subpop_size*1.0);
				vec_pop_b_avg /= (subpop_size*1.0);
			}// for every dimension
			m_v_dim_diver.assign(num_dims,0.0);
			m_dpso_stat.v_dim_divers[a].assign(num_dims,0.0);
			m_dpso_stat.v_dim_divers[b].assign(num_dims,0.0);
			for ( j=0;j<num_dims;j++ )
			{
				for ( i=0;i<pop_size;i++ )
				{
					m_v_dim_diver[j] += abs(m_alg_stat.delta_x[i][j]-vec_avg[j]);
					if ( i>=subpop_size )
					{
						if ( i==subpop_size )
							m_dpso_stat.v_dim_divers[a][j]=m_v_dim_diver[j];
						m_dpso_stat.v_dim_divers[b][j] += abs(m_alg_stat.delta_x[i][j]-vec_avg[j]);
					}
				}// for every particle
				m_v_dim_diver[j] /=  pop_size;
				m_dpso_stat.v_dim_divers[a][j] /= (subpop_size*1.0);
				m_dpso_stat.v_dim_divers[b][j] /= (subpop_size*1.0);
			}// for every dimension
			m_vel_diver=0.0;
			for ( j=0;j<num_dims;j++ )
			{
				m_vel_diver += m_v_dim_diver[j];
				m_dpso_stat.vel_divers[a] += m_dpso_stat.v_dim_divers[a][j];
				m_dpso_stat.vel_divers[b] += m_dpso_stat.v_dim_divers[b][j];
			}
			m_vel_diver /= (num_dims*1.0);
			m_dpso_stat.vel_divers[a] /= (num_dims*1.0);
			m_dpso_stat.vel_divers[b] /= (num_dims*1.0);
		}// end function update_diversity

		void dpso_alg::update_momentum(population &pop)
		{
			size_t pop_size=pop.size();
			size_t num_dims=m_ppara->get_dim();
			unsigned i,j;
			m_moment.assign(m_moment.size(),0.0);// reinitialize to 0
			for ( i=0;i<pop_size;i++ )
			{
				double diff;
				for ( j=0;j<num_dims;j++ )
				{
					diff=m_alg_stat.delta_x[i][j];
					m_moment[i] += diff*diff;
				}
				m_moment[i] *= 0.5;
			}
		}// end function update_momentum

		// update temperature(average momentum of all particle) of dpso
		void dpso_alg::update_temp()
		{
			size_t pop_size=m_ppara->get_pop_size();
			unsigned i;
			m_t=0.0;
			for ( i=0;i<pop_size;i++ )
				m_t += m_moment[i]; 
			m_t /=  pop_size;
		}

		void dpso_alg::update_pool(population &pop)
		{
			size_t pop_size=m_ppara->get_pop_size();
			size_t subpop_size=pop_size/2;
			unsigned i;
			double diff_prob=0.0;// diffusion probability
			uniform_01<> dist;
			variate_generator<mt19937&, uniform_01<> > rnd_num(gen, dist);

			m_diff_pool[a].clear();
			for ( i=0;i<subpop_size;i++ )
			{
				diff_prob = 1.0-1.0/exp(abs(m_moment[i]/m_t));
				if ( rnd_num() <= diff_prob )
					m_diff_pool[a].push_back(pop[i]);
			}

			m_diff_pool[b].clear();
			for ( i=subpop_size;i<pop_size;i++ )
			{
				diff_prob = 1.0-1.0/exp(abs(m_moment[i]/m_t));
				if ( rnd_num() <= diff_prob )
					m_diff_pool[b].push_back(pop[i]);
			}

		}// end function update_pool

		void dpso_alg::diffuse(population &pop,func_base &func_eval)
		{
			size_t pool_size;
			// pool a
			pool_size=m_diff_pool[a].size();
			if ( pool_size>= 2 )
			{
				// generate two different index
				uniform_int<> dist(0,pool_size-1);
				variate_generator<mt19937&, uniform_int<> > rnd_num(gen, dist);
				int rnd_1,rnd_2;
				rnd_1=rnd_num();
				while ( (rnd_2=rnd_num()) == rnd_1 )
					;

				individual trial_ind;
				size_t num_dims=m_ppara->get_dim();
				allocate_ind(trial_ind,num_dims);

				trial_ind.x = m_dpso_stat.gbest_inds[a].x+( m_diff_pool[a][rnd_1].x-m_diff_pool[a][rnd_2].x );
				func_eval(trial_ind);
				m_alg_stat.all_eval_num++;
				m_alg_stat.eval_num++;
				if ( trial_ind < m_dpso_stat.gbest_inds[b] )// better than gbest of population b
				{
					pop[m_dpso_stat.gbest_idx[b]]=trial_ind;
					m_dpso_stat.gbest_inds[b]=trial_ind;
					m_dpso_stat.all_exch_num++;
				}
			}
			// pool b
			pool_size=m_diff_pool[b].size();
			if ( pool_size>= 2 )
			{
				// generate two different index
				uniform_int<> dist(0,pool_size-1);
				variate_generator<mt19937&, uniform_int<> > rnd_num(gen, dist);
				int rnd_1,rnd_2;
				rnd_1=rnd_num();
				while ( (rnd_2=rnd_num()) == rnd_1 )
					;

				individual trial_ind;
				size_t num_dims=m_ppara->get_dim();
				allocate_ind(trial_ind,num_dims);
				trial_ind.x = m_dpso_stat.gbest_inds[b].x+( m_diff_pool[b][rnd_1].x-m_diff_pool[b][rnd_2].x );
				func_eval(trial_ind);
				m_alg_stat.all_eval_num++;
				m_alg_stat.eval_num++;

				if ( trial_ind < m_dpso_stat.gbest_inds[a] )// better than gbest of population a
				{
					pop[m_dpso_stat.gbest_idx[a]]=trial_ind;
					m_dpso_stat.gbest_inds[a]=trial_ind;
					m_dpso_stat.all_exch_num++;
				}
			}
		}// end function diffuse

		inline void dpso_alg::print_Exch_stat(ostream &os)
		{
			os<<"\n"
				<<"total_exchange_num:"
				<<m_dpso_stat.all_exch_num
				<<"\t"
				<<"exchange_per_run:"
				<<m_dpso_stat.all_exch_num/(m_ppara->get_max_run()*1.0);
		}

		int dpso_alg::run()
		{
			if ( !m_ppara )
				return -1;

			timer elapsed_t;

			// retrieve algorithm parameters
			double vtr=m_ppara->get_vtr();
			size_t pop_size=m_ppara->get_pop_size();
			int num_dims=m_ppara->get_dim();
			int vmax_type=m_ppara->get_vmax_type();

			// allocate original pop
			population pop(pop_size);
			allocate_pop(pop,num_dims);
			m_cur_gen=1;

			int m_cur_run;
			int max_run=m_ppara->get_max_run();// run/trial number
			//shared_ptr<progress_display> pprog_dis;// algorithm progress indicator from boost
			// alloc_prog_indicator(pprog_dis);

			// algorithm statistics output file
			ofstream stat_file(m_com_out_path.stat_path.c_str());
			// allocate stop condition object dynamically
			alloc_stop_cond();
			int ob_type=m_ppara->get_ob_type();
			alloc_bound_handle(ob_type,m_pBH);

			// iteration start
			for ( m_cur_run=0;m_cur_run<max_run;m_cur_run++ )
			{
				reset_run_stat();

				set_orig_pop(pop);
				set_ini_max_velocity();
				update_diversity(pop);

				record_gen_vals(m_alg_stat,m_cur_run);

				print_run_times(stat_file,m_cur_run+1);
				print_run_title(stat_file);
				// output original population statistics
				print_gen_stat(stat_file,1,m_alg_stat);

				m_cur_gen=1;
				while ( false==(*m_pstop_cond) ) // for every generation
				{
					update_omega(m_cur_gen);
					update_speed(pop);
					update_search_radius();
					update_pop(pop);
					if ( dynamic_rng==vmax_type )
					{
						update_x_rng(pop);
						set_dynamic_vmax();
					}

					eval_pop(pop,*m_pfunc,m_alg_stat);
					stat_pop(pop,m_alg_stat);
					update_diversity(pop);

					record_gen_vals(m_alg_stat,m_cur_run);

					// dpso specific
					update_momentum(pop);
					update_temp();
					update_pool(pop);
					diffuse(pop,*m_pfunc);

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
			print_Exch_stat(stat_file);// dpso specific
			m_alg_stat.run_avg_time=elapsed_t.elapsed();
			m_alg_stat.run_avg_time/=(max_run*1.0);
			print_avg_time(stat_file,m_alg_stat.run_avg_time);
			
			print_best_x(stat_file,m_alg_stat.bst_ind);
			write_stat_vals();

			cout<<endl;// flush cout output
			return 0;
		}// end function Run

	}// namespace dpso
}// namespace pso
