#include <boost/random/uniform_int.hpp>

#include "PSOBC_Alg.h"

using boost::timer;
using boost::progress_display;
using boost::shared_ptr;

using std::cout;
using std::endl;
using std::ofstream;
using std::vector;
using std::numeric_limits;

using boost::uniform_01;
using boost::mt19937;
using boost::variate_generator;

extern boost::mt19937 gen;

namespace pso
{
	namespace psobc
	{

		void psobc_stat::alloc_spaces(size_t pop_size,int num_dims)
		{
			pre_pop_val.resize(pop_size);

			pworst.resize(pop_size);
			allocate_pop(pworst,num_dims);
			allocate_ind(gworst,num_dims);
		}

		bool operator>(individual &lhs,double val) { return lhs.obj[0]>val; }

		int copy_obj_val_from_pop(vector<double> &pop_val,const population &pop)
		{
			size_t pop_size=pop.size();
			if ( pop_val.size()!= pop_size )
				return -1;
			size_t i;
			for ( i=0;i<pop_size;i++ )
				pop_val[i]=pop[i].obj[0];
			return 0;
		}

		void psobc_alg::initialize()
		{
			size_t pop_size=m_ppara->get_pop_size();
			int num_dims=m_ppara->get_dim();
			arpso_alg::initialize();
			m_psobc_stat.alloc_spaces(pop_size,num_dims);

			m_prev_dir.resize(pop_size);
			m_ob_handled.resize(pop_size);
			size_t i;
			for ( i=0;i<pop_size;i++ )
			{
				m_prev_dir[i].assign(num_dims,0.0);
				m_ob_handled[i].assign(num_dims,false);
			}
		}

		void psobc_alg::stat_pop(population &pop,alg_stat &alg_stat)
		{
			size_t pop_size=pop.size();
			unsigned i;

			m_alg_stat.avg=0.0;
			m_alg_stat.cen_ind.assign(m_alg_stat.cen_ind.size(),0.0);
			bool stag_flag=true;
			for ( i=0;i<pop_size;i++ )
			{
				// stat the pop step by step
				m_alg_stat.cen_ind = m_alg_stat.cen_ind+pop[i].x;
				m_alg_stat.avg += pop[i].obj[0];
				m_alg_stat.eval_num++;
				m_alg_stat.all_eval_num++;

				if ( pop[i] < m_alg_stat.pbest[i] )
					m_alg_stat.pbest[i]=pop[i];// operator=
				// psobc specific
				if ( pop[i] > m_psobc_stat.pworst[i] )
					m_psobc_stat.pworst[i]=pop[i];// operator=
				if ( pop[i] < m_alg_stat.gbest )
				{
					m_alg_stat.gbest=pop[i];
					/*if ( m_alg_stat.gbest.obj[0]<1.0E-6 )
					cout<<" ";*/
					stag_flag=false;
					m_alg_stat.num_stag=0;
				}
				// psobc specific
				if ( pop[i] > m_psobc_stat.gworst )
					m_psobc_stat.gworst=pop[i];
			}// for every particle
			//  end evaluation

			if (stag_flag)
				m_alg_stat.num_stag++;// increment stagnation count
			m_alg_stat.cen_ind=m_alg_stat.cen_ind/(pop_size*1.0);// update centroid
			m_alg_stat.avg /=  pop_size;// average fitness value

			// calclate standard deviation
			m_alg_stat.std=0.0;
			for ( i=0;i<pop_size;i++ )
			{
				double diff=pop[i].obj[0]-m_alg_stat.avg;
				m_alg_stat.std += diff*diff;
			}
			m_alg_stat.std/=(pop_size*1.0);
			m_alg_stat.std=sqrt(m_alg_stat.std);
		}// end function stat_pop

		void psobc_alg::stat_ini_pop( population &pop,alg_stat &alg_stat )
		{
			size_t pop_size=pop.size();
			// pbest,gbest first time initialization
			double d_max=INF;
			m_alg_stat.gbest.obj[0]=d_max;
			// psobc specific
			double d_min=numeric_limits<double>::lowest();
			m_psobc_stat.gworst.obj[0]=d_min;
			m_psobc_stat.pre_pop_val.assign(pop_size,d_min);

			m_alg_stat.avg=0.0;
			m_alg_stat.cen_ind.assign(m_alg_stat.cen_ind.size(),0.0);
			bool stag_flag=true;

			size_t subpop_size=pop_size/2;
			size_t num_dims=pop[0].x.size();
			size_t i;
			for ( i=0;i<pop_size;i++ )
			{
				// stat the pop step by step
				m_alg_stat.cen_ind = m_alg_stat.cen_ind+pop[i].x;
				m_alg_stat.avg += pop[i].obj[0];
				m_alg_stat.eval_num++;
				m_alg_stat.all_eval_num++;

				if ( pop[i] < alg_stat.gbest )
				{
					alg_stat.gbest=pop[i];
					stag_flag=false;
					alg_stat.num_stag=0;
				}

				// psobc specific
				if ( pop[i] > m_psobc_stat.gworst )
					m_psobc_stat.gworst=pop[i];
			}// for every particle
			// end evaluation

			if (stag_flag)
				m_alg_stat.num_stag++;// increment stagnation count
			m_alg_stat.cen_ind=m_alg_stat.cen_ind/(pop_size*1.0);// update centroid
			m_alg_stat.avg /=  pop_size;// average fitness value

			// calclate standard deviation
			m_alg_stat.std=0.0;
			for ( i=0;i<pop_size;i++ )
			{
				double diff=pop[i].obj[0]-m_alg_stat.avg;
				m_alg_stat.std += diff*diff;
			}
			m_alg_stat.std/=(pop_size*1.0);
			m_alg_stat.std=sqrt(m_alg_stat.std);

			alg_stat.pbest=pop;
			// psobc specific
			m_psobc_stat.pworst=pop;
		}// end function stat_ini_pop
	

		void psobc_alg::update_speed(population &pop)
		{
			double CR1,CR2;
			CR1=m_ppara->get_cr1();
			CR2=m_ppara->get_cr2();
			size_t num_dims=m_ppara->get_dim();
			// random 01 generator
			uniform_01<> dist;
			variate_generator<mt19937&, uniform_01<> > rnd_num(gen, dist);

			size_t pop_size=pop.size();
			size_t i,j;
			for ( i=0;i<pop_size;i++ )
			{
				for ( j=0;j<num_dims;j++ )
				{
					if ( !m_outbound[i][j] ) // if particle didn't outbound
					{
						double r1,r2;
						r1=rnd_num();
						r2=rnd_num();
						if ( m_dir_accl > 0 ) // attraction phase, m_dir_accl=1.0
						{
							if ( pop[i] > m_psobc_stat.pre_pop_val[i] ) // current particle is worse than previous generation
							{
								// hence update velocity
								m_alg_stat.delta_x[i][j]= m_omega*m_alg_stat.delta_x[i][j] + 
									m_dir_accl*(
									r1*CR1*(m_alg_stat.pbest[i].x[j]-pop[i].x[j]) +
									r2*CR2*(m_alg_stat.gbest.x[j]-pop[i].x[j])
									);
							}
						}
						else // repulsion phase, m_dir_accl= -1.0
						{
							m_alg_stat.delta_x[i][j]= m_omega*m_alg_stat.delta_x[i][j] + 
								m_dir_accl*(
								r1*CR1*(m_psobc_stat.pworst[i].x[j]-pop[i].x[j]) +
								r2*CR2*(m_psobc_stat.gworst.x[j]-pop[i].x[j])
								);
						}
						// boundaries check
						if ( abs(m_alg_stat.delta_x[i][j]) > m_vmax[j] )
						{
							if ( m_alg_stat.delta_x[i][j] > 0 )
								m_alg_stat.delta_x[i][j]=m_vmax[j];
							else
								m_alg_stat.delta_x[i][j]=-1.0*m_vmax[j];
						}

						// psobc specific:keep previous direction to escape from boundary steadily(NO bounce back for one period)
						if ( m_ob_handled[i][j] )
						{
							m_alg_stat.delta_x[i][j] *= m_prev_dir[i][j];
							m_ob_handled[i][j]=false;
						}
					}
					else
					{
						m_pBH->handle(m_alg_stat.delta_x[i][j]);
						m_outbound[i][j]=false;
						m_ob_handled[i][j]=true;
						m_prev_dir[i][j]=((m_alg_stat.delta_x[i][j]>0 ? 1.0:-1.0));// record current direction for steady escaping
					}
				}// for every dimension
			}// for every particle
		}// end function update_speed

		int psobc_alg::run()
		{
			if ( !m_ppara )
				return -1;

			timer elapsed_t;

			// retrieve algorithm parameters
			int max_gen=m_ppara->get_max_gen();
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
			bool run_once= 1==max_run;
			// alloc_prog_indicator(pprog_dis);

			// algorithm statistics output file
			ofstream stat_file(m_com_out_path.stat_path);

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
				while ( false==(*m_pstop_cond) )
				{
					update_omega(m_cur_gen);
					update_dir(pop);// ARPSO-likewise algo specific procedure
					stat_rep_gen(m_cur_gen);// ARPSO specific procedure for stat
					update_speed(pop);
					update_search_radius();
					update_pop(pop);
					if ( dynamic_rng==vmax_type )
					{
						update_x_rng(pop);
						set_dynamic_vmax();
					}
					// psobc specific:record previous population fitness values
					copy_obj_val_from_pop(m_psobc_stat.pre_pop_val,pop);

					eval_pop(pop,*m_pfunc,m_alg_stat);
					stat_pop(pop,m_alg_stat);
					update_diversity(pop);

					record_gen_vals(m_alg_stat,m_cur_run);

					print_gen_stat(stat_file,m_cur_gen+1,m_alg_stat);
					update_conv_stat(vtr);

					/*if ( run_once )
						++(*pprog_dis);*/

					m_cur_gen++;
				}// while single run stop_condition is false

				// single run end

				// arpso specific output
				calc_rep_gen_perc();
				print_rep_stat(stat_file);

				stat_run(pop,m_cur_run);// stat single run for algorithm analysis
				stat_rep_run();// arpso specific
				if ( is_final_run(m_cur_run,max_run) )
				{
					print_run_stat(stat_file,m_alg_stat,max_run);
					calc_run_avg_rep_stat();
					print_rep_run_stat(stat_file);// arpso specific
				}

				/*if ( !run_once )
					++(*pprog_dis);*/
			}// for every run

			print_avg_gen(stat_file,m_alg_stat.run_avg_gen);
			// stat and output average time per run by second
			m_alg_stat.run_avg_time=elapsed_t.elapsed();
			m_alg_stat.run_avg_time/=(max_run*1.0);
			print_avg_time(stat_file,m_alg_stat.run_avg_time);

			write_stat_vals();

			cout<<endl;// flush cout output
			return 0;
		}// end function Run

	}// namespace psobc
}// namespace pso
