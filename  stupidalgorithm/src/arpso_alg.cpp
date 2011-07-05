#include "arpso_alg.h"

using boost::timer;
//// using boost::progress_display;
using boost::shared_ptr;

using std::cout;
using std::endl;
using std::ofstream;
using std::ostream;
using std::ios;

using boost::uniform_01;
using boost::mt19937;
using boost::variate_generator;

extern boost::mt19937 gen;

namespace pso
{
	namespace arpso
	{
		void arpso_alg::initialize()
		{
			pso_alg::initialize();
			m_d_l=m_ppara->get_dh();
			m_d_h=m_ppara->get_dh();
		}

		void arpso_alg::update_speed(population &pop)
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
			size_t i,j;
			for ( i=0;i<pop_size;i++ )
			{
				for ( j=0;j<num_dims;j++ )
				{
					if ( !m_outbound[i][j] ) // if particle didn't outbound
					{
						r1=rnd_num();
						r2=rnd_num();
						m_alg_stat.delta_x[i][j]= m_omega*m_alg_stat.delta_x[i][j] + 
							m_dir_accl*(
							r1*CR1*(m_alg_stat.pbest[i].x[j]-pop[i].x[j]) +
							r2*CR2*(m_alg_stat.gbest.x[j]-pop[i].x[j])
							);

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

				}// for every dimension
			}// for every particle
		}// update_speed

		void arpso_alg::update_dir(population &pop)
		{
			size_t pop_size=pop.size();
			size_t num_dims=m_ppara->get_dim();
			size_t i,j;

			for ( i=0;i<pop_size;i++ )
			{
				for ( j=0;j<num_dims;j++ )
				{
					if ( m_alg_stat.pos_diver < m_d_l && m_dir_accl==1.0 )
					{
						m_dir_accl=-1.0;
					}
					if ( m_alg_stat.pos_diver > m_d_h && m_dir_accl==-1.0 )
					{
						m_dir_accl=1.0;
					}
				}// for every dimension
			}// for every particle
		}// end function update_dir

		void arpso_alg::stat_rep_gen(int cur_gen)
		{
			if ( m_dir_accl < 0 )
			{
				m_arpso_stat.rep_gen_num++;
				if ( m_arpso_stat.first_rep_flag )
				{
					m_arpso_stat.first_rep_gen=cur_gen;
					m_arpso_stat.first_rep_flag=false;
				}
			}
		}

		void arpso_alg::print_rep_stat(ostream &os)
		{
			os<<"\n"
				<<"first_repulsion_generation:";
			bool had_rep=m_arpso_stat.first_rep_gen != -1;
			if ( had_rep )
				os<<m_arpso_stat.first_rep_gen+1;
			else
				os<<"NONE";

			os<<"\n"
				<<"repulsion_generations:"
				<<m_arpso_stat.rep_gen_num;
			if ( had_rep )
			{
				int stoptype=m_ppara->get_stop_type();
				os<<"/";
				if ( stoptype_stag!=stoptype && stoptype_delta!=stoptype )
					os<<m_ppara->get_max_gen();
				else
					os<<m_cur_gen;
			}
			// percent number format
			os.setf(ios::fixed,ios::floatfield);
			os.precision(2);
			os<<"\t"
				<<"percentage:"
				<<m_arpso_stat.rep_gen_perc*100.0
				<<"%";
		}

		void arpso_alg::stat_rep_run()
		{
			if ( m_arpso_stat.rep_gen_num>0 )
				m_arpso_stat.rep_run++;
			m_arpso_stat.avg_rep_gen_num += m_arpso_stat.rep_gen_num;
			m_arpso_stat.avg_rep_gen_perc += m_arpso_stat.rep_gen_perc;
		}

		void arpso_alg::print_rep_run_stat(ostream &os)
		{
			os<<"\n"
				<<"repulsion_run:"
				<<m_arpso_stat.rep_run
				<<"/"
				<<m_ppara->get_max_run()
				<<"\t"
				<<"avg_rep_gen:"
				<<m_arpso_stat.avg_rep_gen_num;
			int stoptype=m_ppara->get_stop_type();
			if ( stoptype_stag!=stoptype && stoptype_delta!=stoptype )
			{
				os<<"/"
					<<m_ppara->get_max_gen();
			}
			os<<"\t"
				<<"avg_rep_gen_perc:"
				<<m_arpso_stat.avg_rep_gen_perc*100.0
				<<"%";
		}// end function print_rep_run_stat

		void arpso_alg::calc_rep_gen_perc()
		{
			m_arpso_stat.rep_gen_perc=m_arpso_stat.rep_gen_num/(1.0*m_cur_gen);
		}

		void arpso_alg::calc_run_avg_rep_stat()
		{
			int rep_run=m_arpso_stat.rep_run;
			if ( rep_run>0 )
			{
				m_arpso_stat.avg_rep_gen_num /= rep_run;
				m_arpso_stat.avg_rep_gen_perc /= rep_run;
			}
		}

		int arpso_alg::run()
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
			// shared_ptr<progress_display> pprog_dis;// algorithm progress indicator from boost
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
				m_arpso_stat.reset();// arpso specific

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
					update_dir(pop);// ARPSO specific procedure
					stat_rep_gen(m_cur_gen);// ARPSO specific procedure for repulsion stat
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
					print_rep_run_stat(stat_file);
				}
				/*if ( !run_once )
					++(*pprog_dis);*/
			}// for every run

			print_avg_gen(stat_file,m_alg_stat.run_avg_gen);
			// stat and output average time per run by second
			m_alg_stat.run_avg_time=elapsed_t.elapsed();
			m_alg_stat.run_avg_time/=(max_run*1.0);
			print_avg_time(stat_file,m_alg_stat.run_avg_time);
			
			print_best_x(stat_file,m_alg_stat.bst_ind);
			write_stat_vals();

			cout<<endl;// flush cout output
			return 0;
		}// end function Run

	}// namespace arpso
}// namespace pso
