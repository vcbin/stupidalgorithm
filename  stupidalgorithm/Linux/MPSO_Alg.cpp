#include "MPSO_Alg.h"

using std::cout;
using std::endl;
using std::ofstream;
using std::ostream;
using std::ios;

using boost::timer;
// using boost::progress_display;
using boost::shared_ptr;

using boost::uniform_01;
using boost::mt19937;
using boost::variate_generator;

extern boost::mt19937 gen;

namespace pso
{
	namespace mpso
	{
		void mpso_alg::initialize()
		{
			pso_alg::initialize();
			int max_run=m_ppara->get_max_run();
			m_mpso_stat.alloc(max_run);
			// intialize acceleration array
			size_t pop_size=m_ppara->get_pop_size();
			accle.assign(pop_size,1.0);
			// set attraction/repulsion parameters
			m_dl_cof=m_ppara->get_dh();
			m_dh_cof=m_ppara->get_dh();
			set_diver_bound();
		}

		void mpso_alg::record_attr_perc(Mpso_stat &mpso_stat,int cur_run)
		{
			mpso_stat.gen_attr_val[cur_run].push_back((1.0-mpso_stat.rep_perc)*100);// percentage number
		}

		void mpso_alg::record_gen_vals(alg_stat &alg_stat,int cur_run)
		{
			pso_alg::record_gen_vals(alg_stat,cur_run);
			record_attr_perc(m_mpso_stat,cur_run);
		}

		void mpso_alg::update_rep_stat(int cur_gen)
		{
			size_t pop_size=m_ppara->get_pop_size();
			size_t rep_size=m_mpso_stat.rep_size;
			m_mpso_stat.rep_perc=rep_size/(1.0*pop_size);
			m_mpso_stat.avg_rep_size += rep_size;
			m_mpso_stat.avg_rep_perc += m_mpso_stat.rep_perc;
			if ( rep_size>0 )
			{
				if ( m_mpso_stat.first_rep_flag )
				{
					m_mpso_stat.first_rep_gen=cur_gen;
					m_mpso_stat.first_rep_flag=false;
				}
				m_mpso_stat.rep_gen_num++;
			}
		}

		void mpso_alg::calc_run_rep_stat()
		{
			m_mpso_stat.rep_gen_perc=m_mpso_stat.rep_gen_num/(1.0*m_cur_gen);
			int rep_gen_num=m_mpso_stat.rep_gen_num;
			if ( rep_gen_num>0 )
			{
				m_mpso_stat.avg_rep_size /= rep_gen_num;
				m_mpso_stat.avg_rep_perc /= rep_gen_num;
			}
		}

		void mpso_alg::stat_rep_run()
		{
			if ( m_mpso_stat.rep_gen_num>0 )
				m_mpso_stat.rep_run++;

			m_mpso_stat.avg_rep_gen_num += m_mpso_stat.rep_gen_num;
			m_mpso_stat.avg_rep_gen_perc += m_mpso_stat.rep_gen_perc;
			m_mpso_stat.avg_avg_rep_size += m_mpso_stat.avg_rep_size;
			m_mpso_stat.avg_avg_rep_size_perc += m_mpso_stat.avg_rep_perc;
		}

		void mpso_alg::calc_run_avg_rep_stat()
		{
			int rep_run=m_mpso_stat.rep_run;
			if ( rep_run>0 )
			{
				m_mpso_stat.avg_rep_gen_num /= rep_run;
				m_mpso_stat.avg_rep_gen_perc /= rep_run;
				m_mpso_stat.avg_avg_rep_size /= rep_run;
				m_mpso_stat.avg_avg_rep_size_perc /= rep_run;
			}
		}

		void mpso_alg::calc_gen_avg_vals()
		{
			// find min gen
			int max_run=m_ppara->get_max_run();
			int min_gen;
			int i;
			int cur_gen_num;
			min_gen=m_alg_stat.gen_div_val[0].size();
			for ( i=1;i<max_run;i++ )
			{
				cur_gen_num=m_alg_stat.gen_div_val[i].size();
				if ( cur_gen_num < min_gen )
					min_gen=cur_gen_num;
			}

			m_alg_stat.gen_avg_div_val.resize(min_gen,0.0);
			m_alg_stat.gen_avg_rad_val.resize(min_gen,0.0);
			m_gen_avg_vel_div_val.resize(min_gen,0.0);
			m_mpso_stat.gen_avg_attr_val.resize(min_gen,0.0);
			int j;
			for ( j=0;j<min_gen;j++ ) // for every generation
				for ( i=0;i<max_run;i++ ) // for every run
				{
					m_alg_stat.gen_avg_div_val[j] += m_alg_stat.gen_div_val[i][j];// calculate gen avg diversity
					m_alg_stat.gen_avg_rad_val[j] += m_alg_stat.gen_rad_val[i][j];// calculate gen avg radius
					m_gen_avg_vel_div_val[j] += m_gen_vel_div_val[i][j];// calculate gen avg vel diversity
					m_mpso_stat.gen_avg_attr_val[j] += m_mpso_stat.gen_attr_val[i][j];// calculate gen avg attraction percentage
				}
				m_alg_stat.gen_avg_div_val /= max_run;
				m_alg_stat.gen_avg_rad_val /= max_run;
				m_gen_avg_vel_div_val /= max_run;
				m_mpso_stat.gen_avg_attr_val /= max_run;

		}// end function calc_gen_avg_div

		// mpso algo specific,record average attraction phrase percentage of all runs in every generation
		void mpso_alg::write_avg_attr_per_gen()
		{
			ofstream avg_attr_file(m_avg_attr_path.c_str());
			print_gen_avg_val(avg_attr_file,m_mpso_stat.gen_avg_attr_val);
		}

		// record average values of all runs in every generation
		void mpso_alg::write_stat_vals()
		{
			pso_alg::write_stat_vals();
			write_avg_attr_per_gen();
		}

		// update acceleration for every particle
		// TIME-CONSUMING part of mpso algorithm compared to other pso algos!
		void mpso_alg::update_accle(population &pop,int cur_gen)
		{
			size_t pop_size=pop.size();
			size_t num_dims=m_ppara->get_dim();
			size_t i,j;

			for ( i=0;i<pop_size;i++ )
			{
				for ( j=0;j<num_dims;j++ )
				{
					double dis=abs(vec_distance(pop[i].x,m_alg_stat.cen_ind));
					if ( dis < m_d_l && accle[i]==1.0 )
					{
						accle[i]=-1.0;
						update_rep_size();
					}
					if ( dis > m_d_h && accle[i]==-1.0 )
					{
						accle[i]=1.0;
					}
				}// for every dimension
			}// for every particle
		}// end function update_accle

		void mpso_alg::update_speed(population &pop)
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
							accle[i]*(
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
		}// end function update_speed

		void mpso_alg::print_run_title(ostream &os)
		{
			output_base::print_run_title(os);
			os<<"\t"
				<<"rep_particle_num"
				<<"\t"
				<<"repulsion_percent";
		}

		void mpso_alg::print_rep_stat(std::ostream &os)
		{
			os<<"\n"
				<<"first_repulsion_generation:";
			bool had_rep=m_mpso_stat.first_rep_gen != -1;
			if ( had_rep )
				os<<m_mpso_stat.first_rep_gen+1;
			else
				os<<"NONE";

			os<<"\n"
				<<"repulsion_generations:"
				<<m_mpso_stat.rep_gen_num;
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
				<<m_mpso_stat.rep_gen_perc*100.0
				<<"%";

			os<<"\n"
				<<"avg_rep_size:"
				<<m_mpso_stat.avg_rep_size;
			if ( had_rep )
				os<<"/"
				<<m_ppara->get_pop_size();
			os<<"\t"
				<<"avg_rep_perc:"
				<<m_mpso_stat.avg_rep_perc*100.0
				<<"%";
		}// end function print_rep_stat

		void mpso_alg::print_gen_stat(ostream &os,int cur_gen,const alg_stat &alg_stat)
		{
			pso_alg::print_gen_stat(os,cur_gen,alg_stat);
			os<<"\t"
				<<"\t";
			// percent number format
			os<<m_mpso_stat.rep_size;

			os.setf(ios::fixed,ios::floatfield);
			os.precision(2);
			os<<"\t"
				<<m_mpso_stat.rep_perc*100.0
				<<"%";
		}

		void mpso_alg::print_rep_run_stat(ostream &os)
		{
			os<<"\n"
				<<"repulsion_run:"
				<<m_mpso_stat.rep_run
				<<"/"
				<<m_ppara->get_max_run()
				<<"\t"
				<<"avg_rep_gen:"
				<<m_mpso_stat.avg_rep_gen_num;
			int stoptype=m_ppara->get_stop_type();
			if ( stoptype_stag!=stoptype && stoptype_delta!=stoptype )
			{
				os<<"/"
					<<m_ppara->get_max_gen();
			}
			os.precision(2);
			os.setf(ios::fixed,ios::floatfield);
			os<<"\t"
				<<"avg_rep_gen_perc:"
				<<m_mpso_stat.avg_rep_gen_perc*100.0
				<<"%";

			os<<"\n"
				<<"avg_avg_rep_size:"
				<<m_mpso_stat.avg_avg_rep_size
				<<"/"
				<<m_ppara->get_pop_size()
				<<"\t"
				<<"avg_avg_rep_size_perc:"
				<<m_mpso_stat.avg_avg_rep_size_perc*100.0
				<<"%";
		}// end function print_rep_run_stat

		int mpso_alg::run()
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
				m_mpso_stat.reset_run_stat();// mpso specific

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
					m_mpso_stat.reset_gen_stat();// MPSO specific procedure

					update_omega(m_cur_gen);
					update_accle(pop,m_cur_gen);// MPSO specific procedure
					update_rep_stat(m_cur_gen);// stat repulsion particle percentage
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

				// mpso specific output
				calc_run_rep_stat();
				print_rep_stat(stat_file);

				stat_run(pop,m_cur_run);// stat single run for algorithm analysis
				stat_rep_run();// mpso specific
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

	}// namespace MPSO

}// namespace pso
