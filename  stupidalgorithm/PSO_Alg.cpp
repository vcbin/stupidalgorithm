#include <boost/progress.hpp>
#include <boost/timer.hpp>

#include "pso_alg.h"

using std::vector;
using std::cout;
using std::string;
using std::ofstream;
using std::ostream;
using std::ios;
using std::endl;
using std::max_element;
using std::swap;

using boost::shared_ptr;
using boost::progress_display;
using boost::timer;

using boost::uniform_01;
using boost::mt19937;
using boost::variate_generator;

extern boost::mt19937 gen;

namespace pso
{
	void alloc_bound_handle(int ob_type,shared_ptr<bound_handle> &pBH)
	{
		switch (ob_type)
		{
		case obtype_absorb:
			pBH=shared_ptr<bound_handle>(new bh_absorb);
			break;
		case obtype_damp:
			pBH=shared_ptr<bound_handle>(new bh_damp);
			break;
		default:
			pBH=shared_ptr<bound_handle>(new bh_reflect);
		}
	}

		void pso_alg::initialize()
		{
			com_alg::initialize();

			// initialize particle velocity 2D vector
			size_t pop_size=m_ppara->get_pop_size();
			int num_dims=m_ppara->get_dim();
			m_alg_stat.delta_x.resize(pop_size);
			m_outbound.resize(pop_size);
			m_v_dim_diver.resize(num_dims);
			int max_run=m_ppara->get_max_run();
			m_gen_vel_div_val.resize(max_run);
			size_t i;
			for ( i=0;i<pop_size;i++ )
				m_outbound[i].assign(num_dims,false);
		}// end function initialize

		void pso_alg::reinit_x_range()
		{
			size_t num_dims=m_ppara->get_dim();
			const val_range &val_bnd=m_ppara->get_val_bnd();
			m_x_rng.resize(num_dims);

			size_t j;
			for ( j=0;j<num_dims;j++ )
			{
				m_x_rng[j].max_val=val_bnd[j].min_val;
				m_x_rng[j].min_val=val_bnd[j].max_val;
			}
		}// end function reinit_x_range

		void pso_alg::set_dynamic_vmax()
		{
			// set initial vector vmax for vmax of DYNAMIC type
			size_t num_dims=m_ppara->get_dim();
			double range_cof=m_ppara->get_vmax_cof();
			size_t j;
			for ( j=0;j<num_dims;j++ )
				m_vmax[j]=range_cof*(m_x_rng[j].max_val-m_x_rng[j].min_val);
		}

		void pso_alg::set_static_vmax()
		{
			// set initial vector vmax for vmax of STATIC type
			size_t num_dims=m_ppara->get_dim();
			val_range &val_bound=m_ppara->get_val_bnd();// NOT initialization range!!
			double range_cof=m_ppara->get_vmax_cof();
			size_t j;
			for ( j=0;j<num_dims;j++ )
				m_vmax[j]=range_cof*(val_bound[j].max_val-val_bound[j].min_val);
		}

		// invoked AFTER initial population has generated
		void pso_alg::set_ini_max_velocity()
		{
			size_t num_dims=m_ppara->get_dim();
			int vmax_type=m_ppara->get_vmax_type();
			m_vmax.resize(num_dims);
			if ( dynamic_rng==vmax_type )
				set_dynamic_vmax();
			else
				set_static_vmax();
		}// end function set_ini_max_velocity

		void pso_alg::update_speed(population &pop)
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
						m_alg_stat.delta_x[i][j]= m_omega*m_alg_stat.delta_x[i][j] + 
							r1*CR1*(m_alg_stat.pbest[i].x[j]-pop[i].x[j]) +
							r2*CR2*(m_alg_stat.gbest.x[j]-pop[i].x[j]);

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

		void pso_alg::update_pop(population &pop)
		{
			size_t pop_size=pop.size();
			size_t i,j;
			size_t num_dims=m_ppara->get_dim();
			val_range val_bounds=m_ppara->get_val_bnd();
			int vmax_type=m_ppara->get_vmax_type();

			for ( i=0;i<pop_size;i++ )
			{
				for ( j=0;j<num_dims;j++ )
				{
					double prev_x=pop[i].x[j]; // debugging purpose
					pop[i].x[j] += m_alg_stat.delta_x[i][j];
					/*if ( abs(pop[i].x[j])<1.0E-3 ) // debugging purpose
					cout<<" ";*/
					// boundaries check
					bool high_bnd,low_bnd;
					low_bnd=pop[i].x[j] < val_bounds[j].min_val;
					if ( low_bnd )
						pop[i].x[j]=val_bounds[j].min_val;
					high_bnd=pop[i].x[j] > val_bounds[j].max_val;
					if ( high_bnd )
						pop[i].x[j]=val_bounds[j].max_val;
					if ( low_bnd || high_bnd )
						m_outbound[i][j]=true;// record O.B state for next speed update procedure
				}// for every dimenson
			}// for every particle
		}// end function update_pop

		void pso_alg::update_x_rng(const population &pop)
		{
			size_t pop_size=m_ppara->get_pop_size();
			size_t num_dims=m_ppara->get_dim();

			reinit_x_range();
			size_t i,j;
			for ( i=0;i<pop_size;i++ )
			{
				for ( j=0;j<num_dims;j++ )
				{
					// update dynamic x range
					if ( pop[i].x[j] < m_x_rng[j].min_val )
						m_x_rng[j].min_val=pop[i].x[j];
					if ( pop[i].x[j] > m_x_rng[j].max_val )
						m_x_rng[j].max_val=pop[i].x[j];
				}// for every dimenson
			}// for every particle
		}// end function update_x_rng

		void pso_alg::update_omega(int cur_gen)
		{
			m_omega=calc_omega(cur_gen);
			if ( 3==m_ppara->get_stop_type() ) // omega value lower bound check for stop on stagnation(no maximum generation value is specified)
			{
				double omega_min=m_ppara->get_omega_min();
				if ( m_omega < omega_min ) 
					m_omega=omega_min;
			}
		}

		void pso_alg::print_run_title(ostream &os)
		{
			os<<"\n"
				<<"gen"
				<<"\t"
				<<"best_so_far"
				<<"\t"
				<<"avg"
				<<"\t"
				<<"pos_diversity"
				<<"\t"
				<<"vel_diversity"
				<<"\t"
				<<"std"
				<<"\t"
				<<"stag_gen"
				<<"\t"
				<<"eval_count"
				<<"\t"
				<<"omega"
				<<"\t"
				<<"\t"
				<<"avg_radius";
		}// end function print_run_title

		void pso_alg::print_gen_stat(ostream &os,int cur_gen,const alg_stat &alg_stat)
		{
			// // set float-point precision for output
			// os.precision(8);
			// os.setf(ios::scientific,ios::floatfield);
			os<<"\n"
				<<cur_gen
				<<" "
				<<alg_stat.gbest.obj[0];
			os.precision(8);
			// os.setf(ios::fixed,ios::floatfield);
			os<<"\t"
				<<alg_stat.avg;

			os.precision(4);
			os<<"\t"
				<<alg_stat.pos_diver
				<<"\t"
				<<m_vel_diver;

			// os.setf(ios::scientific,ios::floatfield);
			os<<"\t"
				<<alg_stat.std
				<<"\t"
				<<alg_stat.num_stag
				<<"\t"
				<<alg_stat.all_eval_num;

			// os.setf(ios::fixed,ios::floatfield);
			os.precision(4);
			os<<"\t"
				<<m_omega;

			os.precision(8);
			os<<"\t"
				<<"\t"
				<<m_alg_stat.avg_radius;
		}// end function print_gen_stat

		void pso_alg::update_diversity(const population &pop)
		{
			size_t pop_size=m_ppara->get_pop_size();
			int num_dims=m_ppara->get_dim();

			size_t i;
			int j;
			// update centroid individual
			m_alg_stat.cen_ind.assign(num_dims,0.0);
			for ( i=0;i<pop_size;i++ )
			{
				m_alg_stat.cen_ind += pop[i].x;
			}// for every particle
			m_alg_stat.cen_ind /=  pop_size;

			m_alg_stat.pos_diver=0.0;
			for ( i=0;i<pop_size;i++ )
			{
				m_alg_stat.pos_diver += vec_distance(pop[i].x,m_alg_stat.cen_ind);
			}// for every particle
			m_alg_stat.pos_diver /=  pop_size;
			m_alg_stat.pos_diver /= (m_diag_len*1.0);
			// calculate velocity diversity
			vector<double> vec_avg(num_dims,0.0);
			for ( j=0;j<num_dims;j++ )
			{
				for ( i=0;i<pop_size;i++ )
					vec_avg[j] += m_alg_stat.delta_x[i][j];// for every particle

				vec_avg[j] /=  pop_size;
			}// for every dimension
			m_v_dim_diver.assign(num_dims,0.0);
			for ( j=0;j<num_dims;j++ )
			{
				for ( i=0;i<pop_size;i++ )
					m_v_dim_diver[j] += abs(m_alg_stat.delta_x[i][j]-vec_avg[j]);// for every particle

				m_v_dim_diver[j] /=  pop_size;
			}// for every dimension
			m_vel_diver=0.0;
			for ( j=0;j<num_dims;j++ )
				m_vel_diver += m_v_dim_diver[j];
			m_vel_diver /= (num_dims*1.0);
		}// end function update_diversity

		void pso_alg::set_orig_pop(population &pop)
		{
			com_alg::set_orig_pop(pop);
			int vmax_type=m_ppara->get_vmax_type();
			if ( dynamic_rng==vmax_type )
				update_x_rng(pop);
		}

		double pso_alg::calc_omega(int cur_gen)
		{
			int max_gen=m_ppara->get_max_gen();
			double w_max=m_ppara->get_omega_max();
			double w_min=m_ppara->get_omega_min();
			return ( w_max - (m_cur_gen/(max_gen*1.0))*(w_max-w_min) );
		}

		void pso_alg::record_vel_div(int cur_run)
		{
			m_gen_vel_div_val[cur_run].push_back(m_vel_diver);
		}

		void pso_alg::record_gen_vals(alg_stat &alg_stat,int cur_run)
		{
			problem_base::record_gen_vals(alg_stat,cur_run);
			record_vel_div(cur_run);
		}

		void pso_alg::write_avg_vel_div_per_gen()
		{
			ofstream avg_vel_div_file(m_avg_vel_div_path);
			print_gen_avg_val(avg_vel_div_file,m_gen_avg_vel_div_val);
		}

		void pso_alg::write_stat_vals()
		{
			com_alg::write_stat_vals();
			write_avg_vel_div_per_gen();
		}

		void pso_alg::calc_gen_avg_vals()
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
			// calculate generational avg position diversity/velocity diversity/search radius
			m_alg_stat.gen_avg_div_val.resize(min_gen,0.0);
			m_gen_avg_vel_div_val.resize(min_gen,0.0);
			m_alg_stat.gen_avg_rad_val.resize(min_gen,0.0);
			int j;
			for ( j=0;j<min_gen;j++ ) // for every generation
				for ( i=0;i<max_run;i++ ) // for every run
				{
					m_alg_stat.gen_avg_div_val[j] += m_alg_stat.gen_div_val[i][j];// calculate gen avg diversity
					m_gen_avg_vel_div_val[j] += m_gen_vel_div_val[i][j];// calculate gen avg vel diversity
					m_alg_stat.gen_avg_rad_val[j] += m_alg_stat.gen_rad_val[i][j];// calculate gen avg radius
				}
				m_alg_stat.gen_avg_div_val /= max_run;
				m_gen_avg_vel_div_val /= max_run;
				m_alg_stat.gen_avg_rad_val /= max_run;
		}// end function calc_gen_avg_vals

		void pso_alg::reset_run_stat()
		{
			com_alg::reset_run_stat();
			double omega_max=m_ppara->get_omega_max();
			m_omega=omega_max;
		}

		int pso_alg::run()
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
			// shared_ptr<progress_display> pprog_dis;// algorithm progress indicator from boost
			// alloc_prog_indicator(pprog_dis);

			// algorithm statistics output file
			ofstream stat_file(m_com_out_path.stat_path);
			// allocate stop condition object dynamically
			bool run_once=1==max_run;
			alloc_stop_cond();
			int ob_type=m_ppara->get_ob_type();
			alloc_bound_handle(ob_type,m_pBH);

			// iteration start
			for ( m_cur_run=0;m_cur_run<max_run;m_cur_run++ )
			{
				reset_run_stat();

				set_orig_pop(pop);
				update_diversity(pop);

				set_ini_max_velocity();

				record_gen_vals(m_alg_stat,m_cur_run);

				print_run_times(stat_file,m_cur_run+1);
				print_run_title(stat_file);
				// output original population statistics
				print_gen_stat(stat_file,1,m_alg_stat);

				m_cur_gen=1;
				while ( false==(*m_pstop_cond) ) // for every iteration
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
			m_alg_stat.run_avg_time/=(max_run*1.0);
			print_avg_time(stat_file,m_alg_stat.run_avg_time);

			print_best_x(stat_file,m_alg_stat.bst_ind);
			write_stat_vals();

			cout<<endl;// flush cout output
			return 0;
		}// end function Run

}// namespace pso

