#include "mode_alg.h"

#include "mosade_alg.h"
#include <cassert>

using std::vector;
using std::cout;
using std::string;
using std::ofstream;
using std::ostream;
using std::swap;
using std::endl;
using std::ostream;
using std::ios;
using std::copy;
using std::back_inserter;
using std::advance;
using std::random_shuffle;

using boost::shared_ptr;
// using boost::progress_display;
using boost::timer;

using boost::uniform_01;
using boost::uniform_int;
using boost::normal_distribution;
using boost::mt19937;
using boost::variate_generator;

using boost::mutex;
using boost::thread;

extern boost::mt19937 gen;

using namespace mop;

namespace de
{
	namespace mode
	{
		void mode_alg::record_gen_vals(alg_stat &alg_stat,int cur_run)
		{
			// record stat values of run
			record_diver(alg_stat,cur_run);
			record_radius(alg_stat,cur_run);
		}

		void mode_alg::write_stat_vals()
		{
			write_avg_diver_per_gen();
			write_avg_rad_per_gen();
		}

		// SPECIAL boundaries check
		void mode_alg::bound_check(double &tri_x,double orig_x,int dim)
		{
			const val_range& val_bounds=m_ppara->get_val_bnd();
			bool out_up_bnd,out_low_bnd;
			double low_bnd,high_bnd;
			low_bnd=val_bounds[dim].min_val;
			high_bnd=val_bounds[dim].max_val;
			out_low_bnd=tri_x < low_bnd;
			out_up_bnd=tri_x > high_bnd;
			if ( out_up_bnd || out_low_bnd )
			{
				m_alg_stat.all_ob_num++;
				/*normal_distribution<> norm_dist;
				variate_generator<mt19937&, normal_distribution<> > rnd_norm(gen, norm_dist);
				x=m_alg_stat.cen_ind[dim]+rnd_norm();*/

				//// uniformly randomize x within [lower_bound,upper_bound]
				//uniform_01<> dist;
				//variate_generator<mt19937&, uniform_01<> > rnd_num(gen, dist);
				//x=low_bnd+rnd_num()*(high_bnd-low_bnd);

				if ( out_low_bnd )
					tri_x=(orig_x+low_bnd)/2;
				else
					tri_x=low_bnd+(tri_x-high_bnd)/2;
				bound_check(tri_x,orig_x,dim);// RECURSIVE invocation
			}
		}// end function bound_check

		inline bool mode_alg::is_learn_gen(int cur_gen,int learn_p)
		{
			return ( 0==(cur_gen+1)%learn_p );
		}

		int mode_alg::run()
		{
			if ( !m_ppara )
				return -1;

			timer elapsed_t;
			// retrieve algorithm parameters
			int pop_size=m_ppara->get_pop_size();
			int max_pop_size=2*pop_size;
			int num_dims=m_ppara->get_dim();
			int num_obj=m_ppara->get_obj_num();
			int min_gen=m_ppara->get_max_gen();
			double pr_val,f_val;
			pr_val=m_ppara->get_pr();
			f_val=m_ppara->get_f();
			// my_sDE SPECIFIC
			// jDE SPECIFIC
			double f_low_bnd,f_up_bnd;
			f_low_bnd=m_ppara->get_f_low_bnd();
			f_up_bnd=m_ppara->get_f_up_bnd();
			double tau_1,tau_2;
			tau_1=m_ppara->get_tau_1();
			tau_2=m_ppara->get_tau_2();
			int out_interval=m_ppara->get_out_interval();
			int trunc_type=m_ppara->get_trunc_type();
			bool plot=m_ppara->get_plot_flag();
			string plot_script=m_ppara->get_plot_script();

			int m_cur_run;
			int max_run=m_ppara->get_max_run();// run/trial number
			// shared_ptr<progress_display> pprog_dis;// algorithm progress indicator from boost
			// alloc_prog_indicator(pprog_dis);

			// allocate original pop and trial pop
			population pop(pop_size);
			allocate_pop(pop,num_dims,stra_num,num_obj);
			population trial_pop;

			// generate algorithm statistics output file name
			ofstream stat_file(m_com_out_path.stat_path.c_str());
			// allocate stop condition object dynamically
			alloc_stop_cond();

			idx_array pop_idx(pop_size-1);

			// random U(0,1) generator
			uniform_01<> dist_01;
			variate_generator<mt19937&, uniform_01<> > rnd_01(gen, dist_01);

			// generator for random DIMENSION index
			uniform_int<> dist_dim(0,num_dims-1);
			variate_generator<mt19937&, uniform_int<> > rnd_dim_idx(gen, dist_dim);

			individual trial_ind;
			allocate_ind(trial_ind,num_dims,stra_num,num_obj);
			// iteration start
			for ( m_cur_run=0;m_cur_run<max_run;m_cur_run++ )
			{
				bool has_stag;
				int stag_gen;
				reset_run_stat();
				m_de_stat.reset();
				/*m_succ_f.clear();
				m_succ_cr.clear();*/
				int z;
				// f,pr initialial value
				for ( z=0;z<pop_size;z++ )
				{
					pop[z].stra[f].assign(1,f_val);
					pop[z].stra[pr].assign(1,pr_val);
				}
				set_orig_pop(pop);
				update_diversity(pop);

				calc_de_para_stat(pop);
				record_de_para_stat(m_cur_run);

				print_run_times(stat_file,m_cur_run+1);
				print_run_title(stat_file);
				// output original population statistics

				m_cur_gen=1;
				stag_gen=0;
				shared_ptr<ofstream> ppop_file;
				shared_ptr<mutex> ppop_mut;
				
				while ( false==(*m_pstop_cond) ) // for every iteration
				{
					m_de_stat.reset();
					has_stag=true;
					int rnd_dim;
					double dim_mut_chance;
					int i,j,k;
					/*double f_i;
					double pr_i;*/
					if ( is_output_gen(m_cur_gen,out_interval) )
					{
						ppop_file=shared_ptr<ofstream>(new ofstream("Output//all_pop.out"));
						ppop_mut=shared_ptr<mutex>((new mutex));
					}
					trial_pop.clear();
					trial_pop.reserve(max_pop_size);// operator =
					////  recalculate succ_f_mean and succ_f_sigma periodically
     //               if ( is_learn_gen(m_cur_gen,learn_p) )
					//{
					//	int f_size=m_succ_f.size();
					//	int pr_size=m_succ_cr.size();
					//	if ( f_size && pr_size ) 
					//		calc_bi_norm_var(m_succ_f,m_succ_cr,m_bi_norm_var);
					//	else
					//	{

					//	}
					//	m_succ_f.clear();
					//	m_succ_cr.clear();
					//}// if ( is_learn_gen(m_cur_gen,learn_p) )

					for ( i=0;i<pop_size;i++ )
					{
						// generating three mutually different individual index using random shuffle
						// initialize index vector

						for ( k=0;k<pop_size-1;k++ )
						{
							if ( k<i )
								pop_idx[k]=k;
							else
								pop_idx[k]=(k+1)%pop_size;// EXCLUDE i
						}
						random_shuffle(pop_idx.begin(),pop_idx.end());
						int i1,i2,i3;
						// int i4,i5;// i!=i1!=i2!=i3!=i4!=i5
						i1=pop_idx[0];
						i2=pop_idx[1];
						i3=pop_idx[2];
						/*i4=arc_idx[3];
						i5=arc_idx[4];*/

						//double pr_chance=rnd_01();
						//if ( pr_chance<=tau_2 )
						//{
						//	pr_i=rnd_01();
						//	trial_ind.stra[pr][0]=pr_i;
						//}
						//else
						//	pr_i=trial_ind.stra[pr][0];

						//// scaling factor F self-adaptive update equation
						//double f_chance=rnd_01();
						//if ( f_chance<=tau_1 )
						//{
						//	f_i=f_low_bnd+rnd_01()*f_up_bnd;
						//	trial_ind.stra[f][0]=f_i;
						//}
						//else
						//	f_i=trial_ind.stra[f][0];// keep unchanged at thsi iteration

						rnd_dim=rnd_dim_idx();// choose a random dimension as the mutation target

						for ( j=0;j<num_dims;j++ )
						{
							dim_mut_chance=rnd_01();
							
							if ( rnd_dim==j || dim_mut_chance<=pr_val )
							{
								// insufficent elitist size,generate perturbation from current population rather than external elitist archive
								trial_ind.x[j]=pop[i1].x[j]+f_val*(pop[i2].x[j]-pop[i3].x[j]);
								//+f_i*(trial_pop[i4].x[j]-trial_pop[i5].x[j]);
								// boundaries check
								bound_check(trial_ind.x[j],pop[i].x[j],j);
							}
							else
								trial_ind.x[j]=pop[i].x[j];
						}// for every dimension

						eval_ind(trial_ind,*m_pfunc,m_alg_stat);
						int comp_res=check_dominance(trial_ind,pop[i]);
						if ( worse!=comp_res )
						{
							if ( better==comp_res )
								trial_pop.push_back(trial_ind);
							else
							{
								trial_pop.push_back(trial_ind);
								trial_pop.push_back(pop[i]);
							}
						}
						else
							trial_pop.push_back(pop[i]);
					}// for every point
					// evaluate pop
					
					fill_nondominated_sort(trial_pop,pop,pop_size,trunc_type);

					if ( is_output_gen(m_cur_gen,out_interval) )
					{
						update_search_radius();
						update_diversity(pop);

						calc_de_para_stat(pop);
						record_gen_vals(m_alg_stat,m_cur_run);
						record_de_para_stat(m_cur_run);

						/*if ( run_once )
						++(*pprog_dis);*/

						// plot current population and external archive
						output_collection(*ppop_file,pop.begin(),pop.end());
						*ppop_file<<"\n"<<"\n";// output seperator
						output_if(*ppop_file,pop.begin(),pop.end(),front_pred());
						ppop_file->flush();
						if ( plot )
							thread(fwd_plot_fun,ppop_mut,m_cur_gen,
							0,
							0,is_final_out_gen(m_cur_gen,out_interval,min_gen),
							plot_script);
						// system("gnuplot plot_all_point_2d.p");
					}

					m_cur_gen++;
				}// while single run termination criterion is not met
				perf_indice p_ind,nsga2_p_ind;
				d_mat best_pop;
				copy_obj_if(pop.begin(),pop.end(),best_pop,front_pred());

				d_mat nsga2_best_pop;
				load_pop("nsga2_zdt3_best_pop.out",2,nsga2_best_pop);
				zdt3_assess(nsga2_best_pop,1000,point(11,11),nsga2_p_ind);
				cout<<"\n"
					<<"nsga2 results:"
					<<"\n"
					<<"convergence metric gamma="<<nsga2_p_ind.gamma
					<<"\n"
					<<"frontier diversity metric delta="<<nsga2_p_ind.delta
					<<"\n"
					<<"dominance metric hyper-volume="<<nsga2_p_ind.hv
					<<"\n";
				zdt3_assess(best_pop,1000,point(11,11),p_ind);
				cout<<"\n"
					<<"outbound count="<<m_alg_stat.all_ob_num
					<<"\n"
					<<"population diversity="<<m_alg_stat.pos_diver
					<<"\n"
					<<"mean search radius="<<m_alg_stat.avg_radius
					<<"\n"
					<<"stagnation indicator stag_gen="<<stag_gen
					<<"\n"
					<<"convergence metric gamma="<<p_ind.gamma
					<<"\n"
					<<"frontier diversity metric delta="<<p_ind.delta
					<<"\n"
					<<"dominance metric hyper-volume="<<p_ind.hv
					<<"\n";
				// single run end
				
				/*if ( !run_once )
				++(*pprog_dis);*/
			}// for every run

			print_avg_gen(stat_file,m_alg_stat.run_avg_gen);
			// stat and output average time per run by second
			m_alg_stat.run_avg_time=elapsed_t.elapsed();
			m_alg_stat.run_avg_time /= (max_run*1.0);
			print_avg_time(stat_file,m_alg_stat.run_avg_time);
			write_stat_vals();
			cout<<endl;// flush cout output
			return 0;
		}// end function Run

	}// end namespace mode
}// end namespace de
