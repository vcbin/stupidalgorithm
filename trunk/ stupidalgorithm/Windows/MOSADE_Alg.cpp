#include "MOSADE_Alg.h"
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

using boost::shared_ptr;
using boost::progress_display;
using boost::timer;

using boost::uniform_01;
using boost::uniform_int;
using boost::normal_distribution;
using boost::mt19937;
using boost::variate_generator;

using boost::mutex;
using boost::thread;
using boost::lock_guard;

extern boost::mt19937 gen;

using namespace mop;

namespace de
{
	namespace mosade
	{
		void mosade_alg::crowd_tour(const elite_archive &ext_archive,individual &trial,individual &ind,int index,int trunc_type,int neighbor_num)
		{
			// choose the less crowded one with respect to elitist archive
			switch (trunc_type)
			{
			case crowd_harmonic:
				calc_crowd_harmonic(ext_archive,trial,neighbor_num);
				calc_crowd_harmonic(ext_archive,ind,neighbor_num);
				break;
			case crowd_entropy:
				calc_crowd_entropy(ext_archive,trial);
				calc_crowd_entropy(ext_archive,ind);
				break;
			default:
				calc_crowd_dist(ext_archive,trial);
				calc_crowd_dist(ext_archive,ind);
			}
			
			if ( trial.crowd_dist!=ind.crowd_dist )
			{
				if ( trial.crowd_dist<ind.crowd_dist )// choose CLOSER point to the current frontier
				{
					m_alg_stat.delta_x[index]=trial.x-ind.x;// update step size delta_x
					ind=trial;
				}
				else
					m_alg_stat.delta_x[index].assign(m_ppara->get_dim(),0.0);
			}
			else
			{
				// crowding indicator the same,choose one randomly
				uniform_01<> dist_uni_01;
				variate_generator<mt19937&, uniform_01<> > rnd_uni_01(gen, dist_uni_01);
				if ( rnd_uni_01()<=0.5 )
				{
					m_alg_stat.delta_x[index]=trial.x-ind.x;// update step size delta_x
					ind=trial;
				}
				else
					m_alg_stat.delta_x[index].assign(m_ppara->get_dim(),0.0);
			}
		}// end function crowd tour

		void mosade_alg::record_gen_vals(alg_stat &alg_stat,int cur_run)
		{
			// record stat values of run
			record_diver(alg_stat,cur_run);
			record_radius(alg_stat,cur_run);
		}

		void mosade_alg::write_stat_vals()
		{
			write_avg_diver_per_gen();
			write_avg_rad_per_gen();
		}

		// SPECIAL boundaries check
		void mosade_alg::bound_check(double &tri_x,double orig_x,int dim)
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

		int mosade_alg::run()
		{
			if ( !m_ppara )
				return -1;

			timer elapsed_t;
			// retrieve algorithm parameters
			int pop_size=m_ppara->get_pop_size();
			int num_dims=m_ppara->get_dim();
			int num_obj=m_ppara->get_obj_num();
			int min_gen=m_ppara->get_max_gen();
			double vtr=m_ppara->get_vtr();
			double pr_val,f_val;
			pr_val=m_ppara->get_pr();
			f_val=m_ppara->get_f();
			int max_archive_size=m_ppara->get_max_archive();
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
			int trunc_count;
			// shared_ptr<progress_display> pprog_dis;// algorithm progress indicator from boost
			// alloc_prog_indicator(pprog_dis);

			// allocate original pop and trial pop
			population pop(pop_size);
			allocate_pop(pop,num_dims,stra_num,num_obj);
			population trial_pop(pop_size);
			allocate_pop(trial_pop,num_dims,stra_num,num_obj);

			// generate algorithm statistics output file name
			ofstream stat_file(m_com_out_path.stat_path);
			// allocate stop condition object dynamically
			bool run_once=(1==max_run);
			alloc_stop_cond();

			idx_array arc_idx;
			idx_array pop_idx(pop_size-1);
			int shuffle_size;

			// random U(0,1) generator
			uniform_01<> dist_01;
			variate_generator<mt19937&, uniform_01<> > rnd_01(gen, dist_01);

			// generator for random DIMENSION index
			uniform_int<> dist_dim(0,num_dims-1);
			variate_generator<mt19937&, uniform_int<> > rnd_dim_idx(gen, dist_dim);

			// iteration start
			for ( m_cur_run=0;m_cur_run<max_run;m_cur_run++ )
			{
				bool has_stag;
				int stag_gen;
				reset_run_stat();
				m_de_stat.reset();
				trunc_count=0;
				int z;
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

				// initialize elitist archive
				elite_archive ext_archive;
				int i;
				for ( i=0;i<pop_size;i++ )
					update_archive(ext_archive,pop[i]);
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
					double f_i;
					double pr_i;
					if ( is_output_gen(m_cur_gen,out_interval) )
					{
						ppop_file=shared_ptr<ofstream>(new ofstream("Output\\all_pop.out"));
						ppop_mut=shared_ptr<mutex>((new mutex));
					}

					trial_pop=pop;// operator =
					for ( i=0;i<pop_size;i++ )
					{
						// generating three mutually different individual index using random shuffle
						// initialize index vector
						shuffle_size=ext_archive.size()-1;
						arc_idx.resize(shuffle_size);
						init_idx_array(arc_idx,shuffle_size);
						
						// random shuffle
						std::random_shuffle(arc_idx.begin(),arc_idx.end());
						for ( k=0;k<pop_size-1;k++ )
						{
							if ( k<i )
								pop_idx[k]=k;
							else
								pop_idx[k]=(k+1)%pop_size;// EXCLUDE i
						}
						std::random_shuffle(pop_idx.begin(),pop_idx.end());
						int i1,i2,i3;
						// int i4,i5;// i!=i1!=i2!=i3!=i4!=i5
						i1=arc_idx[0];
						i2=pop_idx[0];
						i3=pop_idx[1];
						/*i2=arc_idx[1];
						i3=arc_idx[2];*/
						/*i4=arc_idx[3];
						i5=arc_idx[4];*/

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

						rnd_dim=rnd_dim_idx();// choose a random dimension as the mutation target

						for ( j=0;j<num_dims;j++ )
						{
							dim_mut_chance=rnd_01();
							
							if ( rnd_dim==j || dim_mut_chance<=pr_i )
							{
								if ( !ext_archive.empty() )
								{
									// generator for random elitist archive member index
									/*uniform_int<> dist_uni_arc(0,ext_archive.size()-1);
									variate_generator<mt19937&, uniform_int<> > rnd_arc_idx(gen, dist_uni_arc);*/
									elite_archive::iterator itr_arc_i1;
									elite_archive::iterator itr_arc_i2,itr_arc_i3,itr_arc_i4,itr_arc_i5;
									itr_arc_i4=itr_arc_i5=itr_arc_i2=itr_arc_i3=itr_arc_i1=ext_archive.begin();
									advance(itr_arc_i1,i1);// forward the list pointer to the destination
									/*advance(itr_arc_i2,i2);
									advance(itr_arc_i3,i3);*/
									/*advance(itr_arc_i4,i4);
									advance(itr_arc_i5,i5);*/
									trial_pop[i].x[j]=itr_arc_i1->x[j]+f_i*(trial_pop[i2].x[j]-trial_pop[i3].x[j]);
									//+f_i*(itr_arc_i4->x[j]-itr_arc_i5->x[j]);
								}
								else
									// insufficent elitist size,generate perturbation from current population rather than external elitist archive
									trial_pop[i].x[j]=trial_pop[i1].x[j]+f_i*(trial_pop[i2].x[j]-trial_pop[i3].x[j]);
								//+f_i*(trial_pop[i4].x[j]-trial_pop[i5].x[j]);
								// boundaries check
								bound_check(trial_pop[i].x[j],pop[i].x[j],j);
							}
						}// for every dimension

						eval_ind(trial_pop[i],*m_pfunc,m_alg_stat);
						int comp_res=check_dominance(trial_pop[i],pop[i]);
						if ( worse!=comp_res )
						{
							if ( true==update_archive(ext_archive,trial_pop[i]) )
								has_stag=false;

							if ( indiff==comp_res )
								crowd_tour(ext_archive,trial_pop[i],pop[i],i,trunc_type,2);
						}
					}// for every point
					// evaluate pop
					trunc_external(ext_archive,trunc_type,max_archive_size,trunc_count,2);
					if ( has_stag )
						stag_gen++;// update stagnation generations
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
						output_collection(*ppop_file,ext_archive.begin(),ext_archive.end());
						ppop_file->flush();
						if ( plot )
							thread(fwd_plot_fun,ppop_mut,m_cur_gen,
							ext_archive.size(),
							stag_gen,is_final_out_gen(m_cur_gen,out_interval,min_gen),
							plot_script);
						// system("gnuplot plot_all_point_2d.p");
					}

					m_cur_gen++;
				}// while single run termination criterion is not met
				perf_indice p_ind;
				d_mat best_pop;
				copy_obj(ext_archive.begin(),ext_archive.end(),best_pop);
				zdt3_assess(best_pop,1000,point(11,11),p_ind);
				cout<<"\n"
					<<"outbound count="<<m_alg_stat.all_ob_num
					<<"\n"
					<<"population diversity="<<m_alg_stat.pos_diver
					<<"\n"
					<<"mean search radius="<<m_alg_stat.avg_radius
					<<"\n"
					<<"truncation number="<<trunc_count
					<<"\n"
					<<"\n"
					<<"elitist points number="<<ext_archive.size()
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

	}// end namespace mosade
}// end namespace de