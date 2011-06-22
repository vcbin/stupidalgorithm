#include "GA_Alg.h"
#include "initializer.h" 
#include "Rand_Val.h" 

using boost::timer; 
using boost::shared_ptr;
using boost::progress_display;
using boost::uniform_int;
using boost::uniform_01;
using boost::normal_distribution;
using boost::mt19937; 
using boost::variate_generator; 

using std::cout; 
using std::endl;
using std::ofstream;
using std::nth_element;
using std::random_shuffle;
using std::vector;
using std::copy;
using std::back_inserter;
//using std::make_heap;
//using std::pop_heap;

extern boost::mt19937 gen;

namespace ga
{
	namespace dgea
	{
		// generate breeding pool with binary tournament
		void dgea_alg::tour_breed(const population& pop, population& new_pop)
		{
			int *a1, *a2;
			size_t i;
			individual parent1, parent2;
			size_t pop_size=m_ppara->get_pop_size();

			a1=new int[pop_size];
			a2=new int[pop_size];
			for(i=0; i<pop_size; ++i)
			{
				a1[i]=a2[i]=i;
			}
			random_shuffle(a1, a1+pop_size);
			random_shuffle(a2, a2+pop_size);

			for(i=0; i<pop_size; i += 4)
			{
				parent1=tournament(pop[a1[i]], pop[a1[i+1]]);
				parent2=tournament(pop[a1[i+2]], pop[a1[i+3]]);
				crossover(parent1, parent2, new_pop[i], new_pop[i+1]);
				parent1=tournament(pop[a1[i]], pop[a1[i+1]]);
				parent2=tournament(pop[a1[i+2]], pop[a1[i+3]]);
				crossover(parent1, parent2, new_pop[i+2], new_pop[i+3]);
			}

			delete[] a1;
			delete[] a2;

			return;
		}// end function selection

		void dgea_alg::update_mode(int& mode,double d_l,double d_h)
		{
			if( m_alg_stat.pos_diver > d_h )
				mode=exploit;
			else if( m_alg_stat.pos_diver < d_l )
				mode=explore;
		}

		void dgea_alg::gen_child(int mode,const population& pop,population& child_pop)
		{
			if ( exploit==mode )
				tour_breed(pop, child_pop);
			else
			{
				child_pop=pop;
				mutation(child_pop);
			}
		}

		int dgea_alg::run()
		{ 	
			timer elapsed_t; 
			// retrieve algorithm parameters 
			size_t pop_size=m_ppara->get_pop_size(); 
			size_t num_dims=m_ppara->get_dim();
			double vtr=m_ppara->get_vtr();
			double d_h=m_ppara->get_dh();
			double d_l=m_ppara->get_dh();

			int m_cur_run; 
			int max_run=m_ppara->get_max_run();
			// run/trial number 
			shared_ptr<progress_display> pprog_dis;
			//// algorithm progress indicator from boost 
			//alloc_prog_indicator(pprog_dis);

			// allocate pop 	
			population pop(pop_size); 
			allocate_pop(pop,num_dims);
			population child_pop(pop_size);
			allocate_pop(child_pop,num_dims);
			// generate algorithm statistics output file name 	
			ofstream stat_file(m_com_out_path.stat_path); 
			// alloc stop condition 
			bool run_once=(1==max_run); 
			alloc_stop_cond(); 
			for(m_cur_run=0; m_cur_run<max_run; ++m_cur_run)
			{
				reset_run_stat();
				set_orig_pop(pop);
				update_diversity(pop);

				print_run_times(stat_file,m_cur_run+1);
				print_run_title(stat_file); 
				// output original population statistics 
				print_gen_stat(stat_file,1,m_alg_stat); 
				record_gen_vals(m_alg_stat,m_cur_run);
				m_cur_gen=1;
				int mode=exploit;
				while ( false==(*m_pstop_cond) )	// for every iteration 	
				{
					update_mode(mode,d_l,d_h);
					gen_child(mode,pop,child_pop);
					eval_pop(child_pop, *m_pfunc, m_alg_stat);
					select(pop,child_pop);
					stat_pop(pop, m_alg_stat);
					update_search_radius();
					update_diversity(pop);

					print_gen_stat(stat_file,m_cur_gen+1,m_alg_stat);

					record_gen_vals(m_alg_stat,m_cur_run);
					update_conv_stat(vtr);

					m_cur_gen++; 
				}// while single run termination criterion is not met 

				// single run end
				stat_run(pop,m_cur_run);// stat single run for algorithm analysis
				if ( is_final_run(m_cur_run,max_run) )
					print_run_stat(stat_file,m_alg_stat,max_run);
				/*if ( !run_once )
				++(*pprog_dis);	*/	
			}

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


		// perform (\lambda+\mu) -> (\lambda) selection
		void dgea_alg::select(population& pop, population& child_pop)
		{
			int pop_size=pop.size();
			int num_dims=m_ppara->get_dim();
			int i;
			population pop_merged(pop_size);
			vector<vector<double> > prev_x(pop_size);
			for ( i=0;i<pop_size;i++ )
			{
				prev_x[i].resize(num_dims);
				pop_merged[i]=pop[i];
				// record previous x
				prev_x[i]=pop[i].x;
			}// for every individual

			// push child population at the tail of merged population
			copy( child_pop.begin(),child_pop.end(),back_inserter(pop_merged) );

			nth_element(pop_merged.begin(), pop_merged.begin()+pop_size, pop_merged.end());
			//// heap sorting
			//make_heap(pop_merged.begin(),pop_merged.end());
			//for ( i=0;i<pop_size;i++ )
			//{
			//	// record delta x for stat purpose
			//	pop_heap(pop_merged.begin(),pop_merged.end());
			//	pop_merged.pop_back();
			//	const individual& new_ind=pop_merged.front();
			//	pop[i]=new_ind;
			//	m_alg_stat.delta_x[i]=pop[i].x-prev_x[i];// pop_merged:{parents,offsprings}
			//}

			for ( i=0;i<pop_size;i++ )
			{
				pop[i]=pop_merged[i];
				// record delta x
				m_alg_stat.delta_x[i] = (pop[i].x-prev_x[i]);
			}// for every individual
		}// end function select

		/* Routine for binary tournament */
		const individual& dgea_alg::tournament(const individual& ind1,const individual& ind2)
		{
			if ( ind1.obj[0]<ind2.obj[0] )
			{
				return (ind1);
			}
			else if ( ind1.obj[0]==ind2.obj[0] )
			{
				uniform_int<> unint_dist;
				variate_generator<mt19937, uniform_int<> > rnd_01(gen, unint_dist);
				bool ind1_flag= (rnd_01()==0?true:false);
				return ( ind1_flag?ind1:ind2 );
			}
			else 
			{
				return (ind2);
			}
		}

		/* Routine for real variable SBX crossover */
		void dgea_alg::crossover(const individual& ind1,const individual& ind2, individual& ind3, individual& ind4)
		{
			size_t i;
			int ncross=0;
			double rand;
			double y1, y2, yl, yu;
			double alpha, beta, betaq;
			double pc=m_ppara->get_pr();
			size_t ind_size=ind1.x.size();
			val_range bound=m_ppara->get_val_bnd(); 

			uniform_01<> uni01_dist;
			variate_generator<mt19937&, uniform_01<> > uni_01(gen, uni01_dist);
			if (uni_01() <= pc)
			{
				ncross++;
				for (i=0; i<ind_size; ++i)
				{
					if (uni_01() <= 0.5 )
					{
						if (ind1.x[i]<ind2.x[i])
						{
							y1 = ind1.x[i];
							y2 = ind2.x[i];
						}
						else
						{
							y1 = ind1.x[i];
							y2 = ind2.x[i];
						}
						if (fabs(ind1.x[i]-ind2.x[i]) > 1.0e-14)
						{
							yl=bound[i].min_val;
							yu=bound[i].max_val;
							rand=uni_01();

							const double eta_c=20.0;	//distribution index
							beta = 1.0 + (2.0*(y1-yl)/(y2-y1));
							alpha = 2.0 - pow(beta,-(eta_c+1.0));
							if (rand <= (1.0/alpha))
							{
								betaq = pow ((rand*alpha),(1.0/(eta_c+1.0)));
							}
							else
							{
								betaq = pow ((1.0/(2.0 - rand*alpha)),(1.0/(eta_c+1.0)));
							}
							ind3.x[i] = 0.5*((y1+y2)-betaq*(y2-y1));
							beta = 1.0 + (2.0*(yu-y2)/(y2-y1));
							alpha = 2.0 - pow(beta,-(eta_c+1.0));
							if (rand <= (1.0/alpha))
							{
								betaq = pow ((rand*alpha),(1.0/(eta_c+1.0)));
							}
							else
							{
								betaq = pow ((1.0/(2.0 - rand*alpha)),(1.0/(eta_c+1.0)));
							}
							ind4.x[i] = 0.5*((y1+y2)+betaq*(y2-y1));
							if (ind3.x[i]<yl)
							{
								ind3.x[i]=yl;
							}
							if (ind3.x[i]>yu)
							{
								ind3.x[i]=yu;
							}
							if (ind4.x[i]<yl)
							{
								ind4.x[i]=yl;
							}
							if (ind4.x[i]>yu)
							{
								ind4.x[i]=yu;
							}
						}
						else
						{
							ind3.x[i] = ind1.x[i];
							ind4.x[i] = ind2.x[i];
						}
					}
					else
					{
						ind3.x[i] = ind1.x[i];
						ind4.x[i] = ind2.x[i];
					}
				}
			}
			else
			{
				for (i=0; i<ind_size; i++)
				{
					ind3.x[i] = ind1.x[i];
					ind4.x[i] = ind2.x[i];
				}
			}
			return;
		}

		/* Function to perform mutation in a population */
		void dgea_alg::mutation(population& pop)
		{
			size_t i;
			size_t pop_size=m_ppara->get_pop_size();
			for (i=0; i<pop_size; i++)
			{
				mutation_ind(pop[i]);
			}
			return;
		}

		/* Function to perform mutation of an individual */
		void dgea_alg::mutation_ind(individual& ind)
		{
			size_t j;
			double rnd, delta1, delta2, mut_pow, deltaq;
			double delta;
			double y, yl, yu, val, xy;
			double pm=m_ppara->get_pm();
			val_range bound=m_ppara->get_val_bnd(); 
			const double eta_m=20.0;	//distribution index
			// random
			uniform_01<> uni01_dest;
			variate_generator<mt19937&, uniform_01<> > uni01_rnd(gen, uni01_dest);
			for (j=0; j<ind.x.size(); j++)
			{
				rnd = uni01_rnd();
				if (rnd <= pm)
				{
					y = ind.x[j];
					yl = bound[j].min_val;
					yu = bound[j].max_val;
					delta1 = (y-yl)/(yu-yl);
					delta2 = (yu-y)/(yu-yl);
					if(delta1<delta2)
						delta=delta1;
					else
						delta=delta2;
					rnd=uni01_rnd();
					mut_pow = 1.0/(eta_m+1.0);
					if (rnd <= 0.5)
					{
						xy = 1.0-delta1;
						val = 2.0*rnd+(1.0-2.0*rnd)*(pow(xy,(eta_m+1.0)));
						deltaq =  pow(val,mut_pow) - 1.0;
					}
					else
					{
						xy = 1.0-delta2;
						val = 2.0*(1.0-rnd)+2.0*(rnd-0.5)*(pow(xy,(eta_m+1.0)));
						deltaq = 1.0 - (pow(val,mut_pow));
					}
					y = y + deltaq*(yu-yl);
					if (y<yl)
					{
						y = yl;
					}
					if (y>yu)
					{
						y = yu;
					}
					ind.x[j] = y;
				}
			}
			return;
		}
	}// end namespace dgea
}// end namespace ea
