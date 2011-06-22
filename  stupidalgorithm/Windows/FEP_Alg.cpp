#include "FEP_Alg.h"

using std::vector;
using std::cout;
using std::string;
using std::ofstream;
using std::ostream;
using std::endl;
using std::ostream;
using std::ios;

using std::make_heap;
using std::pop_heap;

using boost::shared_ptr;
using boost::timer;

using boost::uniform_int;
using boost::uniform_01;
using boost::cauchy_distribution;
using boost::normal_distribution;
using boost::mt19937;
using boost::variate_generator;

extern boost::mt19937 gen;

namespace ep
{
	namespace fep
	{
		void fep_alg::initialize()
		{
			com_alg::initialize();
			set_tau_1();
			set_tau_2();
			int num_dims=m_ppara->get_dim();
			int union_size=2*m_ppara->get_pop_size();
			m_sort_pop.reserve(union_size);
			m_sort_pop.resize(union_size);
			int k;
			for ( k=0;k<union_size;k++ )
				m_sort_pop[k].alloc(num_dims,stra_num);
			double ini_eta=m_ppara->get_ini_eta();
			int max_run=m_ppara->get_max_run();
			m_gen_eta_val.resize(max_run);
			int i;
			for ( i=0;i<max_run;i++ )
				m_gen_eta_val[i].clear();
			m_gen_avg_eta_val.resize(max_run);
		}

		void fep_alg::record_dim_eta(int cur_run)
		{
			m_gen_eta_val[cur_run].push_back(m_eta_mean);
		}

		void fep_alg::record_gen_vals(alg_stat &alg_stat,int cur_run)
		{
			problem_base::record_gen_vals(alg_stat,cur_run);
			record_dim_eta(cur_run);
		}

		void fep_alg::update_mean_eta(const population& pop)
		{
			int num_dims=m_ppara->get_dim();
			int pop_size=m_ppara->get_pop_size();
			int i,j;
			m_eta_mean.assign(num_dims,0.0);
			for ( j=0;j<num_dims;j++ ) // for every dimension
				for ( i=0;i<pop_size;i++ ) // for every individual
					m_eta_mean[j] += pop[i].stra[eta][j];
			m_eta_mean /= pop_size;
		}

		void fep_alg::calc_gen_avg_vals()
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
			int num_dims=m_ppara->get_dim();
			m_alg_stat.gen_avg_div_val.assign(min_gen,0.0);
			m_alg_stat.gen_avg_rad_val.assign(min_gen,0.0);
			m_gen_avg_eta_val.assign(min_gen,vector<double>(num_dims,0.0));
			int j,k;
			for ( j=0;j<min_gen;j++ ) // for every generation
				for ( i=0;i<max_run;i++ ) // for every run
				{
					m_alg_stat.gen_avg_div_val[j] += m_alg_stat.gen_div_val[i][j];// calculate gen avg diversity
					m_alg_stat.gen_avg_rad_val[j] += m_alg_stat.gen_rad_val[i][j];// calculate gen avg radius
					for ( k=0;k<num_dims;k++ )// calculate gen avg eta value of every dim
						m_gen_avg_eta_val[j][k] += m_gen_eta_val[i][j][k];
				}

				m_alg_stat.gen_avg_div_val /= max_run;
				m_alg_stat.gen_avg_rad_val /= max_run;
				for ( j=0;j<min_gen;j++ )
					m_gen_avg_eta_val[j] /= max_run;
		}// end function calc_gen_avg_vals

		void fep_alg::write_avg_eta_per_gen()
		{
			ofstream avg_eta_file(m_avg_eta_path);
			print_gen_avg_vec(avg_eta_file,m_gen_avg_eta_val);
		}

		void fep_alg::write_stat_vals()
		{
			com_alg::write_stat_vals();
			write_avg_eta_per_gen();
		}

		double fep_alg::generate_std_rnd_num()
		{
			int mut_type=m_ppara->get_mut_type();
			double rnd_num;
			if ( 1==mut_type )
				rnd_num=generate_cauchy_num();
			else
				rnd_num=generate_normal_num();
			return rnd_num;
		}

		void fep_alg::generate_offspring(population &parent,population &offspring)
		{
			int pop_size=m_ppara->get_pop_size();
			int num_dims=m_ppara->get_dim();
			int i;
			int j;
			for ( i=0; i<pop_size; i++ )
			{
				double ind_rnd=generate_normal_num();
				for ( j=0; j<num_dims; j++ )
				{
					double rnd_num_std=generate_std_rnd_num();
					double dim_norm_rnd=generate_normal_num();

					update_eta(offspring[i].stra[eta][j],ind_rnd,dim_norm_rnd);

					offspring[i].x[j]=parent[i].x[j]+(offspring[i].stra[eta][j]*rnd_num_std);
					// offspring[i].stra=parent[i].stra;

					bool out_bound=bound_check(offspring[i].x[j],j);
					if ( out_bound )
					{
						m_alg_stat.all_ob_num++;
						double ini_eta=m_ppara->get_ini_eta();
						offspring[i].stra[eta][j]=ini_eta;// reset eta to initial value
					}
				}// for every dimension
			}// for every individual
		}// end function generate_offspring

		/* Routine to merge two populations into one */
		void fep_alg::merge(const population &parent_pop, const population &off_pop, population &dst_pop)
		{
			int pop_size=parent_pop.size();
			int i, k;
			for ( i=0; i<pop_size; i++ )
				dst_pop[i]=parent_pop[i];
			for ( i=0, k=pop_size; i<pop_size; i++, k++ )
				dst_pop[k]=off_pop[i];
			return;
		}

		double fep_alg::generate_cauchy_num(double mean,double sigma)
		{
			// generator for random SHUFFLE VECTOR index
			cauchy_distribution<> cau_dist(mean,sigma);
			variate_generator<mt19937&, cauchy_distribution<> > cauc_rnd(gen, cau_dist);
			return cauc_rnd();
		}

		double fep_alg::generate_normal_num(double mean,double sigma)
		{
			// generator for random SHUFFLE VECTOR index
			normal_distribution<> norm_dist(mean,sigma);
			variate_generator<mt19937&, normal_distribution<> > norm_rnd(gen, norm_dist);
			return norm_rnd();
		}

		void fep_alg::set_tau_1() 
		{
			int num_dims=m_ppara->get_dim();
			m_tau_1= 1.0/sqrt(2.0*sqrt(num_dims*1.0));
		}

		void fep_alg::set_tau_2() 
		{
			int num_dims=m_ppara->get_dim();
			m_tau_2= 1.0/sqrt(2.0*num_dims);
		}

		void fep_alg::update_eta(double &eta,double ind_rnd,double rnd_dim)
		{
			double tmp=exp(m_tau_2*ind_rnd+m_tau_1*rnd_dim);
			eta=eta*tmp;
			/*if ( eta>1.0 )
			eta=1.0;*/
		}

		// OBSOLETE function
		void fep_alg::update_eta(population& pop)
		{
			int pop_size=m_ppara->get_pop_size();
			int num_dims=m_ppara->get_dim();
			int i;
			int j;
			for ( i=0; i<pop_size; i++ )
			{
				double ind_rnd=generate_normal_num();
				for ( j=0; j<num_dims; j++ )
				{
					double rnd_dim=generate_normal_num();
					update_eta(pop[i].stra[eta][j],ind_rnd,rnd_dim);
				}
			}
		}// end function update_eta

		void fep_alg::tour_select(const population &merged_pop,int tour_size,population &new_pop)
		{
			int union_size=merged_pop.size();
			int pop_size=union_size/2;
			int k,l;

			// generator for random merged population index
			uniform_int<> pop_rnd(0,union_size-1);
			variate_generator<mt19937&, uniform_int<> > rnd_pop_idx(gen, pop_rnd);
			// tournament selection
			int num_dims=m_ppara->get_dim();

			// tournament ranking
			vector<int> ind_score(union_size,0);
			for ( k=0;k<union_size;k++ )
			{
				for ( l=0;l<tour_size;l++ )
				{
					if ( merged_pop[k] <= merged_pop[rnd_pop_idx()] )
						ind_score[k]++;
				}
			}

			m_sort_pop.resize(union_size);
			//for ( k=pop_size;k<union_size;k++ ) // previous pop_heap operation reduced m_sort_pop size,reallocate to fit
			//{
			//	m_sort_pop[k].alloc(num_dims,stra_num);
			//	m_sort_pop[k].ind.stra.resize(num_dims);
			//}

			//double ini_eta=m_ppara->get_ini_eta();
			for ( k=0;k<union_size;k++ )
			{
				m_sort_pop[k].set_all_prop(merged_pop[k],ind_score[k]);
				/*if ( k<pop_size )
				m_sort_pop[k].set_all_prop(merged_pop[k],m_eta[k],ind_score[k]);
				else
				m_sort_pop[k].set_all_prop(merged_pop[k],vector<double>(num_dims,ini_eta),ind_score[k]);*/
			}

			int i;
			// SLOWER than heap sort
			//vector<int> sort_res_idx(union_size);
			//for ( i=0;i<union_size;i++ )
			//	sort_res_idx[i]=i;
			// // sort the merged pop according to individual score
			// quicksort_ind_score(m_sort_pop,sort_res_idx);
			//
			//for ( i=0;i<pop_size;i++ )
			//{
			//	// record delta x for stat purpose
			//	new_pop[i]=merged_pop[sort_res_idx[i]];
			//	m_alg_stat.delta_x[i]=new_pop[i].x-merged_pop[i].x;
			//	
			//}

			// heap sorting
			make_heap(m_sort_pop.begin(),m_sort_pop.end());
			for ( i=0;i<pop_size;i++ )
			{
				// record delta x for stat purpose
				pop_heap(m_sort_pop.begin(),m_sort_pop.end());
				m_sort_pop.pop_back();
				const ind_info& new_ind=m_sort_pop.front();
				new_pop[i]=new_ind.ind;
				m_alg_stat.delta_x[i]=new_pop[i].x-merged_pop[i].x;// merged_pop:{parents,offsprings}
			}
			return;
		}// end function tour_select

		// SPECIAL boundaries check
		bool fep_alg::bound_check(double &x,int dim)
		{
			const val_range& val_bounds=m_ppara->get_val_bnd();
			bool out_up_bnd,out_low_bnd;
			out_low_bnd=x < val_bounds[dim].min_val;
			out_up_bnd=x > val_bounds[dim].max_val;
			double low_bnd,high_bnd;
			low_bnd=val_bounds[dim].min_val;
			high_bnd=val_bounds[dim].max_val;
			if ( out_up_bnd || out_low_bnd )
			{
				//x=m_alg_stat.cen_ind[dim]+generate_normal_num();
				//bound_check(x,dim);// RECURSIVE invocation

				// uniformly randomize x within [lower_bound,upper_bound]
				uniform_01<> dist;
				variate_generator<mt19937&, uniform_01<> > rnd_num(gen, dist);
				x=low_bnd+rnd_num()*(high_bnd-low_bnd);
				bound_check(x,dim);// RECURSIVE invocation
				return true;
			}
			return false;
		}// end function bound_check

		void fep_alg::print_run_title(ostream &os)
		{
			output_base::print_run_title(os);
			os<<"\t"
				<<"\t"
				<<"avg_radius";
		}// end function print_run_title

		void fep_alg::print_gen_stat(ostream &os,int cur_gen,const alg_stat &alg_stat)
		{
			output_base::print_gen_stat(os,cur_gen,alg_stat);
			// os.setf(ios::fixed,ios::floatfield);
			os.precision(8);
			os<<"\t"
				<<"\t"
				<<m_alg_stat.avg_radius;
		}// end function print_gen_stat

		void fep_alg::initialize_eta(population& pop,double eta_val)
		{
			// initialize 2D vector of eta value
			int pop_size=m_ppara->get_pop_size();
			int num_dims=m_ppara->get_dim();
			double ini_eta=m_ppara->get_ini_eta();
			int i;
			for ( i=0;i<pop_size;i++ )
				pop[i].stra[eta].assign(num_dims,eta_val);
		}

		int fep_alg::run()
		{
			if ( !m_ppara )
				return -1;

			timer elapsed_t;
			// retrieve algorithm parameters
			int pop_size=m_ppara->get_pop_size();
			int num_dims=m_ppara->get_dim();
			double vtr=m_ppara->get_vtr();
			int max_gen=m_ppara->get_max_gen();
			max_gen=max_gen/2;// ifep specific
			int tour_size=m_ppara->get_tour_size();
			double ini_eta=m_ppara->get_ini_eta();

			int m_cur_run;
			int max_run=m_ppara->get_max_run();// run/trial number
			// shared_ptr<progress_display> pprog_dis;// algorithm progress indicator from boost
			// alloc_prog_indicator(pprog_dis);

			// allocate original pop and trial pop
			population pop(pop_size);
			allocate_pop(pop,num_dims,stra_num);
			population off_pop(pop_size);
			allocate_pop(off_pop,num_dims,stra_num);
			int union_size=2*pop_size;
			population uni_pop(union_size);
			allocate_pop(uni_pop,num_dims,stra_num);
			reallocate_stra(pop,off_pop,eta,num_dims);
			reallocate_stra(uni_pop,eta,num_dims);

			// generate algorithm statistics output file name
			ofstream stat_file(m_com_out_path.stat_path);
			// allocate stop condition object dynamically
			bool run_once=(1==max_run);
			alloc_stop_cond();

			// iteration start
			for ( m_cur_run=0;m_cur_run<max_run;m_cur_run++ )
			{
				reset_run_stat();
				initialize_eta(pop,ini_eta);
				initialize_eta(off_pop,ini_eta);
				m_eta_mean.assign(num_dims,ini_eta);

				set_orig_pop(pop);
				update_diversity(pop);

				record_gen_vals(m_alg_stat,m_cur_run);

				print_run_times(stat_file,m_cur_run+1);
				print_run_title(stat_file);
				// output original population statistics
				print_gen_stat(stat_file,1,m_alg_stat);

				m_cur_gen=1;
				while ( false==(*m_pstop_cond) ) // for every iteration
				{
					generate_offspring(pop,off_pop);
					// evaluate offspring population
					eval_pop(off_pop,*m_pfunc,m_alg_stat);
					merge(pop,off_pop,uni_pop);
					tour_select(uni_pop,tour_size,pop);
					stat_pop(pop,m_alg_stat);

					update_mean_eta(pop);
					update_search_radius();
					update_diversity(pop);

					record_gen_vals(m_alg_stat,m_cur_run);

					print_gen_stat(stat_file,m_cur_gen+1,m_alg_stat);
					update_conv_stat(vtr);

					/*if ( run_once )
					++(*pprog_dis);*/

					m_cur_gen++;
				}// while single run termination criterion is not met

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
		// end class fep_algo definition
	}// namespace fep
}// end namespace ep