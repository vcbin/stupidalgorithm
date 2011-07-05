#include "ifep_alg.h"

using std::rand;
using std::vector;
using std::cout;
using std::string;
using std::ofstream;
using std::ostream;
using std::endl;
using std::ostream;
using std::ios;
using std::ifstream;
using std::getline;
using std::stringstream;
using std::logic_error;

using boost::shared_ptr;
using boost::timer;

using boost::uniform_01;
using boost::mt19937;
using boost::variate_generator;

extern boost::mt19937 gen;

namespace ep
{
	namespace ifep
	{
		int ifep_alg::load_ini_pop(population &pop,std::string ini_pop_path)
		{
			int pop_size=m_ppara->get_pop_size();
			ifstream ini_pop_data(ini_pop_path.c_str());
			int num_dims=m_ppara->get_dim();
			int i=0;// item read number
			int j;
			string inbuf;
			vector<double> x_tmp(num_dims);
			while ( !ini_pop_data.eof() )
			{
				getline(ini_pop_data,inbuf);
				if ( false==is_dataLine(inbuf) )
					continue;
				if ( i>=pop_size ) 
					break;
				stringstream line(inbuf);
				for ( j=0;j<num_dims;j++ )
				{
					line>>x_tmp[j];
					bound_check(x_tmp[j],j);
					pop[i].x[j]=x_tmp[j];
				}
				if ( ini_pop_data.fail() ) // read error
				{
					stringstream str_err("ERROR:Read data from file");
					str_err<<" "
						<<ini_pop_path.c_str()
						<<" "
						<<"failed."
						<<"Check your data file for data validity."
						<<"\n";
					throw logic_error(str_err.str());
				}
				i++;
			}// while !eof
			// ALLOW entry count in ini_pop_file inconsistent with pop_size in config file,ifep SPECIFIC
			return 0;
		}// end function load_orig_pop

		void ifep_alg::generate_offspring(population &parent,population &offspring)
		{
			int pop_size=m_ppara->get_pop_size();
			int num_dims=m_ppara->get_dim();
			int i;
			int j;
			// allocate space for two candidate offspring
			individual offspring1,offspring2;
			allocate_ind(offspring1,num_dims,stra_num);
			allocate_ind(offspring2,num_dims,stra_num);
			offspring1.stra[eta].resize(num_dims);
			offspring2.stra[eta].resize(num_dims);
			for ( i=0; i<pop_size; i++ )
			{
				double ind_rnd=generate_normal_num();
				bool spr1_ob,spr2_ob;
				int spr1_ob_num=0,spr2_ob_num=0;
				for ( j=0; j<num_dims; j++ )
				{
					double dim_cauc_rnd=generate_cauchy_num();
					double dim_norm_rnd=generate_normal_num();

					update_eta(offspring[i].stra[eta][j],ind_rnd,dim_norm_rnd);

					offspring1.x[j]=parent[i].x[j]+(offspring[i].stra[eta][j]*dim_cauc_rnd);
					offspring1.stra[eta][j]=offspring[i].stra[eta][j];
					// boundary checking of X and eta
					spr1_ob=bound_check(offspring1.x[j],j);
					if ( spr1_ob )
					{
						spr1_ob_num++;
						double ini_eta=m_ppara->get_ini_eta();
						offspring1.stra[eta][j]=ini_eta;
					}

					offspring2.x[j]=parent[i].x[j]+(offspring[i].stra[eta][j]*dim_norm_rnd);
					offspring2.stra[eta][j]=offspring[i].stra[eta][j];
					// boundary checking of X and eta
					spr2_ob=bound_check(offspring2.x[j],j);
					if ( spr2_ob )
					{
						spr2_ob_num++;
						double ini_eta=m_ppara->get_ini_eta();
						offspring2.stra[eta][j]=ini_eta;
					}
				}// for every dimension
				// evaluate two candidate offsprings
				eval_ind(offspring1,*m_pfunc,m_alg_stat);
				eval_ind(offspring2,*m_pfunc,m_alg_stat);
				// choose the better one between two offspring
				if ( offspring1 < offspring2 )
				{
					offspring[i]=offspring1;
					m_alg_stat.all_ob_num += spr1_ob_num;
				}
				else if ( offspring1 == offspring2 )
				{
					uniform_01<> uni_01;
					variate_generator<mt19937&, uniform_01<> > rnd_01(gen, uni_01);
					bool sel_1=rnd_01()<=0.5;
					offspring[i]=( sel_1 ? offspring1:offspring2 );// uniform random selection
					m_alg_stat.all_ob_num += ( sel_1 ? spr1_ob_num:spr2_ob_num );
				}// if offspring1==offspring2
				else
				{
					offspring[i]=offspring2;
					m_alg_stat.all_ob_num += spr2_ob_num;
				}
			}// for every individual
		}// end function generate_offspring

		int ifep_alg::run()
		{
			if ( !m_ppara )
				return -1;

			timer elapsed_t;
			// retrieve algorithm parameters
			int pop_size=m_ppara->get_pop_size();
			int num_dims=m_ppara->get_dim();
			int tour_size=m_ppara->get_tour_size();
			double ini_eta=m_ppara->get_ini_eta();
			double vtr=m_ppara->get_vtr();

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
			int i;
			for ( i=0;i<pop_size;i++ )
			{
				pop[i].stra[eta].resize(num_dims);
				off_pop[i].stra[eta].resize(num_dims);
			}
			for ( i=0;i<union_size;i++ )
				uni_pop[i].stra[eta].resize(num_dims);

			// generate algorithm statistics output file name
			ofstream stat_file(m_com_out_path.stat_path.c_str());
			// allocate stop condition object dynamically
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
					// NO need to evaluate offspring population
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

	}// end namespace ifep
}// end namespace ep
