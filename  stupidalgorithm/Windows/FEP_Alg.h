#ifndef STUPIDALGO_FAST_EVOLUTIONARY_PROGRAMMING
#define STUPIDALGO_FAST_EVOLUTIONARY_PROGRAMMING

#include "com_alg.h"
#include "ep_para.h"

namespace ep
{
	namespace fep
	{
		struct ind_info // struct for score sorting
		{
			void set_all_prop(const individual &rhs,int scr_val) { ind=rhs;score=scr_val; }
			void alloc(int num_dims,int stra_num) { allocate_ind(ind,num_dims,stra_num); }
			bool operator<(const ind_info &rhs) { return score<rhs.score; }
			bool operator<=(const ind_info &rhs) { return score<=rhs.score; }
			bool operator>(const ind_info &rhs) { return score>rhs.score; }
			bool operator>=(const ind_info &rhs) { return score>=rhs.score; }

			individual ind;
			int score;
		};

		/* 
		Randomized quick sort routine to sort a population based on a particular objective chosen,
		in DEscending order
		*/
		template<typename Elem_Type>
		void quicksort_ind_score(std::vector<Elem_Type> &pop,std::vector<int> &idx_res)
		{
			q_sort_ind_score (pop,idx_res, 0, idx_res.size()-1);
			return;
		}

		/* Actual implementation of the randomized quick sort used to sort a population based on a particular objective chosen */
		template<typename Elem_Type>
		void q_sort_ind_score(std::vector<Elem_Type> &pop,std::vector<int> &idx_res, int left, int right)
		{
			int index;
			int temp;
			int i, j;
			Elem_Type pivot;
			if ( left<right )
			{
				uniform_int<> dist_uni(left,right);
				variate_generator<mt19937&, uniform_int<> > rnd_idx(gen, dist_uni);
				index = rnd_idx();
				temp = idx_res[right];
				idx_res[right] = idx_res[index];
				idx_res[index] = temp;
				pivot = pop[idx_res[right]];
				i = left-1;
				for ( j=left; j<right; j++ )
				{
					if ( pop[idx_res[j]] > pivot )
					{
						i+=1;
						temp = idx_res[j];
						idx_res[j] = idx_res[i];
						idx_res[i] = temp;
					}
				}
				index=i+1;
				temp = idx_res[index];
				idx_res[index] = idx_res[right];
				idx_res[right] = temp;
				q_sort_ind_score (pop, idx_res, left, index-1);
				q_sort_ind_score (pop, idx_res, index+1, right);
			}
			return;
		}

		class fep_alg:public com_alg<ep_para>
		{
		public:
			// ctor
			fep_alg(std::string conf_path):
			  com_alg(conf_path) {}
			  fep_alg(boost::shared_ptr<ep_para> ppara):// set algorithm parameters pointer directly
			  com_alg(ppara) {}
			  // dtor
			  virtual ~fep_alg() { }

			  void initialize();
			  int run();
			  inline void  set_gen_eta_file_name( std::string str_filename ){ m_avg_eta_path=str_filename; }
		protected:
			static const int stra_num=1;
			enum stra_name{eta};
			void set_tau_1();
			void set_tau_2();
			double m_tau_1;
			double m_tau_2;

			void update_mean_eta(const population& pop);

			std::string m_avg_eta_path;

			void record_dim_eta(int cur_run);
			void record_gen_vals(alg_stat &alg_stat,int cur_run);
			void calc_gen_avg_vals();
			void write_avg_eta_per_gen();
			void write_stat_vals();
			d_array m_eta_mean;
			std::vector<d_mat > m_gen_eta_val;
			d_mat m_gen_avg_eta_val;

			virtual void generate_offspring(population &parent,population &offspring);
			double generate_normal_num(double mean=0.0,double sigma=1.0);
			double generate_cauchy_num(double mean=0.0,double sigma=1.0);
			void update_eta(population& pop);
			void update_eta(double &eta,double ind_rnd,double rnd_dim);
			void merge(const population &parent_pop, const population &off_pop, population &dst_pop);
			void tour_select(const population &merged_pop,int tour_size,population &new_pop);

			bool bound_check(double &x,int dim);// SPECIAL boundaries check

			void print_run_title(std::ostream &os);
			void print_gen_stat(std::ostream &os,int cur_gen,const alg_stat &alg_stat);

			double generate_std_rnd_num();
			void initialize_eta(population& pop,double eta_val);
			std::vector<ind_info> m_sort_pop;// real population for tournament sorting
		};// end class fep_algo declaration
	}// end namespace fep
}// end namespace ep


#endif