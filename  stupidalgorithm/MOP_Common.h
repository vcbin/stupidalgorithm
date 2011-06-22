#ifndef JERRYLIU_MOP_PUBLIC_INTERFACE
#define JERRYLIU_MOP_PUBLIC_INTERFACE

#include "Algo_DataStruct.h"
//#include "plot.h"
#include <list>
#include <boost/shared_array.hpp>
#include <boost/thread.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/lexical_cast.hpp>
#include <string>
//#include <deque>

namespace mop
{
	typedef std::list<individual> elite_archive;

	enum comp{worse=-1,better=1,indiff=0};
	enum {invalid=-1};
	enum {first,second,third};// objective
	enum {y1,y2,y3};

	void q_sort_front_obj(const population& pop, int obj_idx,std::vector<int>& rank, int left, int right);
	void quicksort_ind_score(const population& pop, int obj_idx,std::vector<int>& rank, int rank_size);
	void quicksort_dist(const population &pop,idx_array &d_idx, int front_size);
	void q_sort_dist(const population &pop,idx_array &d_idx, int left, int right);

	typedef void (*pplot_fun)(boost::shared_ptr<boost::mutex> pmut,int gen_idx,int arc_size,
		int stag_gen,bool max_plot_gen,std::string cmd_file);
	void plot_fun(boost::shared_ptr<boost::mutex> pmut,int gen_idx,int arc_size,
		int stag_gen,bool max_plot_gen,std::string cmd_file);

	class front_pred:public std::unary_function<individual,bool>
	{
	public:
		bool operator()(individual &ind);
	};
	class first_obj_le:public std::binary_function<d_array,d_array,bool>
	{
	public:
		bool operator() (const d_array &lhs,const d_array &rhs);
	};

	template <typename ForwardIterator>
	void copy_obj(ForwardIterator first,ForwardIterator last,d_mat &coor)
	{
		ForwardIterator cur_itr;
		for ( cur_itr=first;cur_itr!=last;++cur_itr )
			coor.push_back(cur_itr->obj);
	}
	template <typename ForwardIterator,typename Pred>
	void copy_obj_if(ForwardIterator first,ForwardIterator last,d_mat &coor,Pred pred)
	{
		ForwardIterator cur_itr;
		for ( cur_itr=first;cur_itr!=last;++cur_itr )
		{
			if ( pred(*cur_itr) )
				coor.push_back(cur_itr->obj);
		}
	}

	bool load_pop(std::string path,int col_num,d_mat &coor);

	class mop_common
	{
	public:
		int check_dominance(const individual &a,const individual &b);
		bool update_archive(elite_archive &ext_archive,const individual &ind);
		// wholesome truncation
		void trunc_external(elite_archive &ext_archive,int trunc_type,int max_size,int &trunc_count,int neighbor_num);
		void trunc_archive_cd(elite_archive &ext_archive,int max_size,int &trunc_count,int neighbor_num=2);
		void trunc_archive_entropy(elite_archive &ext_archive,int max_size,int &trunc_count);
		void trunc_archive_harmonic(elite_archive &ext_archive,int max_size,
			int &trunc_count,int neighbor_num=2);
		// step-by-step calculation and truncation,precise and TIME-CONSUMING
		void trunc_archive_entropy_exact(elite_archive &ext_archive,int max_size,int &trunc_count);
		void trunc_archive_harmonic_exact(elite_archive &ext_archive,int max_size,
			int &trunc_count,int neighbor_num=2);

		void calc_crowd_dist(const elite_archive &ext_archive,individual &ind,int neighbor_num=2);// neighbor_num argument IGNORED
		void calc_crowd_entropy(const elite_archive &ext_archive,individual &ind,int neighbor_num=2);// neighbor_num argument IGNORED
		void calc_crowd_harmonic(const elite_archive &ext_archive,individual &ind,int neighbor_num=2);

		void assign_crowd_indices(population &pop,const idx_array& idx_arr,int trunc_type);
		void assign_crowd_range(population &pop,int left,int right,int trunc_type);
		typedef std::list<int> front_list;
		void assign_crowd_list(population &pop,const front_list &list,int front_size,int trunc_type);

		void assign_crowd_dist_entropy(population &pop,const idx_array &idx_arr,int assi_size);
		bool assign_crowd_dist_harmonic(population &pop,const idx_array &idx_arr,int assi_size,int neighbor_num=2);
		void assign_crowd_dist(population &pop,const idx_array &idx_arr,int assi_size);
		void crowding_fill_list(population &mixed_pop,population &new_pop,int pop_size,
			int begin,int front_size,const front_list &f_idx,
			int trunc_type);
		void crowding_fill_indices(population &mixed_pop,population &new_pop,int pop_size,
			int begin,int front_size,idx_array &idx_arr,
			int trunc_type);

		void fill_nondominated_sort (population &mixed_pop, population &new_pop,int pop_size,int trunc_type);
		template <typename ForwardIterator>
		void output_collection(std::ostream &os ,const ForwardIterator first,const ForwardIterator last)
		{
			ForwardIterator itr_cur=first;
			int num_obj=first->obj.size();
			for ( ;itr_cur!=last;itr_cur++ )
			{
				os<<"\n";
				int i;
				for ( i=0;i<num_obj;i++ )
				{
					os<<itr_cur->obj[i];
					if ( i!=num_obj-1 )
						os<<" ";
				}
			}// for every member
		}// end function output_collection

		template <typename ForwardIterator,typename Pred>
		void output_if(std::ostream &os ,const ForwardIterator first,const ForwardIterator last,Pred pred)
		{
			ForwardIterator itr_cur=first;
			int num_obj=first->obj.size();
			for ( ;itr_cur!=last;itr_cur++ )
			{
				if ( pred(*itr_cur) )
				{
					os<<"\n";
					int i;
					for ( i=0;i<num_obj;i++ )
					{
						os<<itr_cur->obj[i];
						if ( i!=num_obj-1 )
							os<<" ";
					}
				}
			}// for every member
		}// end function output_frontier

		struct perf_indice// performance indicator
		{
			double gamma;// convergence indicator
			double delta;// spread,diversity indicator
			double hv;// hyper-volume
		};
		struct point
		{
			point(double x1_val,double x2_val):
		y1(x1_val),y2(x2_val){}
		double y1;
		double y2;
		};
		bool is_final_out_gen(int gen_idx,int out_interval,int max_gen);
		double zdt3_calc_gamma(d_mat &opt_coor,int sample_size);
		void zdt3_assess(const d_mat &coor,int sample_size,const point& ref_point,perf_indice &p_ind);

	};// end class mop_public

}// end namespace mop

#endif