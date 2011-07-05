#ifdef WIN32
#define _USE_MATH_DEFINES // for math constant eg. pi,e...
#endif

#include <cmath>

#include "mop_common.h"
#include "rand_val.h"
#include "cmdline_para_type.h"

//#include <boost/lambda/lambda.hpp>
#include <string>
#include <fstream>
#include <sstream>
#include <iterator>
// #include <boost/type_traits.hpp>
// #include <boost/assign.hpp>
using std::list;
using std::vector;
using std::copy;
using std::string;
using std::back_inserter;

using std::numeric_limits;
using std::stringstream;
using std::ifstream;
using std::logic_error;
using std::nth_element;

using boost::variate_generator;
using boost::uniform_int;
using boost::uniform_real;
using boost::mt19937;
using boost::lock_guard;
using boost::shared_ptr;
using boost::lexical_cast;
using std::for_each;
using boost::bind;
extern boost::mt19937 gen;

using namespace mop;

namespace mop
{
	void fwd_plot_fun(shared_ptr<boost::mutex> pmut,int gen_idx,int arc_size,
		int stag_gen,bool max_plot_gen,string cmd_file)
	{
		lock_guard<boost::mutex> lock(*pmut);
		string gp_cmd;
		gp_cmd="gnuplot -e \"gen=\'"+lexical_cast<string>(gen_idx+1)+"\'"
			+";"+"arc_size=\'"+lexical_cast<string>(arc_size)+"\'"+";"
			+"stag_gen=\'"+lexical_cast<string>(stag_gen)+"\'"
			+";"+"max_gen="+lexical_cast<string>(max_plot_gen)+"\""
			+" "+cmd_file;
		system(gp_cmd.c_str());
	}// end function fwd_plot_fun

	bool front_pred::operator()(individual &ind) { return ind.rank==1; }

	bool first_obj_le::operator() (const d_array &lhs,const d_array &rhs) { return lhs[first]<=rhs[first]; }

	bool load_pop(std::string path,int col_num,d_mat &coor)
	{
		ifstream opt_coor(path.c_str());
		int i=0;// entry read number
		int j;
		string inbuf;
		d_array x_tmp(col_num);
		if ( opt_coor.fail() )
		{
			stringstream str_err;
			str_err<<"\nopen file "
				<<path
				<<" ERROR!\n"
				<<"please check file path!\n";
			throw logic_error(str_err.str());
		}
		while ( !opt_coor.eof() )
		{
			getline(opt_coor,inbuf);
			if ( false==is_dataLine(inbuf) )
				continue;

			stringstream line(inbuf);
			for ( j=0;j<col_num;j++ )
			{
				line>>x_tmp[j];
			}

			coor.push_back(x_tmp);
			if ( opt_coor.fail() ) // read error
			{
				stringstream str_err("ERROR:Read data from file");
				str_err<<" "
					<<path.c_str()
					<<" "
					<<"failed."
					<<"Check your data file for data validity."
					<<"\n";
				throw logic_error(str_err.str());
			}
			i++;
		}// while !eof
		return true;
	}// load_pop

	int mop_common::check_dominance(const individual &a,const individual &b)
	{
		int i;
		int flag1;
		int flag2;
		flag1 = 0;
		flag2 = 0;
		if ( a.constr_viol<0 && b.constr_viol<0 )
		{
			if ( a.constr_viol > b.constr_viol )
			{
				return (1);
			}
			else
			{
				if ( a.constr_viol < b.constr_viol )
				{
					return (-1);
				}
				else
				{
					return (0);
				}
			}
		}
		else
		{
			if ( a.constr_viol < 0 && b.constr_viol == 0 )
			{
				return (-1);
			}
			else
			{
				if ( a.constr_viol == 0 && b.constr_viol <0 )
				{
					return (1);
				}
				else
				{
					int num_obj=a.obj.size();
					for (i=0; i<num_obj; i++)
					{
						if (a.obj[i] < b.obj[i])
						{
							flag1 = 1;

						}
						else
						{
							if (a.obj[i] > b.obj[i])
							{
								flag2 = 1;
							}
						}
					}
					if ( flag1==1 && flag2==0 )
					{
						return (1);
					}
					else
					{
						if ( flag1==0 && flag2==1 )
						{
							return (-1);
						}
						else
						{
							return (0);
						}
					}
				}
			}
		}
	}// end function check_dominance

	/* Actual implementation of the randomized quick sort used to sort a population based on a particular objective chosen */
	// Ascending order
	void q_sort_front_obj(const population& pop,int obj_idx,idx_array& rank,int obj_left,int obj_right)
	{
		int index;
		int temp;
		int i, j;
		double pivot;
		if (obj_left<obj_right)
		{
			uniform_int<> dist_uni(obj_left,obj_right);
			variate_generator<mt19937&, uniform_int<> > rnd_idx(gen, dist_uni);
			index = rnd_idx();
			temp = rank[obj_right];
			rank[obj_right] = rank[index];
			rank[index] = temp;
			pivot =pop[rank[obj_right]].obj[obj_idx];
			i = obj_left-1;
			for (j=obj_left; j<obj_right; j++)
			{
				if (pop[rank[j]].obj[obj_idx] <= pivot)
				{
					i+=1;
					temp = rank[j];
					rank[j] = rank[i];
					rank[i] = temp;
				}
			}
			index=i+1;
			temp = rank[index];
			rank[index] = rank[obj_right];
			rank[obj_right] = temp;
			q_sort_front_obj (pop, obj_idx, rank, obj_left, index-1);
			q_sort_front_obj (pop, obj_idx, rank, index+1, obj_right);
		}
		return;
	}

	/* Randomized quick sort routine to sort a population based on a particular objective chosen */
	void quicksort_ind_score(const population& pop,int obj_idx,idx_array& rank,int rank_size)
	{
		q_sort_front_obj (pop, obj_idx, rank, 0, rank_size-1);
		return;
	}

	/* Actual implementation of the randomized quick sort used to sort a population based on a particular objective chosen */
	// Ascending order
	void q_sort_array(const d_array& vec,idx_array &rank,int left,int right)
	{
		int index;
		int temp;
		int i, j;
		double pivot;
		if (left<right)
		{
			uniform_int<> dist_uni(left,right);
			variate_generator<mt19937&, uniform_int<> > rnd_idx(gen, dist_uni);
			index = rnd_idx();
			temp = rank[right];
			rank[right] = rank[index];
			rank[index] = temp;
			pivot =vec[rank[right]];
			i = left-1;
			for (j=left;j<right; j++)
			{
				if (vec[rank[j]] <= pivot)
				{
					i+=1;
					temp = rank[j];
					rank[j] = rank[i];
					rank[i] = temp;
				}
			}
			index=i+1;
			temp = rank[index];
			rank[index] = rank[right];
			rank[right] = temp;
			q_sort_array(vec,rank,left,index-1);
			q_sort_array (vec,rank,index+1,right);
		}
		return;
	}// end function q_sort_array

	/* Randomized quick sort routine to sort a population based on a particular objective chosen */
	void quicksort_array(const d_array& vec,idx_array &rank,int size)
	{
		q_sort_array (vec,rank,0,size-1);
		return;
	}

	/* Actual implementation of the randomized quick sort used to sort a population based on crowding distance */
	void q_sort_dist(const population &pop,idx_array &d_idx,int left,int right)
	{
		int index;
		int temp;
		int i, j;
		double pivot;
		if (left<right)
		{
			uniform_int<> dist_uni(left,right);
			variate_generator<mt19937&, uniform_int<> > rnd_idx(gen, dist_uni);
			index = rnd_idx();
			temp = d_idx[right];
			d_idx[right] = d_idx[index];
			d_idx[index] = temp;
			pivot = pop[d_idx[right]].crowd_dist;
			i = left-1;
			for (j=left; j<right; j++)
			{
				if (pop[d_idx[j]].crowd_dist <= pivot)
				{
					i+=1;
					temp = d_idx[j];
					d_idx[j] = d_idx[i];
					d_idx[i] = temp;
				}
			}
			index=i+1;
			temp = d_idx[index];
			d_idx[index] = d_idx[right];
			d_idx[right] = temp;
			q_sort_dist (pop, d_idx, left, index-1);
			q_sort_dist (pop, d_idx, index+1, right);
		}
		return;
	}

	/* Randomized quick sort routine to sort a population based on crowding distance */
	void quicksort_dist(const population &pop,idx_array &d_idx,int front_size)
	{
		q_sort_dist (pop, d_idx, 0, front_size-1);
		return;
	}

	/* Actual implementation of the randomized quick sort used to sort a population based on crowding distance */
	template <typename Pred>
	void q_sort_dmat(const d_mat &vec,idx_array &d_idx,int left,int right,Pred le_pred)
	{
		int index;
		int temp;
		int i, j;
		d_array pivot;
		if ( left<right )
		{
			uniform_int<> dist_uni(left,right);
			variate_generator<mt19937&, uniform_int<> > rnd_idx(gen, dist_uni);
			index = rnd_idx();
			temp = d_idx[right];
			d_idx[right] = d_idx[index];
			d_idx[index] = temp;
			pivot = vec[d_idx[right]];
			i = left-1;
			for (j=left; j<right; j++)
			{
				if ( le_pred(vec[d_idx[j]],pivot) )
				{
					i+=1;
					temp = d_idx[j];
					d_idx[j] = d_idx[i];
					d_idx[i] = temp;
				}
			}
			index=i+1;
			temp = d_idx[index];
			d_idx[index] = d_idx[right];
			d_idx[right] = temp;
			q_sort_dmat(vec, d_idx, left, index-1,le_pred);
			q_sort_dmat(vec, d_idx, index+1, right,le_pred);
		}
		return;
	}

	/* Randomized quick sort routine to sort a population based on crowding distance */
	template <typename Pred>
	void quicksort_dmat(const d_mat &vec,idx_array &d_idx,int front_size,Pred le_pred)
	{
		q_sort_dmat(vec, d_idx, 0, front_size-1,le_pred);
		return;
	}

	bool mop_common::is_final_out_gen(int gen_idx,int out_interval,int max_gen)
	{
		int cur_gen=gen_idx+1;
		if ( cur_gen==max_gen ) return true;
		else if ( cur_gen<max_gen && (cur_gen+out_interval)>max_gen ) return true;
		else return false;
	}

	bool mop_common::update_archive(elite_archive &ext_archive,const individual &ind)
	{
		elite_archive::iterator itr_elite=ext_archive.begin();
		bool erase=false;
		int comp_res;
		while ( itr_elite!=ext_archive.end() )
		{
			comp_res=check_dominance(*itr_elite,ind);
			if ( worse==comp_res )
				erase=true;
			else if ( better==comp_res )
				return false;
			if ( !erase )
				++itr_elite;
			else
			{
				itr_elite=ext_archive.erase(itr_elite);
				erase=false;
			}
		}
		ext_archive.push_back(ind);
		return true;
	}// end function update_archive

	void mop_common::assign_crowd_list(population &pop,const front_list &lst,int front_size,int trunc_type)
	{
		front_list::const_iterator cur;
		cur = lst.begin();
		if (front_size==1)
		{
			pop[*cur].crowd_dist = INF;
			return;
		}
		if (front_size==2)
		{
			pop[*cur].crowd_dist = INF;
			pop[*(++cur)].crowd_dist = INF;
			return;
		}

		idx_array d_idx;
		copy(lst.begin(),lst.end(),back_inserter<idx_array>(d_idx));
		switch (trunc_type)
		{
		case crowd_dist:
			{
				assign_crowd_dist(pop,d_idx,front_size);
				break;
			}
		case crowd_harmonic:
			{
				assign_crowd_dist_harmonic(pop,d_idx,front_size);
				break;
			}
		default:
			assign_crowd_dist_entropy(pop,d_idx,front_size);
			break;
		}
	}// end function assign_crowd_list

	void mop_common::assign_crowd_indices(population &pop,const idx_array& idx_arr,int trunc_type)// non-sequential indices
	{
		int front_size=idx_arr.size();

		if (front_size==1)
		{
			pop[idx_arr[0]].crowd_dist = INF;
			return;
		}
		if (front_size==2)
		{
			pop[idx_arr[0]].crowd_dist = INF;
			pop[idx_arr[1]].crowd_dist = INF;
			return;
		}
		switch (trunc_type)
		{
		case crowd_dist:
			{
				assign_crowd_dist(pop,idx_arr,front_size);
				break;
			}
		case crowd_harmonic:
			{
				assign_crowd_dist_harmonic(pop,idx_arr,front_size);
				break;
			}
		default:
			assign_crowd_dist_entropy(pop,idx_arr,front_size);
			break;
		}
		return;
	}// end function assign_crowd_indices

	void mop_common::assign_crowd_range(population &pop,int left,int right,int trunc_type)// sequential range
	{
		int front_size;
		front_size = right-left+1;

		if (front_size==1)
		{
			pop[left].crowd_dist = INF;
			return;
		}
		if (front_size==2)
		{
			pop[left].crowd_dist = INF;
			pop[right].crowd_dist = INF;
			return;
		}
		idx_array d_idx(front_size);
		init_idx_array(d_idx,front_size,left);
		switch (trunc_type)
		{
		case crowd_dist:
			{
				assign_crowd_dist(pop,d_idx,front_size);
				break;
			}
		case crowd_harmonic:
			{
				assign_crowd_dist_harmonic(pop,d_idx,front_size);
				break;
			}
		default:
			assign_crowd_dist_entropy(pop,d_idx,front_size);
			break;
		}
		return;
	}// end function assign_crowd_indices

	void mop_common::crowding_fill_list (population &mixed_pop,population &new_pop,int pop_size,
		int begin,int front_size,const front_list &f_idx,
		int trunc_type)
	{
		assign_crowd_list(mixed_pop,f_idx,front_size,trunc_type);
		idx_array d_idx;
		copy(f_idx.begin(),f_idx.end(),back_inserter<idx_array>(d_idx));
		quicksort_dist(mixed_pop,d_idx,front_size);
		int i, j;
		for ( i=begin,j=front_size-1; i<pop_size; i++, j-- )
			new_pop[i]=mixed_pop[d_idx[j]];
	}// end function crowding_fill_list

	void mop_common::crowding_fill_indices(population &mixed_pop,population &new_pop,int pop_size,
		int begin,int front_size,idx_array &idx_arr,
		int trunc_type)
	{
		assign_crowd_indices(mixed_pop,idx_arr,trunc_type);

		quicksort_dist(mixed_pop,idx_arr,front_size);
		int i, j;
		for ( i=begin,j=front_size-1; i<pop_size; i++, j-- )
			new_pop[i]=mixed_pop[idx_arr[j]];
	}// end function crowding_fill_indices

	/* Routine to compute crowding distances */
	void mop_common::assign_crowd_dist(population &pop,const idx_array &idx_arr,int front_size)
	{
		int i, j;
		int num_obj=pop[0].obj.size();
		const int last=front_size-1;
		idx_mat rank(num_obj,idx_array(front_size));// sorting order index of every objective
		for ( i=0; i<num_obj; i++ )
		{
			rank[i]=idx_arr;
			quicksort_ind_score(pop,i,rank[i],front_size);
		}
		for ( j=0; j<front_size; j++ )
			pop[idx_arr[j]].crowd_dist = 0.0;
		for ( i=0; i<num_obj; i++ )
			pop[rank[i][first]].crowd_dist = INF;
		for (i=0; i<num_obj; i++)
		{
			for (j=1; j<front_size-1; j++)
			{
				if (pop[rank[i][j]].crowd_dist != INF)
				{
					if (pop[rank[i][last]].obj[i] == pop[rank[i][first]].obj[i])
					{
						pop[rank[i][j]].crowd_dist += 0.0;
					}
					else
					{
						pop[rank[i][j]].crowd_dist += (pop[rank[i][j+1]].obj[i] - pop[rank[i][j-1]].obj[i])/
							(pop[rank[i][last]].obj[i] - pop[rank[i][first]].obj[i]);
					}
				}
			}
		}
		for ( j=0; j<front_size; j++ )
		{
			if (pop[idx_arr[j]].crowd_dist != INF)
			{
				pop[idx_arr[j]].crowd_dist /= num_obj;
			}
		}
		return;
	}// end function assign_crowding_distance

	bool mop_common::assign_crowd_dist_harmonic(population &pop,const idx_array &idx_arr,int assi_size,int neighbor_num)
	{
		if ( pop.size()<static_cast<size_t>(neighbor_num+1) ) return false;// not enough points for k neighbors
		d_mat dist(assi_size,d_array(assi_size));
		int i;
		for ( i=0;i<assi_size;i++ )
			pop[idx_arr[i]].crowd_dist=0.0;
		// calc mutual distance between any two initial member of elitist
		int j;
		for ( i=0;i<assi_size;i++ )
		{
			for ( j=0;j<assi_size;j++ )
			{
				if ( j>i )
				{
					dist[i][j]=vec_distance(pop[i].obj,pop[j].obj);
					dist[j][i]=dist[i][j];// symmetric
				}
				else if ( j==i )
					dist[i][j]=INF;// assign maximum value to exclude itself
			}
		}
		idx_mat rank(assi_size,idx_array(assi_size));
		for ( i=0;i<assi_size;i++ )
		{
			init_idx_array(rank[i],assi_size);
			quicksort_array(dist[i],rank[i],assi_size);// sorting accords to distance
		}
		// calc harmonic crowding distance of every archive member
		for ( i=0;i<assi_size;i++ )
		{
			pop[i].crowd_dist=0.0;
			int k;
			for ( k=0;k<neighbor_num;k++ )
			{
				int nei_idx=rank[i][k];
				pop[i].crowd_dist += 1.0/dist[i][nei_idx];
			}// for kth neighbor
			pop[i].crowd_dist=1.0/pop[i].crowd_dist;
		}// for every archive member
		return true;
	}// end function assign_crowd_dist_harmonic

	void mop_common::assign_crowd_dist_entropy(population &pop,const idx_array &idx_arr,int assi_size)
	{
		int i;
		for ( i=0;i<assi_size;i++ )
			pop[idx_arr[i]].crowd_dist=0.0;

		int num_obj=pop[0].obj.size();
		idx_mat rank(num_obj,idx_array(assi_size));// sorting order index of every objective
		// calc mutual distances
		int k;
		for ( k=0;k<num_obj;k++ )
		{
			rank[k]=idx_arr;
			quicksort_ind_score(pop,k,rank[k],assi_size);// sorting accords to kth objective,record precedence index in rank
		}
		int k_num=2;

		int obj_low;
		int obj_high;
		int low_idx,high_idx;

		// calc crowding distance
		const double log2=log(2.0);
		for ( k=0;k<num_obj;k++ )
		{
			obj_low=0;
			obj_high=assi_size-1;
			low_idx=rank[k][obj_low];
			high_idx=rank[k][obj_high];
			double obj_rng=pop[high_idx].obj[k]-pop[low_idx].obj[k];// f_max-f_min for normalization
			if ( 0!=obj_rng )
			{
				pop[low_idx].crowd_dist = INF;
				pop[high_idx].crowd_dist = INF;
			}
			else
			{
				pop[low_idx].crowd_dist += 0;
				pop[high_idx].crowd_dist += 0;
			}

			// calc crowding metric of current objective
			int cur_idx;
			for ( i=1;i<assi_size-1;i++ )
			{
				cur_idx=rank[k][i];
				if ( 0!=obj_rng )
				{
					// find k adjacent points index at current objective
					int left_idx,right_idx;
					left_idx=rank[k][i-1];
					right_idx=rank[k][i+1];
					// calc crowding entropy at current objective
					double dl,du,c;
					dl=abs(pop[cur_idx].obj[k]-pop[left_idx].obj[k]);
					du=abs(pop[right_idx].obj[k]-pop[cur_idx].obj[k]);
					c=dl+du;
					double pl;
					double pu;
					if ( 0!=c )
					{
						pl=dl/c;
						pu=du/c;
					}
					else
						pl=pu=0;
					double item_1,item_2;
					if ( pl!=0 )
						item_1=dl*(log(pl)/log2);
					else
						item_1=0;
					if ( pu!=0 )
						item_2=du*(log(pu)/log2);
					else
						item_2=0;
					if ( INF!=pop[cur_idx].crowd_dist )
						pop[cur_idx].crowd_dist += -1.0*( (item_1+item_2)/obj_rng );
				}// if ( 0!=obj_rng )
				else
					pop[cur_idx].crowd_dist += 0;
			}// for every specified non-extreme point
		}// for every objective
	}// end function assign_crowd_dist_entropy

	void mop_common::calc_crowd_dist(const elite_archive &ext_archive,individual &ind,int neighbor_num)// neighbor_num argument IGNORED
	{
		ind.crowd_dist=0.0;

		double dl,du;
		double min_obj,max_obj;
		double cur_obj;
		const double log2=log(2.0);
		int num_obj=ext_archive.front().obj.size();
		elite_archive::const_iterator cur_itr;
		int k;
		for ( k=0;k<num_obj;k++ )
		{
			// find two extreme points index at current objective
			dl=INF;
			du=INF;
			min_obj=INF;
			max_obj=MIN;
			double cur_d;
			for ( cur_itr=ext_archive.begin();cur_itr!=ext_archive.end();++cur_itr )
			{
				cur_obj=cur_itr->obj[k];
				if ( cur_obj<min_obj )
					min_obj=cur_obj;
				if ( cur_obj>max_obj )
					max_obj=cur_obj;
				cur_d=ind.obj[k]-cur_obj;
				if ( cur_d>0 )
				{
					if ( cur_d<dl )
						dl=cur_d;
				}
				else
				{
					cur_d=-1.0*cur_d;
					if ( cur_d<du )
						du=cur_d;
				}
			}// for every archive member
			if ( ind.obj[k]==min_obj || ind.obj[k]==max_obj )
			{
				ind.crowd_dist=INF;
				return;
			}
			double obj_rng=max_obj-min_obj;// f_max-f_min for normalization
			if ( 0!=obj_rng )
			{
				// calc crowding entropy at current objective
				ind.crowd_dist += (dl+du)/obj_rng;
			}// if ( 0!=obj_rng )
			else
				ind.crowd_dist += 0.0;
		}// for every objective
	}// end function calc_crowd_dist

	void mop_common::calc_crowd_entropy(const elite_archive &ext_archive,individual &ind,int neighbor_num)
	{
		ind.crowd_dist=0.0;

		double dl,du;
		double c;
		double min_obj,max_obj;
		double cur_obj;
		const double log2=log(2.0);
		int num_obj=ext_archive.front().obj.size();
		elite_archive::const_iterator cur_itr;
		int k;
		for ( k=0;k<num_obj;k++ )
		{
			// find two extreme points index at current objective
			dl=INF;
			du=INF;
			min_obj=INF;
			max_obj=MIN;
			double cur_d;
			for ( cur_itr=ext_archive.begin();cur_itr!=ext_archive.end();++cur_itr )
			{
				cur_obj=cur_itr->obj[k];
				if ( cur_obj<min_obj )
					min_obj=cur_obj;
				if ( cur_obj>max_obj )
					max_obj=cur_obj;
				cur_d=ind.obj[k]-cur_obj;
				if ( cur_d>0 )
				{
					if ( cur_d<dl )
						dl=cur_d;
				}
				else
				{
					cur_d=-1.0*cur_d;
					if ( cur_d<du )
						du=cur_d;
				}
			}// for every archive member
			if ( ind.obj[k]==min_obj || ind.obj[k]==max_obj )
			{
				ind.crowd_dist=INF;
				return;
			}
			double obj_rng=max_obj-min_obj;// f_max-f_min for normalization
			if ( 0!=obj_rng )
			{
				c=dl+du;
				double pl;
				double pu;
				if ( 0!=c )
				{
					pl=dl/c;
					pu=du/c;
				}
				else
					pl=pu=0;
				double item_1,item_2;
				if ( pl!=0 )
					item_1=dl*(log(pl)/log2);
				else
					item_1=0;
				if ( pu!=0 )
					item_2=du*(log(pu)/log2);
				else
					item_2=0;

				ind.crowd_dist += -1.0*( (item_1+item_2)/obj_rng );
			}// if ( 0!=obj_rng )
			else
				ind.crowd_dist += 0;
		}// for every objective
	}// end function calc_crowd_entropy

	void mop_common::calc_crowd_harmonic(const elite_archive &ext_archive,individual &ind,int neighbor_num)
	{
		ind.crowd_dist=0.0;
		int arc_size=ext_archive.size();

		int num_obj=ext_archive.front().obj.size();
		d_array dist(arc_size,0.0);
		// calc mutual distances
		elite_archive::const_iterator cur_itr;
		int i;
		for ( cur_itr=ext_archive.begin(),i=0;cur_itr!=ext_archive.end();++cur_itr,i++ )
			dist[i]=vec_distance(ind.obj,cur_itr->obj);
		for ( i=0;i<neighbor_num;i++ )
		{
			d_array::iterator nth_itr=dist.begin();
			advance(nth_itr,i);
			nth_element(dist.begin(),nth_itr,dist.end());
			ind.crowd_dist+= 1.0/(*nth_itr);
		}
		ind.crowd_dist=1.0/ind.crowd_dist;
	}// end function calc_crowd_harmonic

	void mop_common::trunc_external(elite_archive &ext_archive,int trunc_type,int max_size,int &trunc_count,int neighbor_num)
	{
		switch (trunc_type)
		{
		case crowd_entropy:
			trunc_archive_entropy(ext_archive,max_size,trunc_count);
			break;
		case crowd_harmonic:
			trunc_archive_harmonic(ext_archive,max_size,trunc_count,neighbor_num);
			break;
		default:
			trunc_archive_cd(ext_archive,max_size,trunc_count);
			break;
		}
	}// end function trunc_external

	void mop_common::trunc_archive_cd(elite_archive &ext_archive,int max_size,int &trunc_count,int neighbor_num)
	{
		int arc_size=ext_archive.size();
		if ( arc_size>max_size )
		{
			population tmp_arc;
			copy(ext_archive.begin(),ext_archive.end(),back_inserter<population>(tmp_arc));

			idx_array idx_arr(arc_size);
			init_idx_array(idx_arr,arc_size);
			crowding_fill_indices(tmp_arc,tmp_arc,max_size,0,arc_size,idx_arr,crowd_dist);
			// copy back
			ext_archive.clear();
			int i;
			for ( i=0;i<max_size;i++ )
				ext_archive.push_back(tmp_arc[i]);

			trunc_count++;
		}
	}// end function trunc_archive_harmonic

	void mop_common::trunc_archive_harmonic(elite_archive &ext_archive,int max_size,int &trunc_count,int neighbor_num)
	{
		int arc_size=ext_archive.size();
		if ( arc_size>max_size )
		{
			population tmp_arc;
			copy(ext_archive.begin(),ext_archive.end(),back_inserter<population>(tmp_arc));

			idx_array idx_arr(arc_size);
			init_idx_array(idx_arr,arc_size);
			crowding_fill_indices(tmp_arc,tmp_arc,max_size,0,arc_size,idx_arr,crowd_harmonic);
			// copy back
			ext_archive.clear();
			int i;
			for ( i=0;i<max_size;i++ )
				ext_archive.push_back(tmp_arc[i]);

			trunc_count++;
		}
	}// end function trunc_archive_harmonic

	void mop_common::trunc_archive_entropy(elite_archive &ext_archive,int max_size,int &trunc_count)
	{
		int arc_size=ext_archive.size();
		if ( arc_size>max_size )
		{
			population tmp_arc;
			copy(ext_archive.begin(),ext_archive.end(),back_inserter<population>(tmp_arc));

			idx_array idx_arr(arc_size);
			init_idx_array(idx_arr,arc_size);
			crowding_fill_indices(tmp_arc,tmp_arc,max_size,0,arc_size,idx_arr,crowd_entropy);
			// copy back
			ext_archive.clear();
			int i;
			for ( i=0;i<max_size;i++ )
				ext_archive.push_back(tmp_arc[i]);

			trunc_count++;
		}
	}// end function trunc_archive_entropy

	void mop_common::trunc_archive_harmonic_exact(elite_archive &ext_archive,int max_size,int &trunc_count,int neighbor_num)
	{
		int arc_size=ext_archive.size();
		int trun_num=arc_size-max_size;

		idx_array valid_ind_idx(arc_size);// erased flag
		int i;

		valid_ind_idx.assign(arc_size,1);
		population tmp_arc;
		copy(ext_archive.begin(),ext_archive.end(),back_inserter<population>(tmp_arc));
		d_mat dist(arc_size,d_array(arc_size));
		// calc mutual distance between any two initial member of elitist
		int j;
		for ( i=0;i<arc_size;i++ )
		{
			for ( j=0;j<arc_size;j++ )
			{
				if ( j>i )
				{
					dist[i][j]=vec_distance(tmp_arc[i].obj,tmp_arc[j].obj);
					dist[j][i]=dist[i][j];// symmetric
				}
				else if ( j==i )
					dist[i][j]=INF;// assign maximum value to exclude itself
			}
		}
		if ( trun_num>0 )
		{
			int cur_arc_size=arc_size;
			int min_idx;
			double min_dist;

			idx_mat rank(arc_size,idx_array(arc_size));
			for ( i=0;i<arc_size;i++ )
			{
				init_idx_array(rank[i],arc_size);
				quicksort_array(dist[i],rank[i],arc_size);// sorting accords to distance
			}
			int valid_count;
			while ( cur_arc_size>max_size )
			{
				// calc harmonic crowding distance of every archive member
				for ( i=0;i<arc_size;i++ )
				{
					if ( invalid!=valid_ind_idx[i] )
					{
						tmp_arc[i].crowd_dist=0.0;
						int k;
						for ( k=1;k<=neighbor_num;k++ )
						{
							valid_count=0;
							int n;
							for ( n=0;n<arc_size;n++ )
							{
								int nei_idx=rank[i][n];
								if ( invalid!=valid_ind_idx[nei_idx] )// choose existing neighbor poinit only
								{
									valid_count++;
									if ( valid_count==k )
									{
										tmp_arc[i].crowd_dist+= 1.0/dist[i][nei_idx];
										break;
									}
								}
							}// for every neighboring point
						}// for kth neighbor
						tmp_arc[i].crowd_dist=1.0/tmp_arc[i].crowd_dist;
					}
				}// for every archive member
				// find crowdest index
				min_dist=INF;
				for ( i=0;i<arc_size;i++ )
				{
					if ( invalid!=valid_ind_idx[i] )
						if ( tmp_arc[i].crowd_dist<min_dist )
						{
							min_dist=tmp_arc[i].crowd_dist;
							min_idx=i;
						}
				}// for every archive member
				valid_ind_idx[min_idx]=invalid;// label this point "erased"
				cur_arc_size--;
			}// while ( cur_arc_size>max_size )
			trunc_count++;
		}// trun_num>0
		// copy back
		ext_archive.clear();// clear current content
		idx_array rank(arc_size);
		init_idx_array(rank,arc_size);
		quicksort_ind_score(tmp_arc,first,rank,arc_size);
		for ( i=0;i<arc_size;i++ )
			if ( invalid!=valid_ind_idx[rank[i]] )
			{ 
				// ascending order for output convienience,sorting accords to first objective
				ext_archive.push_back(tmp_arc[rank[i]]);
			}
	}// end functon trunc_archive_harmonic_exact

	void mop_common::trunc_archive_entropy_exact(elite_archive &ext_archive,int max_size,int &trunc_count)
	{
		int arc_size=ext_archive.size();
		int trun_num=arc_size-max_size;
		population tmp_arc(arc_size);
		copy(ext_archive.begin(),ext_archive.end(),tmp_arc.begin());

		int i;
		idx_array valid_ind_idx(arc_size);// erased flag
		valid_ind_idx.assign(arc_size,1);
		ext_archive.clear();
		if ( trun_num>0 )
		{
			for ( i=0;i<arc_size;i++ )
				tmp_arc[i].crowd_dist=0.0;

			int num_obj=2;
			idx_mat rank(num_obj,idx_array(arc_size));// sorting order index of every objective
			// calc mutual distances
			int k;
			for ( k=0;k<num_obj;k++ )
			{
				init_idx_array(rank[k],arc_size);
				quicksort_ind_score(tmp_arc,k,rank[k],arc_size);// sorting accords to kth objective,record precedence index in rank
			}
			int k_num=2;
			int neighbor_num=k_num/2;

			int cur_arc_size=arc_size;

			int obj_low;
			int obj_high;
			int low_idx,high_idx;
			while ( cur_arc_size>max_size )
			{
				// calc crowding distance
				const double log2=log(2.0);
				int i;
				int k;
				for ( k=0;k<num_obj;k++ )
				{
					obj_low=0;
					obj_high=arc_size-1;
					// find two extreme points index at current objective
					while ( invalid==valid_ind_idx[rank[k][obj_low]] )
						obj_low++;
					while ( invalid==valid_ind_idx[rank[k][obj_high]] )
						obj_high--;
					low_idx=rank[k][obj_low];
					high_idx=rank[k][obj_high];
					double obj_rng=abs(tmp_arc[high_idx].obj[k]-tmp_arc[low_idx].obj[k]);// f_max-f_min for normalization
					tmp_arc[low_idx].crowd_dist = INF;
					tmp_arc[high_idx].crowd_dist = INF;
					// calc crowding metric of current objective
					int cur_obj;
					cur_obj=obj_low+1;
					int cur_idx;
					for ( i=0;i<arc_size && cur_obj<obj_high;i++,cur_obj++ )
					{
						cur_idx=rank[k][cur_obj];
						if ( invalid!=valid_ind_idx[cur_idx] && cur_idx!=low_idx && cur_idx!=high_idx )
						{
							// find k adjacent points index at current objective
							int obj_left,obj_right;
							obj_left=cur_obj-1;
							obj_right=cur_obj+1;
							while ( invalid==valid_ind_idx[rank[k][obj_left]] )
								obj_left++;
							while ( invalid==valid_ind_idx[rank[k][obj_right]] )
								obj_right--;
							int left_idx,right_idx;
							left_idx=rank[k][obj_left];
							right_idx=rank[k][obj_right];
							// calc crowding entropy at current objective
							double dl,du,c;
							dl=abs(tmp_arc[cur_idx].obj[k]-tmp_arc[left_idx].obj[k]);
							du=abs(tmp_arc[right_idx].obj[k]-tmp_arc[cur_idx].obj[k]);
							c=dl+du;
							double pl;
							double pu;
							if ( 0!=c )
							{
								pl=dl/c;
								pu=du/c;
							}
							else
								pl=pu=0;
							double item_1,item_2;
							if ( pl!=0 )
								item_1=dl*(log(pl)/log2);
							else
								item_1=0;
							if ( pu!=0 )
								item_2=du*(log(pu)/log2);
							else
								item_2=0;
							if ( INF!=tmp_arc[cur_idx].crowd_dist )
								tmp_arc[cur_idx].crowd_dist += -1.0*( (item_1+item_2)/obj_rng );
						}
					}// for every non-extreme archive member
				}// for every objective
				// find crowdest point and erase
				double min_dist=INF;
				int min_idx;
				for ( i=0;i<arc_size;i++ )
				{
					if ( invalid!=valid_ind_idx[i] )
					{
						if ( tmp_arc[i].crowd_dist<min_dist )
						{
							min_dist=tmp_arc[i].crowd_dist;
							min_idx=i;
						}
					}
				}
				valid_ind_idx[min_idx]=-1;// label this point "erased"
				cur_arc_size--;
			}// while current archive size > maximum size
			trunc_count++;
			// copy back
			for ( i=0;i<arc_size;i++ )
				if ( invalid!=valid_ind_idx[rank[first][i]] )
				{ 
					// ascending order for output convienience,sorting accords to first objective
					ext_archive.push_back(tmp_arc[rank[first][i]]);
				}
		}// if ( trun_num>0 )
		else
		{
			// copy back
			idx_array rank(arc_size);
			init_idx_array(rank,arc_size);
			quicksort_ind_score(tmp_arc,first,rank,arc_size);
			for ( i=0;i<arc_size;i++ )
				if ( invalid!=valid_ind_idx[rank[i]] )
				{ 
					// ascending order for output convienience,sorting accords to first objective
					ext_archive.push_back(tmp_arc[rank[i]]);
				}
		}
	}// end functon trunc_archive_entropy_exact

	void mop_common::fill_nondominated_sort(population &mixed_pop, population &new_pop,int pop_size,int trunc_type)
	{
		int flag;
		int j;
		int end;
		int front_size;
		int archieve_size;
		int rank=1;
		int mixed_size=mixed_pop.size();
		int cur_pos;
		front_list f_idx;
		front_list elite;

		front_list::iterator cur_pop;
		front_list::iterator cur_arc;

		front_size = 0;
		archieve_size=0;

		int i;
		for ( i=0; i<mixed_size; i++ )
			f_idx.push_back(i);

		cur_pos=0;
		do
		{
			elite.push_back(f_idx.front());
			front_size = 1;
			cur_pop = f_idx.erase(f_idx.begin());

			do
			{
				cur_arc = elite.begin();
				if (cur_pop==f_idx.end())
					break;
				do
				{
					end = 0;
					flag = check_dominance (mixed_pop[*cur_pop],mixed_pop[*cur_arc]);
					if (flag == better)
					{
						f_idx.push_front(*cur_arc);
						cur_arc=elite.erase(cur_arc);
						front_size--;
					}
					if (flag == indiff)
						++cur_arc;
					if (flag == worse)
						end = 1;
				}
				while ( end!=1 && cur_arc!=elite.end() );
				if ( flag!=worse )
				{
					elite.push_back(*cur_pop);
					front_size++;
					cur_pop = f_idx.erase(cur_pop);
				}
				else
					++cur_pop;
			} while ( cur_pop != f_idx.end() );
			cur_arc = elite.begin();
			if ( (archieve_size+front_size) <= pop_size )
			{
				do
				{
					new_pop[cur_pos]=mixed_pop[*cur_arc];
					new_pop[cur_pos].rank = rank;
					archieve_size+=1;
					++cur_arc;
					cur_pos++;
				}
				while ( cur_arc != elite.end() );
				rank+=1;
			}
			else
			{
				crowding_fill_list(mixed_pop,new_pop,pop_size,cur_pos,front_size,elite,trunc_type);
				archieve_size = pop_size;
				for (j=cur_pos; j<pop_size; j++)
					new_pop[j].rank = rank;
			}
			elite.clear();
		} while (archieve_size < pop_size);

		return;
	}// end function fill_nondominated_sort

	double mop_common::zdt3_calc_gamma(d_mat &opt_coor,int sample_size)
	{
		const double front_bnd[5][2]=
		{
			{0,0.0830015349},
			{0.1822287280,0.2577623634},
			{0.4093136748,0.4538821041},
			{0.6183967944,0.6525117038},
			{0.8233317983,0.8518328654}
		};

		int arc_size=opt_coor.size();
		d_array distance(sample_size,0.0);// archive member to frontier member distances
		d_array d(sample_size,0.0);// nearest distance of each archive member to frontier member:D_i
		vector<d_array> frontier(sample_size,d_array(2));// sampled pareto-optimal point
		int i,j;

		// sampling sample_size point from pareto-optimal front
		enum {low,high};
		const int f_seg_num=5;
		for ( i=0;i<f_seg_num;i++ )
		{
			uniform_real<> dist_uni_bnd(front_bnd[i][low],front_bnd[i][high]);
			variate_generator<mt19937&, uniform_real<> > rnd_front_coor(gen,dist_uni_bnd);
			for ( j=0;j<sample_size/f_seg_num;j++ )
			{
				double y1=rnd_front_coor();
				double y2=1-sqrt(y1)-y1*sin(10*M_PI*y1);
				double tmp[2]={y1,y2};
				frontier.push_back(d_array(tmp,tmp+2));
			}
		}
		double gamma=0.0;
		d_mat::iterator cur_itr(opt_coor.begin());
		for ( i=0;cur_itr!=opt_coor.end();++cur_itr,i++ )
		{
			int j;
			d[i]=INF;
			for ( j=0;j<sample_size;j++ )
			{
				// calc d_i
				distance[j]=vec_distance(*cur_itr,frontier[j]);
				if ( distance[j]<d[i] )
					d[i]=distance[j];
			}// for every sampling point at frontier
			gamma += d[i];
		}// for every archive member
		gamma/=arc_size;
		return gamma;
	}// end function zdt3_calc_gamma

	void mop_common::zdt3_assess(const d_mat &coor,int sample_size,const point& ref_point,perf_indice &p_ind)
	{
		const double front_bnd[5][2]=
		{
			{0,0.0830015349},
			{0.1822287280,0.2577623634},
			{0.4093136748,0.4538821041},
			{0.6183967944,0.6525117038},
			{0.8233317983,0.8518328654}
		};

		int arc_size=coor.size();
		d_array distance(sample_size,0.0);// archive member to frontier member distances
		d_array d(sample_size,0.0);// nearest distance of each archive member to frontier member:D_i
		vector<d_array > frontier(sample_size,d_array(2));// sampled pareto-optimal point
		int i,j;

		// sampling sample_size point from pareto-optimal front
		enum {low,high};
		const int f_seg_num=5;
		for ( i=0;i<f_seg_num;i++ )
		{
			uniform_real<> dist_uni_bnd(front_bnd[i][low],front_bnd[i][high]);
			variate_generator<mt19937&, uniform_real<> > rnd_front_coor(gen,dist_uni_bnd);
			for ( j=0;j<sample_size/f_seg_num;j++ )
			{
				double y1=rnd_front_coor();
				double y2=1-sqrt(y1)-y1*sin(10*M_PI*y1);
				double tmp[2]={y1,y2};
				frontier.push_back(d_array(tmp,tmp+2));
			}
		}
		p_ind.gamma=0.0;
		d_mat::const_iterator cur_itr(coor.begin());
		for ( i=0;cur_itr!=coor.end();++cur_itr,i++ )
		{
			int j;
			d[i]=INF;
			for ( j=0;j<sample_size;j++ )
			{
				// calc d_i
				distance[j]=vec_distance(*cur_itr,frontier[j]);
				if ( distance[j]<d[i] )
					d[i]=distance[j];
			}// for every sampling point at frontier
			p_ind.gamma += d[i];
		}// for every archive member
		p_ind.gamma/=arc_size;
		// push elitist to temp array for calculation convenience
		int num_obj=2;
		idx_array rank(arc_size);// sorting order index of the first objective
		double dl,du;// distances between extreme points at pareto-optimal frontier and boundary points at obtained frontier
		int fl_idx,fu_idx;
		// find extreme points and boundary points,respectively
		// calc d_l,d_u index
		fl_idx=0;
		fu_idx=sample_size-1;
		for ( i=0;i<sample_size;i++ )
		{
			if ( i<5 )
				if ( frontier[i][y1]<frontier[fl_idx][y1] )
					fl_idx=i;
			if ( i>=sample_size-5 )
				if ( frontier[i][y1]>frontier[fu_idx][y1] )
					fu_idx=i;
		}
		init_idx_array(rank,arc_size);
		quicksort_dmat(coor,rank,arc_size,first_obj_le());// sorting accords to kth objective,record precedence index in rank

		// calc spread metric:delta
		int dist_num=arc_size-1;
		int first_idx=0;
		int last_idx=dist_num;
		d_array m_d(dist_num,0.0);

		dl=vec_distance(coor[rank[first_idx]],frontier[fl_idx]);
		du=vec_distance(coor[rank[last_idx]],frontier[fu_idx]);
		double d_mean=0.0;
		for ( i=0;i<dist_num;i++ )
		{
			m_d[i]=vec_distance(coor[rank[i+1]],coor[rank[i]]);
			d_mean+=m_d[i];
		}
		d_mean/=dist_num;
		double sum_d=0.0;
		for ( i=0;i<dist_num;i++ )
			sum_d+=abs(m_d[i]-d_mean);
		p_ind.delta=(dl+du+sum_d)/(dl+du+(arc_size-1)*d_mean);
		// calc hv
		p_ind.hv=(ref_point.y1-coor[rank[last_idx]][y1])*(ref_point.y2-coor[rank[last_idx]][y2]);
		for ( i=1;i<arc_size;i++ )
			p_ind.hv += (coor[rank[i]][y1]-coor[rank[i-1]][y1])*(ref_point.y2-coor[rank[i-1]][y2]);
	}//end function zdt3_assess
}// end namespace mop
