#include <vector>
#include <stdexcept> // range_error
#include <ostream>
#include <sstream>
//#include <boost/array.hpp>
//#include <boost/type.hpp>

#include <boost/random.hpp>
#include "algo_datastruct.h"

using std::vector;
using std::stringstream;
using std::logic_error;
using std::string;
//using boost::array;
//using boost::extents;
using boost::normal_distribution;
using boost::mt19937;
using boost::variate_generator;

extern boost::mt19937 gen;

void allocate_ind(individual &ind,int NumDims,int stra_num,int num_obj)
{
	ind.x.resize(NumDims);
	ind.obj.resize(num_obj);
	if ( stra_num )
	{
		ind.stra.resize(stra_num);
		int i;
		for ( i=0;i<stra_num;i++ )
			ind.stra[i].resize(1);
	}
	ind.rank=-1;
	ind.constr_viol=0.0;
	ind.crowd_dist=0.0;
}

void allocate_pop(population &pop,size_t NumDims,int stra_num,int num_obj)
{
	size_t pop_size=pop.size();
	size_t i;
	for ( i=0;i<pop_size;i++ )
		allocate_ind(pop[i],NumDims,stra_num,num_obj);
}

void reallocate_stra(population &pop,int stra_idx,int size,double val)
{
	int pop_size=pop.size();
	int i;
	for ( i=0;i<pop_size;i++ )
		pop[i].stra[stra_idx].assign(size,val);
}

void reallocate_stra(population &pop,population &trial_pop,int stra_idx,int size,double val)
{
	int pop_size=pop.size();
	int tri_pop_size=trial_pop.size();
	if ( pop_size!=tri_pop_size )
	{
		stringstream ss_err;
		ss_err<<"\n"
			<<"Error in reallocate_stra:"
			<<"pop and trial_pop's size mismatch."
			<<"\n";
		throw logic_error(ss_err.str());
	}
	int i;
	for ( i=0;i<pop_size;i++ )
	{
		pop[i].stra[stra_idx].assign(size,val);
		trial_pop[i].stra[stra_idx].assign(size,val);
	}
}

void alg_stat::alloc_space(int max_run,int pop_size,int num_dims)
{
	// allocate centroid
	cen_ind.resize(num_dims,0.0);
	// allocate pbest
	pbest.resize(num_dims);
	allocate_pop(pbest,num_dims);
	// allocate gbest of single run
	allocate_ind(gbest,num_dims);
	// allocate gbest of all runs
	allocate_ind(bst_ind,num_dims);
	gen_bst_val.resize(max_run);
	gen_div_val.resize(max_run);
	gen_rad_val.resize(max_run);
	// initialize stat vector 
	int i;
	radius.assign(num_dims,0.0);
	delta_x.resize(pop_size);
	for ( i=0;i<pop_size;i++ )
		delta_x[i].resize(num_dims,0.0);
}

bool individual::operator< (const individual &rhs) { return obj[0] < rhs.obj[0]; }
bool individual::operator<= (const individual &rhs) { return obj[0] <= rhs.obj[0]; }
bool individual::operator== (const individual &rhs) { return obj[0] == rhs.obj[0]; }
bool operator< (const individual &lhs,const individual &rhs) { return lhs.obj[0] < rhs.obj[0]; }
bool operator<= (const individual &lhs,const individual &rhs) { return lhs.obj[0] <= rhs.obj[0]; }

bool individual::operator> (const individual &rhs) { return obj[0] > rhs.obj[0]; }
bool individual::operator>= (const individual &rhs) { return obj[0] >= rhs.obj[0]; }
bool operator> (const individual &lhs,const individual &rhs) { return lhs.obj[0] > rhs.obj[0]; }
bool operator>= (const individual &lhs,const individual &rhs) { return lhs.obj[0] >= rhs.obj[0]; }
bool operator== (const individual &lhs,const individual &rhs) { return lhs.obj[0] == rhs.obj[0]; }

std::ostream& operator<<(std::ostream &os,const individual& ind){ return os<<ind.obj[0]; }

bool gen_avgbst_state::operator<(const gen_avgbst_state& rhs) { return cur_gen<rhs.cur_gen; }
bool operator<(const gen_avgbst_state& lhs,const gen_avgbst_state& rhs) { return lhs.cur_gen<rhs.cur_gen; }


bool is_dataLine(string line)
{
	if ( ""==line || '#'==line[0]  ) return false; // use lazy evaluation otherwise will cause subscript out of range error if ( ""==line )
	int i,len,data_count;
	len=line.size();
	data_count=0;
	for ( i=0;i<len;i++ )
	{
		if ( line[i]!=' ' && line[i]!='\t' )
			data_count++;
	}
	return (data_count>0);
}

void init_idx_array(idx_array &idx_arr,int size,int begin_idx)
{
	int i;
	int cur_idx=begin_idx;
	for ( i=0;i<size;i++,cur_idx++ )
		idx_arr[i]=cur_idx;
}

void gen_bi_norm_num(const bi_norm_var &var,double &x_1,double &x_2)
{
	normal_distribution<> std_norm_dist_1;
	normal_distribution<> std_norm_dist_2;
	variate_generator<mt19937&, normal_distribution<> > std_norm_1(gen, std_norm_dist_1);
	variate_generator<mt19937&, normal_distribution<> > std_norm_2(gen, std_norm_dist_2);
	double norm_1=std_norm_1();
	x_1=var.mu_x+var.sig_x*norm_1;
	x_2=var.mu_y+var.sig_y*(norm_1*var.rho+std_norm_2()*sqrt(1-var.rho*var.rho));
}// end function bi_norm_dist

bool calc_bi_norm_var(const d_array &x_vec,const d_array &y_vec,
	bi_norm_var &var)
{
	int x_size=x_vec.size();
	int y_size=y_vec.size();
	if ( 0==x_size || 0==y_size ) return false;
	int i,j;
	var.mu_x=0;
	for ( i=0;i<x_size;i++ )
		var.mu_x += x_vec[i];
	var.mu_x /= x_size;

	var.sig_x=0;
	for ( i=0;i<x_size;i++ )
	{
		double tmp=x_vec[i]-var.mu_x;
		var.sig_x += tmp*tmp;
	}
	var.sig_x /= x_size;
	var.sig_x=sqrt(var.sig_x);

	var.mu_y=0;
	for ( i=0;i<y_size;i++ )
		var.mu_y += y_vec[i];
	var.mu_y /= y_size;

	var.sig_y=0;
	for ( i=0;i<y_size;i++ )
	{
		double tmp=y_vec[i]-var.mu_y;
		var.sig_y += tmp*tmp;
	}
	var.sig_y /= y_size;
	var.sig_y=sqrt(var.sig_y);

	int comb_size=x_size*y_size;
	double mu_xy;
	mu_xy=0;
	for ( i=0;i<x_size;i++ )
		for ( j=0;j<y_size;j++ )
			mu_xy += x_vec[i]*y_vec[j];
	mu_xy /= comb_size;
	var.rho=mu_xy-var.mu_x*var.mu_y;
	return true;
}// end function calc_bi_norm_var
