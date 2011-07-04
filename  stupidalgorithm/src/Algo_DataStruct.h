#ifndef STUPIDALGO_COMMON_DATASTRUCT
#define STUPIDALGO_COMMON_DATASTRUCT

#include <vector>
#include <ostream>
#include <cmath>
#include <limits>
#include <stdexcept>

typedef std::vector<double> d_array;
typedef std::vector<d_array> d_mat;
typedef std::vector<int> idx_array;
typedef std::vector<idx_array> idx_mat;
typedef idx_array i_array;
typedef idx_mat i_mat;

struct min_max
{
	double min_val;
	double max_val;
};
typedef std::vector<min_max> val_range;

struct individual
{
	d_array x;// independent variables
	// std::vector<int> gene;
	// array_type gene; // Gene of every individual
	// array_type bound;// lower/upper bound of each sub component of individual
	d_array obj;// function object value(s)
	d_mat stra;// strategy parameter
	int rank;// pareto-rank
	double crowd_dist;// crowding distance
	double constr_viol;// constraints violation indicator
	bool operator<(const individual &rhs);
	bool operator<=(const individual &rhs);
	bool operator>(const individual &rhs);
	bool operator>=(const individual &rhs);
	bool operator==(const individual &rhs);
	friend std::ostream& operator<<(std::ostream &os,const individual& ind);
};
// for fitness ordering
// const lhs version
bool operator< (const individual &lhs,const individual &rhs);
bool operator<= (const individual &lhs,const individual &rhs);
bool operator>(const individual &lhs,const individual &rhs);
bool operator>= (const individual &lhs,const individual &rhs);
bool operator== (const individual &lhs,const individual &rhs);

typedef std::vector<individual> population;

struct alg_stat
{
	alg_stat():
	stop_stag(false),
	avg(0.0),
	std(0.0),
	pos_diver(0.0),
	num_stag(0),
	eval_num(0),
	avg_radius(0.0),
	all_eval_num(0),
	run_avg_bst(0.0),
	run_bst_std(0.0),
	conv_count(0),
	conv_ratio(0.0),
	run_avg_time(0.0),
	run_avg_gen(0.0),
	all_ob_num(0),
	run_avg_stag_gen(0.0),
	avg_conv_eval_num(0.0) {}
void initialize(int max_run) // reset previous run stat variables to initial value
{
	conv_count=0;
	eval_num=0;
	num_stag=0;
	all_ob_num=0;
	all_eval_num=0;

	run_gbest_val.clear();
	run_gen_val.clear();

	int i;
	for ( i=0;i<max_run;i++ )
	{
		gen_bst_val[i].clear();
		gen_div_val[i].clear();
		gen_rad_val[i].clear();
	}
}
void alloc_space(int max_run,int pop_size,int num_dims);
void reset_run_stat() {num_stag=0;eval_num=0;}

// variables for each run
bool stop_stag;// whether config choose to stop on stagnation
double avg;// average fitness of population
double std;// population standard deviation
double pos_diver;// index of population position diversity
d_array cen_ind;// swarm centroid
int num_stag; // count of stagnation generations
int eval_num; // current evaluation number
double avg_radius;// mean search radius of population
population pbest;
individual gbest;

// variables for all runs/trials
int all_eval_num;// total function evaluation number
individual bst_ind;// best gbest solution of all runs
individual wst_ind;// worst gbest solution of all runs
double run_avg_bst;// average best fitness of all runs
double run_bst_std;// standard deviation of all gbest fitness
int conv_count;// convergence count
double conv_ratio;// convergence ratio
double run_avg_time;// average elapsed time of algorithm in second
d_array run_gbest_val;// gbest value of all runs
std::vector<int> run_gen_val;// aux variable
double run_avg_gen;
int all_ob_num;// total out of bounds number
double run_avg_stag_gen;// average stagnation generations ON STOP
double avg_conv_eval_num;// average function evaluation till convergence
// vector stat variables
d_mat gen_bst_val;// aux var,generational best_so_far of all runs
d_array gen_avg_bst_val;// average best_so_far of all runs in every generation
d_mat gen_div_val;// aux var,generational diversity of all runs
d_array gen_avg_div_val;// average diversity of all runs in every generation
d_mat gen_rad_val;// aux var,generational search radius of all runs
d_array gen_avg_rad_val;// average search radius of all runs in every generation
d_mat delta_x;// search step size,stat variable
d_array radius;// search radius of every individual,stat variable
};

// vector subtraction
template <typename T>
std::vector<T> operator-(const std::vector<T> &lhs,const std::vector<T> &rhs) 
{
	size_t num_dims=lhs.size();
	if ( num_dims != rhs.size() )
		throw std::range_error("two dimensons did not match!");
	size_t i;
	std::vector<T> res(num_dims);
	for (i=0;i<num_dims;i++)
	{
		res[i]=lhs[i]-rhs[i];
	}
	return res;
}

// vector addition
template <typename T>
std::vector<T> operator+(const std::vector<T> &lhs,const std::vector<T> &rhs) 
{
	size_t num_dims=lhs.size();
	if ( num_dims != rhs.size() )
		throw std::range_error( "two dimensons did not match!");
	size_t i;
	std::vector<T> res(num_dims);
	for (i=0;i<num_dims;i++)
	{
		res[i]=lhs[i]+rhs[i];
	}
	return res;
}

template <typename T>
std::vector<T>& operator+=(std::vector<T> &lhs,const std::vector<T> &rhs) {return lhs=operator+(lhs,rhs);}

// vector divided by float value
template <typename T>
std::vector<T> operator/(const std::vector<T> &lhs,double n) 
{
	size_t num_dims=lhs.size();
	if ( 0==n )
		throw std::logic_error("vector divided by ZERO!");
	size_t i;
	std::vector<T> res_vec(num_dims);
	for ( i=0;i<num_dims;i++ )
		res_vec[i]=lhs[i]/n;
	return res_vec;
}

template <typename T>
std::vector<T>& operator/=(std::vector<T> &lhs,double n) {return lhs=operator/(lhs,n);}

// vector comparison
template <typename type>
bool operator==(const std::vector<type> &lhs,const std::vector<type> &rhs)
{
	size_t num_dims=lhs.size();
	if ( num_dims != rhs.size() )
		return false;
	size_t i;
	for (i=0;i<num_dims;i++)
	{
		if ( lhs[i]!=rhs[i] )
			return false;
	}
	return true;
}

template <typename type>
bool operator!=(const std::vector<type> &lhs,const std::vector<type> &rhs){ return !operator==(lhs,rhs); }

template <typename T>
double vec_distance(const std::vector<T>& lhs,const std::vector<T>& rhs)
{
	size_t num_dims=lhs.size();
	if ( num_dims != rhs.size() )
		throw std::range_error("two dimensons did not match!");
	size_t i;
	double dis=0.0;
	double diff;
	for (i=0;i<num_dims;i++)
	{
		diff=lhs[i]-rhs[i];
		dis += diff*diff;
	}
	return sqrt(dis);
}

struct gen_avgbst_state
{
	int cur_gen;
	double gbest;
	bool operator<(const gen_avgbst_state& rhs);// for std::max_element
};
// const lhs version
bool operator<(const gen_avgbst_state& lhs,const gen_avgbst_state& rhs); 
// stop on stagnation state-record variable
typedef std::vector<gen_avgbst_state> Avgbst_states;


// Space Allocation methods
void allocate_ind(individual &ind,int NumDims,int stra_num=0,int num_obj=1);
void allocate_pop(population &pop,size_t NumDims,int stra_num=0,int num_obj=1);
void reallocate_stra(population &pop,int stra_idx,int size,double val=0.0);
void reallocate_stra(population &pop,population &trial_pop,int stra_idx,int size,double val=0.0);

struct common_path
{
	std::string all_bst_val_path;
	std::string avg_bst_path;
	std::string avg_div_path;
	std::string avg_rad_path;
	std::string stat_path;
};

inline bool is_int(double val) { return ( 0==(val-floor(val)) ); }
bool is_dataLine(std::string line);

const double INF=std::numeric_limits<double>::max();
#ifdef WIN32
const double MIN=std::numeric_limits<double>::lowest();
#else
const double MIN=-std::numeric_limits<double>::max();
#endif

void init_idx_array(idx_array &idx_arr,int size,int begin_idx=0);

struct bi_norm_var
{
	double mu_x;
	double mu_y;
	double sig_x;
	double sig_y;
	double rho;
	inline bi_norm_var& operator+=(const bi_norm_var& rhs) 
	{
		mu_x += rhs.mu_x;
		sig_x += rhs.sig_x;
		mu_y += rhs.mu_y;
		sig_y += rhs.sig_y;
		rho += rhs.rho;
		return *this;
	};
	inline bi_norm_var& operator/=(int n)
	{
		mu_x /= n;
		sig_x /= n;
		mu_y /= n;
		sig_y /= n;
		rho /= n;
		return *this;
	};
	friend std::ostream& operator<<(std::ostream &os,const bi_norm_var& rhs);
};
inline std::ostream& operator<<(std::ostream &os,const bi_norm_var& rhs)
{
	os<<rhs.mu_x<<" "
		<<rhs.sig_x<<" "
		<<rhs.mu_y<<" "
		<<rhs.sig_y<<" "
		<<rhs.rho;
	return os;
}
bool calc_bi_norm_var(const d_array &x_vec,const d_array &y_vec,bi_norm_var &var);
void gen_bi_norm_num(const bi_norm_var &var,double &x_1,double &x_2);

#endif
