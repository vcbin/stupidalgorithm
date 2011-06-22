#include "initializer.h"

extern boost::mt19937 gen;

using std::vector;

using boost::uniform_01;
using boost::mt19937;
using boost::variate_generator;

void initializer::ini_real_uni(population &pop,const val_range &bounds)
{
	size_t pop_size=pop.size();
	size_t num_dims=pop[0].x.size();
	size_t i,j;

	for ( i=0;i<pop_size;i++ )
	{
		for ( j=0;j<num_dims;j++ )
		{
			double low_bnd,high_bnd;

			low_bnd=bounds[j].min_val;
			high_bnd=bounds[j].max_val;

			uniform_01<> dist;
			variate_generator<mt19937&, uniform_01<> > rnd_num(gen, dist);

			pop[i].x[j]=low_bnd+rnd_num()*(high_bnd-low_bnd);
		}
	}
}// end function ini_real_uni

double initializer::gen_ortho_x(int q_idx,double low_bnd,double high_bnd,int q_size)
{
	if ( 0==q_idx )
		return low_bnd;
	else if ( q_size-1==q_idx )
		return high_bnd;
	else
		return ( low_bnd+q_idx*(high_bnd-low_bnd)/(q_size-1) );
}// end function gen_ortho_x

// generate orthogonal tab
void initializer::gen_ortho_tab(i_mat &b,int rows,int cols,int q_size)
{
	int i,j;
	for ( i=0;i<rows;i++ )
	{
		b[i][0]=(i/q_size)%q_size;
		b[i][1]=i%q_size;
	}
	for ( j=2;j<cols;j++ )
		for ( i=0;i<rows;i++ )
			b[i][j]=(b[i][0]*(j-1)+b[i][1])%q_size;
}

void initializer::ini_orthogonal(population &pop,int q_size,const val_range &bounds)
{
	size_t pop_size=pop.size();
	size_t num_dims=pop[0].x.size();
	i_mat o_tab(pop_size,i_array(num_dims));// orthogonal table

	gen_ortho_tab(o_tab,pop_size,num_dims,q_size);
	size_t i,j;
	// init real x value with table
	for ( i=0;i<pop_size;i++ )
	{
		for ( j=0;j<num_dims;j++ )
		{
			double low_bnd=bounds[j].min_val;
			double high_bnd=bounds[j].max_val;
			pop[i].x[j]=gen_ortho_x(o_tab[i][j],low_bnd,high_bnd,q_size);
		}
	}
}// end function ini_orthogonal