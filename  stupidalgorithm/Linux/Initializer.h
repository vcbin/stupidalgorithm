#ifndef STUPIDALGO_INITIALIZER
#define STUPIDALGO_INITIALIZER

#include "Rand_Val.h"
#include "Algo_DataStruct.h"

class initializer
{
public:
	void ini_real_uni(population &pop,const val_range &bounds);// Real Variable:initialize by uniform distribution
	void ini_orthogonal(population &pop,int q_size,const val_range &bounds);
protected:
	void gen_ortho_tab(i_mat &b,int rows,int cols,int q_size);
	double gen_ortho_x(int j,double low_bnd,double high_bnd,int q_size);
};

#endif