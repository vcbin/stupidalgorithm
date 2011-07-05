#ifndef STUPIDALGO_SOP_ALLOCATE_OBJECT_DYNAMICALLY
#define STUPIDALGO_SOP_ALLOCATE_OBJECT_DYNAMICALLY

#include "com_alg.h"

int allocate_alg(int alg_type,boost::shared_ptr<alg_base> &pbase_alg,std::string conf_path );
int allocate_all_alg(
						std::vector<boost::shared_ptr<alg_base> > &vec_ptr_alg,
						std::string conf_path
					);
int allocate_func(int func_type,boost::shared_ptr<func_base> &pfunc_base,int num_dim );
void allocate_para(int algo_type,boost::shared_ptr<para_base> &pPara_Base,std::string conf_path );

#endif
