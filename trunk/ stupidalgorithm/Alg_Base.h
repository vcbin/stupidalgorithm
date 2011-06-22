#ifndef JERRYLIU_PSO_ALG_INTERFACE
#define JERRYLIU_PSO_ALG_INTERFACE

#include <boost/shared_ptr.hpp>
#include <boost/thread/mutex.hpp>

#include "FuncDef.h"
#include "para_base.h"
#include "Algo_DataStruct.h"

using namespace benchmark;

class alg_base
{
public:
	virtual void load_parameter( const boost::shared_ptr<para_base> &ppara )=0;
	virtual void set_para_val_desc( std::string str_para_val )=0;
	virtual void set_eval_func( boost::shared_ptr<func_base> &pfunc_Eval )=0;
	virtual void initialize()=0;
	virtual int run()=0;
	virtual void set_com_file_path(const common_path &com_out_path)=0;
	virtual void print_algo_stat( std::ostream &os_dst,int func_type,int algo_type,bool max_gen_set=NULL )=0;
	virtual void print_para_stat( std::ostream &os,bool max_gen_set,boost::shared_ptr<boost::mutex> pmut )=0;
};// end class com_alg

#endif

