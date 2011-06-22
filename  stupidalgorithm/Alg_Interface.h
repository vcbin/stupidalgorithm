#ifndef JERRYLIU_PSO_ALG_INTERFACE
#define JERRYLIU_PSO_ALG_INTERFACE

#include <boost/shared_ptr.hpp>
#include <boost/thread/mutex.hpp>

#include "FuncDef.h"
#include "para_base.h"
#include "Algo_DataStruct.h"

using namespace test_func;

class alg_base
{
public:
	virtual void load_parameter( const boost::shared_ptr<para_base> &pPara )=0;
	virtual void set_para_val_desc( std::string str_para_val )=0;
	virtual void set_eval_func( boost::shared_ptr<func_base> &pFunc_Eval )=0;
	virtual void initialize()=0;
	virtual int run()=0;
	virtual void set_com_file_path(const common_path &com_out_path)=0;
	virtual void print_algo_stat( std::ostream &os_dst,int func_type,int algo_type,bool has_max_gen=NULL )=0;
	virtual void print_para_stat( std::ostream &os,bool has_max_gen,boost::shared_ptr<boost::mutex> pmut )=0;
};// end class alg_base

#endif

