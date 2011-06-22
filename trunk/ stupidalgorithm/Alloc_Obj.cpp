#include "Alloc_Obj.h"
#include "All_Alg.h"
#include "All_Para.h"
#include "CmdLine_Para_Type.h"

using namespace pso;
using namespace pso::mpso;
using namespace pso::arpso;
using namespace pso::dpso;
using namespace pso::dpso_m;
using namespace pso::psobc;
using namespace de;
using namespace de::bbde;
using namespace de::sde;
using namespace de::spde;
using namespace de::jde;
using namespace ep;
using namespace ep::fep;
using namespace ep::ifep;
using namespace ga;
using namespace ga::dgea;
using namespace de::ede;
using namespace de::dmde;
using namespace de::de_eda;
using namespace de::nsde;
using namespace eda;

using namespace de::mosade;
using namespace de::mode;

using boost::shared_ptr;

using std::vector;
using std::string;
using std::logic_error;

int allocate_alg(int alg_type,shared_ptr<alg_base> &pbase_alg,string conf_path )
{
	switch (alg_type)// polymorphism of algorithm class
	{
	case alg_pso:
		{
			pbase_alg=shared_ptr<alg_base>(new pso_alg(conf_path));
			break;
		}
	case alg_mpso:
		{
			pbase_alg=shared_ptr<alg_base>(new mpso_alg(conf_path));
			break;
		}
	case alg_arpso:
		{
			pbase_alg=shared_ptr<alg_base>(new arpso_alg(conf_path));
			break;
		}
	case alg_dpso:
		{
			pbase_alg=shared_ptr<alg_base>(new dpso_alg(conf_path));
			break;
		}
	case alg_dpso_m:
		{
			pbase_alg=shared_ptr<alg_base>(new dpso_m_alg(conf_path));
			break;
		}
	case alg_psobc:
		{
			pbase_alg=shared_ptr<alg_base>(new psobc_alg(conf_path));
			break;
		}
	case alg_de:
		{
			pbase_alg=shared_ptr<alg_base>(new de_alg(conf_path));
			break;
		}
	case alg_bbde:
		{
			pbase_alg=shared_ptr<alg_base>(new bbde_alg(conf_path));
			break;
		}
	case alg_sde:
		{
			pbase_alg=shared_ptr<alg_base>(new sde_alg(conf_path));
			break;
		}
	case alg_my_sde:
		{
			pbase_alg=shared_ptr<alg_base>(new my_sde_alg(conf_path));
			break;
		}
	case alg_spde:
		{
			pbase_alg=shared_ptr<alg_base>(new spde_alg(conf_path));
			break;
		}
	case alg_jde:
		{
			pbase_alg=shared_ptr<alg_base>(new jde_alg(conf_path));
			break;
		}
	case alg_fep:
		{
			pbase_alg=shared_ptr<alg_base>(new fep_alg(conf_path));
			break;
		}
	case alg_ifep:
		{
			pbase_alg=shared_ptr<alg_base>(new ifep_alg(conf_path));
			break;
		}
	case alg_dgea:
		{
			pbase_alg=shared_ptr<alg_base>(new dgea_alg(conf_path));
			break;
		}
	case alg_ede:
		{
			pbase_alg=shared_ptr<alg_base>(new ede_alg(conf_path));
			break;
		}
	case alg_eda:
		{
			pbase_alg=shared_ptr<alg_base>(new eda_alg(conf_path));
			break;
		}
	case alg_deeda:
		{
			pbase_alg=shared_ptr<alg_base>(new de_eda_alg(conf_path));
			break;
		}
	case alg_dmde:
		{
			pbase_alg=shared_ptr<alg_base>(new dmde_alg(conf_path));
			break;
		}
	case alg_nsde:
		{
			pbase_alg=shared_ptr<alg_base>(new nsde_alg(conf_path));
			break;
		}
	case alg_mosade:
		{
			pbase_alg=shared_ptr<alg_base>(new mosade_alg(conf_path));
			break;
		}
	case alg_mode:
		{
			pbase_alg=shared_ptr<alg_base>(new mode_alg(conf_path));
			break;
		}
	default:
		return -1;
	}// end switch
	return 0;
}// end function allocate_alg


// OBSOLETE function
int allocate_all_alg( vector<shared_ptr<alg_base> > &vec_ptr_alg,
	string conf_path )
{
	int i;
	vec_ptr_alg.resize(max_algo_code);
	for ( i=1;i<=max_algo_code;i++ )
		allocate_alg(i,vec_ptr_alg[i-1],conf_path);
	return 0;
}

int allocate_func(int func_type,shared_ptr<func_base> &pfunc_base,int num_dim )
{
	switch (func_type)// polymorphism of test problem
	{
	case sphere:
		{
			pfunc_base=shared_ptr<func_base>(new func_Sphere(num_dim));
			break;
		}
	case griewank:
		{
			pfunc_base=shared_ptr<func_base>(new func_Griewank(num_dim));
			break;
		}
	case rastrigin:
		{
			pfunc_base=shared_ptr<func_base>(new func_Rastrigin(num_dim));
			break;
		}
	case ackey:
		{
			pfunc_base=shared_ptr<func_base>(new func_AckeyF1(num_dim));
			break;
		}
	case f5:
		{
			pfunc_base=shared_ptr<func_base>(new func_F5(num_dim));
			break;
		}
	case rosenbrock:
		{
			pfunc_base=shared_ptr<func_base>(new func_Rosenbrock(num_dim));
			break;
		}
	case step:
		{
			pfunc_base=shared_ptr<func_base>(new func_Step(num_dim));
			break;
		}
	case quartic_with_noise:
		{
			pfunc_base=shared_ptr<func_base>(new func_Quartic_with_Noise(num_dim));
			break;
		}
	case ws_location:
		{
			pfunc_base=shared_ptr<func_base>(new func_WS_Location(num_dim,"original_coor.txt"));
			break;
		}
	case f2:
		{
			pfunc_base=shared_ptr<func_base>(new func_F2(num_dim));
			break;
		}
	case f8:
		{
			pfunc_base=shared_ptr<func_base>(new func_F8(num_dim));
			break;
		}
	case camelback:
		{
			pfunc_base=shared_ptr<func_base>(new func_CamelBack(num_dim));
			break;
		}
	case f12:
		{
			pfunc_base=shared_ptr<func_base>(new func_F12(num_dim));
			break;
		}
	case f13:
		{
			pfunc_base=shared_ptr<func_base>(new func_F13(num_dim));
			break;
		}
	case f3:
		{
			pfunc_base=shared_ptr<func_base>(new func_F3(num_dim));
			break;
		}
	case f4:
		{
			pfunc_base=shared_ptr<func_base>(new func_F4(num_dim));
			break;
		}
	case f14:
		{
			pfunc_base=shared_ptr<func_base>(new func_F14(num_dim));
			break;
		}
	case f15:
		{
			pfunc_base=shared_ptr<func_base>(new func_F15(num_dim));
			break;
		}
	case foxhole:
		{
			pfunc_base=shared_ptr<func_base>(new func_Mod_Shekel(num_dim));
			break;
		}
	case langerman:
		{
			pfunc_base=shared_ptr<func_base>(new func_Mod_Langerman(num_dim));
			break;
		}
	case michaelwicz:
		{
			pfunc_base=shared_ptr<func_base>(new func_Mod_Michalewicz(num_dim));
			break;
		}
	case chebyshev:
		{
			pfunc_base=shared_ptr<func_base>(new func_Chebyshev(num_dim));
			break;
		}
	case zdt3:
		{
			pfunc_base=shared_ptr<func_base>(new func_ZDT3(num_dim));
			break;
		}
	default:
		return -1;
	}// end switch
	return 0;
}// end function allocate_func

void allocate_para(int algo_type,shared_ptr<para_base> &pPara_Base,string conf_path )
{
	if ( is_pso(algo_type) )
		pPara_Base=shared_ptr<para_base>(new pso_para(conf_path));
	else if ( is_de(algo_type) || is_hybrid_de(algo_type) || is_mode(algo_type) || is_nsde(algo_type) )
		pPara_Base=shared_ptr<para_base>(new de_para(conf_path));
	else if ( is_ep(algo_type) )
		pPara_Base=shared_ptr<para_base>(new ep_para(conf_path));
	else if ( is_ea(algo_type) )
		pPara_Base=shared_ptr<para_base>(new ga_para(conf_path));
	else if ( is_hybrid_de(algo_type) )
		pPara_Base=shared_ptr<para_base>(new de_para(conf_path));
	else if ( is_eda(algo_type) )
		pPara_Base=shared_ptr<para_base>(new eda_para(conf_path));
	else
	{
		string str_err("\nERROR:Illegal algo_type!");
		throw logic_error(str_err);
		return;
	}
}// end function allocate_para
