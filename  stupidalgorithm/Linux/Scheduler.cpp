#include <boost/program_options.hpp>
#include <iostream>
#include <fstream>
#include <boost/shared_ptr.hpp>
#include <boost/progress.hpp>
#include <exception>

#include <boost/thread.hpp>
#include <boost/threadpool.hpp> 
#include <boost/bind.hpp>

#include <ctime>

#include "All_Para.h"
#include "All_Alg.h"
#include "Batch_Conf.h"
#include "Alloc_Obj.h"
#include "Output_Base.h"
#include "Rand_Val.h"

#include "Scheduler.h"

using namespace pso;
using namespace pso::mpso;
using namespace de;
using namespace de::nsde;
using namespace ep;
using namespace ep::fep;

using namespace boost::threadpool;// NON-standard boost library currently

using boost::shared_ptr;
using boost::any_cast;
using namespace boost::program_options;

using std::string;
using std::map;
using std::vector;
using std::exception;
using std::cout;
using std::ifstream;
using std::ofstream;
using std::logic_error;
using std::bad_cast;
using std::transform;
using std::toupper;
using std::ostream;
using std::auto_ptr;

using std::stringstream;

using boost::lexical_cast;
using boost::progress_timer;
using boost::dynamic_pointer_cast;
using boost::mutex;
using boost::bind;
using boost::any;

boost::mt19937 gen;// global random number generator

const char* algo_name[]=
{
	"pso_std",
	"mpso",  
	"arpso",
	"dpso", 
	"my_dpso",
	"psobc",
	"de",
	"bbde",   
	"sde",
	"my_sde",
	"spde",
	"jde",
	"fep",  
	"ifep",
	"dgea",
	"ede",
	"eda",
	"de-eda",
	"dmde",
	"nsde",
	"mosade",
	"mode"
};

// mutex run_stat_file_mutex;// common stat file write mutex
mutex para_stat_file_mutex;// para stat file write mutex
mutex io_mutex;// cout mutex

// global common stat files for all algorithm
ofstream run_stat_file("Output//run_stat.txt");
ofstream para_stat_filelist("Output//para_stat_filelist.txt");
ofstream com_filelist("Output//common_filelist.txt");
ofstream para_stat_all_bst_val_filelist("Output//para_stat_all_bst_val_filelist.txt");

// extern const char* algo_name[];

void thread_fun( const shared_ptr<alg_base>& pAlg,int func_type,int algo_type,bool max_gen_set,
	bool has_seq_val,const shared_ptr<ostream>& prun_stat_line, 
	const shared_ptr<ostream>& ppara_stat_file,const shared_ptr<mutex>& pmut
	)
{
	try
	{
		pAlg->run();// actual computing part
		pAlg->print_algo_stat(*prun_stat_line,func_type,algo_type,max_gen_set);// record overall batch run stat
		{
			if ( has_seq_val )
				if ( ppara_stat_file )
					pAlg->print_para_stat(*ppara_stat_file,max_gen_set,pmut);
		}
	}
	catch (exception &e)
	{
		cout<<"\n"<<e.what()<<"\n";
		throw;
	}
}

// common stat filename generator
string scheduler::gen_stat_file_name(string str_common,string str_seq_desc,string str_append)
{
	string str_filename;
	if ( !str_seq_desc.empty() )
		str_filename="Output//"+str_common+"_"+str_seq_desc+"_"+str_append; 
	else
		str_filename="Output//"+str_common+"_"+str_append; 
	return str_filename;
}

string scheduler::gen_all_bst_val_file_name(string str_common,string str_seq_desc) 
{
	return gen_stat_file_name(str_common,str_seq_desc,"all_bst_val.txt"); 
}

string scheduler::gen_stat_file_name(string str_common,string str_seq_desc) 
{ 
	return gen_stat_file_name(str_common,str_seq_desc,"stat.txt"); 
}

string scheduler::gen_avg_bst_file_name(string str_common,string str_seq_desc) 
{
	return gen_stat_file_name(str_common,str_seq_desc,"gen_avg_bst.txt"); 
}

string scheduler::gen_avg_div_file_name(string str_common,string str_seq_desc) 
{ 
	return gen_stat_file_name(str_common,str_seq_desc,"gen_avg_div.txt"); 
}

string scheduler::gen_avg_rad_file_name(string str_common,string str_seq_desc) 
{ 
	return gen_stat_file_name(str_common,str_seq_desc,"gen_avg_rad.txt"); 
}

string scheduler::gen_para_stat_file_name(string str_common)
{
	return gen_stat_file_name(str_common,"","para_stat.txt");
}

// algorithm specific stat filename generator
string scheduler::gen_avg_vel_div_file_name(string str_common,string str_seq_desc)
{ return gen_stat_file_name(str_common,str_seq_desc,"gen_avg_vel_div.txt"); }

string scheduler::gen_avg_attr_file_name(string str_common,string str_seq_desc)
{ return gen_stat_file_name(str_common,str_seq_desc,"gen_avg_attr_perc.txt"); }

string scheduler::gen_avg_f_file_name(string str_common,string str_seq_desc)
{ return gen_stat_file_name(str_common,str_seq_desc,"gen_avg_f.txt"); }

string scheduler::gen_avg_pr_file_name(string str_common,string str_seq_desc)
{ return gen_stat_file_name(str_common,str_seq_desc,"gen_avg_pr.txt"); }

string scheduler::gen_avg_eta_file_name(string str_common,string str_seq_desc)
{ return gen_stat_file_name(str_common,str_seq_desc,"gen_avg_eta.txt"); }

string scheduler::gen_avg_bi_norm_stat_file_name(string str_common,string str_seq_desc)
{ return gen_stat_file_name(str_common,str_seq_desc,"gen_avg_bi_norm_stat.txt"); }

void scheduler::gen_common_file_path(common_path &com_out_path,string str_common,string str_seq_desc)
{
	com_out_path.all_bst_val_path=gen_all_bst_val_file_name(str_common,str_seq_desc);
	com_out_path.avg_bst_path=gen_avg_bst_file_name(str_common,str_seq_desc);
	com_out_path.avg_div_path=gen_avg_div_file_name(str_common,str_seq_desc);
	com_out_path.avg_rad_path=gen_avg_rad_file_name(str_common,str_seq_desc);
	com_out_path.stat_path=gen_stat_file_name(str_common,str_seq_desc);
}

void print_para_stat_All_Best_Val_Path( ostream &os,int func_type,string para_text,string str_path ) { os<<func_type<<"\t"<<para_text<<"\t"<<str_path<<"\n"; }

void scheduler::print_common_filelist(ostream &os,int func_type,int algo_type,const common_path &com_out_path)
{
	os<<func_type<<"\t"
		<<algo_type<<"\t"
		<<com_out_path.all_bst_val_path<<"\t"
		<<com_out_path.avg_bst_path<<"\t"
		<<com_out_path.avg_div_path<<"\t"
		<<com_out_path.avg_rad_path
		<<"\n";
}

void scheduler::algo_spec_output(int func_type,int algo_type,const shared_ptr<alg_base> &pAlg,string str_common,string str_seq_desc)
{
	if ( is_pso(algo_type) )
	{
		static ofstream gen_vel_div_filelist("Output//gen_vel_div_filelist.txt");
		shared_ptr<pso_alg> pPSO_alg=dynamic_pointer_cast<pso_alg>(pAlg);// downcast
		string filename=gen_avg_vel_div_file_name(str_common,str_seq_desc);
		pPSO_alg->set_gen_vel_div_file_name(filename);
		print_file_name(gen_vel_div_filelist,func_type,algo_type,filename);
		if ( is_mpso(algo_type) )
		{
			shared_ptr<mpso_alg> pMPSO_alg=dynamic_pointer_cast<mpso_alg>(pAlg);// downcast
			filename=gen_avg_attr_file_name(str_common,str_seq_desc);
			pMPSO_alg->set_gen_avg_attra_file_name(filename);
			ofstream gen_attr_perc_file("Output//mpso_gen_attr_filename.txt");
			print_file_name(gen_attr_perc_file,func_type,algo_type,filename);
		}
	}
	if ( is_de(algo_type) || is_hybrid_de(algo_type) || is_mode(algo_type) || is_nsde(algo_type) )
	{
		static ofstream gen_f_filelist("Output//gen_f_filelist.txt");
		static ofstream gen_pr_filelist("Output//gen_pr_filelist.txt");
		shared_ptr<de_alg> pDE_alg=dynamic_pointer_cast<de_alg>(pAlg);// downcast
		string filename=gen_avg_f_file_name(str_common,str_seq_desc);
		pDE_alg->set_gen_avg_f_file_name(filename);
		print_file_name(gen_f_filelist,func_type,algo_type,filename);
		filename=gen_avg_pr_file_name(str_common,str_seq_desc);
		pDE_alg->set_gen_avg_pr_file_name(filename);
		print_file_name(gen_pr_filelist,func_type,algo_type,filename);
		if ( is_nsde(algo_type) )
		{
			shared_ptr<nsde_alg> pnsde_alg=dynamic_pointer_cast<nsde_alg>(pAlg);// downcast
			filename=gen_avg_bi_norm_stat_file_name(str_common,str_seq_desc);
			pnsde_alg->set_gen_avg_bi_norm_stat_file_name(filename);
			ofstream gen_stat_file("Output//nsde_gen_bi_norm_stat_filename.txt");
			print_file_name(gen_stat_file,func_type,algo_type,filename);
		}
	}
	if ( is_ep(algo_type) )
	{
		static ofstream gen_eta_filelist("Output//gen_eta_filelist.txt");
		shared_ptr<fep_alg> pFEP_alg=dynamic_pointer_cast<fep_alg>(pAlg);// downcast
		string filename=gen_avg_eta_file_name(str_common,str_seq_desc);
		pFEP_alg->set_gen_eta_file_name(filename);
		print_file_name(gen_eta_filelist,func_type,algo_type,filename);
	}
}// end function algo_spec_output

void scheduler::print_stat_file(const vector<shared_ptr<stringstream> >& pvec_line)
{
	int line_count=pvec_line.size();
	int i;
	for ( i=0;i<line_count;i++ )
	{
		run_stat_file<<(*pvec_line[i]).str();
	}
}

int scheduler::batch_run(int max_thread)
{
	// every single line of run stat file
	vector<shared_ptr<stringstream> > pvec_stat_line;// dummy file buffer to write,this buffer is to keep the run stat data output FIFO ordered for post-computation analysis convenience
	// using pointer to work around the problem that stringstream is non-copyable
	{
		try 
		{
			// load batch run script file
			batch_conf batch_conf;
			batch_conf.load_conf(m_batch_conf_path);

			vector<shared_ptr<func_base> > vec_ptr_func;// function object base class pointer
			vector<shared_ptr<para_base> > vec_ptr_para;// parameter object base class pointer
			vector<shared_ptr<alg_base> > vec_ptr_alg;// algorithm object base class pointer

			int num_alg;// number of algorithm wait to test
			num_alg=batch_conf.tasks.size();
			vec_ptr_func.resize(num_alg);
			vec_ptr_para.resize(num_alg);
			vec_ptr_alg.resize(num_alg);
			// pool algo_pool(num_alg<max_thread?num_alg:max_thread);// thread pool of algorithms
			pool algo_pool(max_thread); // thread pool of algorithms
			int cur_task_count=0;

			// batch run of specified algorithm
			int alloc_algo_res;
			int algo_type;
			string conf_path;
			bool first_task_run=true;// first time run of batch run
			bool seq_task_first_time=true;// first time run of single sequential value task
			bool different_algo=true;
			int pre_algo_type=0;
			int i,j;

			vec_task::iterator itr=batch_conf.tasks.begin();
			for ( i=0;itr!=batch_conf.tasks.end();++itr,i++ )
			{
				algo_type=(*itr).algo_code;
				conf_path=(*itr).conf_path;
				// allocate specific algorithm object dynamically
				alloc_algo_res=allocate_alg(algo_type,vec_ptr_alg[i],conf_path);
				allocate_para(algo_type,vec_ptr_para[i],conf_path);
				int load_res;
				// load configuration file data
				load_res=vec_ptr_para[i]->read_para(conf_path);
				if ( -1==load_res )
				{
					cout << "read config file ERROR!" 
						<< "\n";
					return -1;
				}

				different_algo=(pre_algo_type!=algo_type);
				if ( different_algo )
					seq_task_first_time=true;
				int func_type=vec_ptr_para[i]->get_obj_func();// test problem type
				int num_dim=vec_ptr_para[i]->get_dim();// dimension of test problem
				int stop_type=vec_ptr_para[i]->get_stop_type();// termination criterion value
				// allocate Objective Function object dynamically
				int alloc_fun_res=allocate_func(func_type,vec_ptr_func[i],num_dim);
				bool max_gen_set=has_max_gen(stop_type);

				// seed random number generator for EACH task
				gen.seed(static_cast<uint32_t>(time(0)));

				// generate output file name generic prefix
				string str_out_prefix_common;
				str_out_prefix_common=gen_out_prefix_common(algo_type,vec_ptr_func[i],vec_ptr_para[i]);
				if ( first_task_run ) 				// print overall batch run stat file title once
				{
					print_algo_stat_title(run_stat_file,max_gen_set);
					print_common_filelist_title(com_filelist);
					first_task_run=false;
				}
				// start sequential parameter valus SPECIFIC processing
				bool has_seq_val=vec_ptr_para[i]->has_seq_value();
				if ( has_seq_val )
				{
					// open sequential parameter value stat file BEFORE all combination is scheduled
					string para_stat_path=gen_para_stat_file_name(str_out_prefix_common);
					// using heap memory and smart pointer to shared the data to sub-threads
					locked_para_stat_os lps_os={shared_ptr<ostream>(new ofstream(para_stat_path.c_str())),shared_ptr<mutex>(new mutex)};
					shared_ptr<ostream> &ppara_stat_file=lps_os.pos; 
					if ( seq_task_first_time )
					{
						print_para_stat_title(*ppara_stat_file,gen_seq_var_name_str(vec_ptr_para[i]),gen_seq_var_num_str(vec_ptr_para[i]),max_gen_set);
						print_file_name(para_stat_filelist,func_type,algo_type,para_stat_path); 
						seq_task_first_time=false;
					}
					// start sequential paramter value processing
					const vector<string>& seq_var_name=vec_ptr_para[i]->get_seq_var_name_vec();
					const int seq_var_size=seq_var_name.size();
					const vector<int>& seq_dim_size=vec_ptr_para[i]->get_seq_dim_size_vec();
					int comb_num=1;
					// calc total combination number
					for ( j=0;j<seq_var_size;j++ )
						comb_num *= seq_dim_size[j];
					// generate all parameter value combination index
					vector<vector<int> > para_comb_idx;
					para_comb_idx.push_back(vector<int>(seq_var_size,0));
					int cur_dim;
					int last_dim=seq_var_size-1;
					// using this method since seq_var_size is NOT constant,
					// hence we can't use for loop to iterate
					for ( j=1;j<comb_num;j++ )
					{
						cur_dim=last_dim;
						vector<int> cur_comb_idx=para_comb_idx.back(); 
						// find approriate dim to increment index
						while ( cur_comb_idx[cur_dim] == seq_dim_size[cur_dim]-1 )
							cur_dim--;
						cur_comb_idx[cur_dim]++;// 'carry propagation'
						if ( cur_dim<last_dim )
						{
							int j;
							for ( j=cur_dim+1;j<seq_var_size;j++ )
								cur_comb_idx[j]=0;
						}
						para_comb_idx.push_back(cur_comb_idx);
					}

					// start combination schedule
					string cur_para;
					int val_idx;
					vector<int> cur_comb_idx(seq_var_size);
					int k;
					// allocate EXTRA algorithm objects andd function objects
					// Parameters object are SHARED between different parameter value combination(const)
					vector<shared_ptr<alg_base> > vec_seq_alg;
					vector<shared_ptr<func_base> > vec_seq_func;
					vector<shared_ptr<para_base> > vec_seq_para;
					vec_seq_alg.resize(comb_num);
					vec_seq_func.resize(comb_num);
					vec_seq_para.resize(comb_num);
					vec_seq_alg[0]=vec_ptr_alg[i];
					vec_seq_func[0]=vec_ptr_func[i];
					vec_seq_para[0]=vec_ptr_para[i];
					for ( k=1;k<comb_num;k++ )
					{
						allocate_alg(algo_type,vec_seq_alg[k],conf_path);
						allocate_func(func_type,vec_seq_func[k],num_dim);
						allocate_para(algo_type,vec_seq_para[k],conf_path);
						(*vec_seq_para[k])=(*vec_ptr_para[i]);// copy parameter values
					}
					for ( k=0;k<comb_num;k++ )
					{
						// alter the sequential parameter value for grouped run
						cur_comb_idx=para_comb_idx[k];
						int j;
						for ( j=0;j<seq_var_size;j++ )
						{
							cur_para=seq_var_name[j];
							val_idx=cur_comb_idx[j];
							vec_seq_para[k]->set_cur_seq_val_idx(cur_para,val_idx);
						}
						// all combination values are set
						// set task description text for output
						string str_out_prefix_para;
						str_out_prefix_para=gen_out_prefix_para(vec_seq_para[k],"_");
						vec_seq_alg[k]->set_para_val_desc(gen_seq_val_str(vec_seq_para[k]));
						string str_para_stat=gen_out_prefix_para(vec_seq_para[k],",");
						// output common stat file path list for post-computation plotting and analysis
						common_path com_out_path;// common stat file output
						gen_common_file_path(com_out_path,str_out_prefix_common,str_out_prefix_para);
						vec_seq_alg[k]->set_com_file_path(com_out_path);
						print_common_filelist(com_filelist,func_type,algo_type,com_out_path);
						print_para_stat_All_Best_Val_Path(para_stat_all_bst_val_filelist,func_type,str_para_stat,com_out_path.all_bst_val_path);// sequential parameters tuning SPECIFIC
						//  set algorithm specific stat file name
						algo_spec_output(func_type,algo_type,vec_seq_alg[k],str_out_prefix_common,str_out_prefix_para);

						// run the algo according to single parameter value combination
						vec_seq_alg[k]->set_eval_func(vec_seq_func[k]);
						vec_seq_alg[k]->load_parameter(vec_seq_para[k]);// Loading algorithm parameters from outer source

						vec_seq_alg[k]->initialize();
						// progress_timer p_t;// timer of specified algorithm
						cur_task_count++;
						output_running_prompt(cur_task_count,vec_seq_para[k]);
						pvec_stat_line.push_back(shared_ptr<stringstream>(new stringstream));// expand size by 1
						shared_ptr<stringstream> pcur_line=pvec_stat_line.back();
						shared_ptr<boost::mutex> pmutex=lps_os.pmut;
						algo_pool.schedule( 
							boost::bind(thread_fun,vec_seq_alg[k],func_type,algo_type,
							max_gen_set,has_seq_val,pcur_line,
							ppara_stat_file,pmutex ) 
							);
					}// for every single value combination
				}
				else
				{
					// NO value combination:no need to schedule
					common_path com_out_path;
					gen_common_file_path(com_out_path,str_out_prefix_common,"");
					vec_ptr_alg[i]->set_com_file_path(com_out_path);
					print_common_filelist(com_filelist,func_type,algo_type,com_out_path);
					algo_spec_output(func_type,algo_type,vec_ptr_alg[i],str_out_prefix_common,"");
					// run the algo according to specified parameter value combination
					vec_ptr_alg[i]->set_eval_func(vec_ptr_func[i]);
					vec_ptr_alg[i]->load_parameter(vec_ptr_para[i]);// Loading algorithm parameters from outer source

					vec_ptr_alg[i]->initialize();

					// progress_timer p_t;// timer of specified algorithm
					cur_task_count++;
					output_running_prompt(cur_task_count,vec_ptr_para[i]);
					pvec_stat_line.push_back(shared_ptr<stringstream>(new stringstream));// expand size by 1,use heap memory to shared the data between main thread and sub-thread 
					shared_ptr<stringstream> pcur_line=pvec_stat_line.back();
					algo_pool.schedule(
							boost::bind(thread_fun, vec_ptr_alg[i], func_type,
									algo_type, max_gen_set, has_seq_val,
									pcur_line, shared_ptr<ostream>(),
									shared_ptr<boost::mutex>()));
				}// if ( has_seq_val )
				pre_algo_type=algo_type;
			}// for every single algorithm task
		}// try
		catch(bad_cast &bc)
		{
			cout<<"bad_cast caught in batch_run:" 
				<< bc.what() << "\n";
			return -1;
		}
		catch(logic_error &le)
		{
			cout<<"logic_error caught in batch_run:" 
				<< le.what() << "\n";
			return -1;
		}
		catch(exception& e) 
		{
			cout<<"exception caught in batch_run:" 
				<<e.what()<< "\n";
			return -1;
		}
	}
	// all task are finished
	print_stat_file(pvec_stat_line);// keep run stat file output data ordered according to FIFO
	return 0;
}// end function batch_run

bool is_int(const any & operand)
{
	return operand.type() == typeid(int);
}

bool is_double(const any & operand)
{
	return operand.type() == typeid(double);
}

string scheduler::gen_func_eval_name(const shared_ptr<func_base> &pFunc)
{
	string str_func;
	str_func+=typeid(*pFunc).name();
	string str_tar("class benchmark::func_");
	int rm_pos=str_func.find(str_tar);
	if ( string::npos!=rm_pos )
		str_func.replace(rm_pos,str_tar.size(),"");// delete "benchmark::class func_"
	return str_func;
}

string scheduler::gen_seq_var_name_str(const shared_ptr<para_base> &ppara)
{
	string str_seq_var_name("");
	if ( ppara->has_seq_value() )
	{
		// find all sequential parameters to append
		const vector<string>& seq_var_name=ppara->get_seq_var_name_vec();
		int seq_var_size=seq_var_name.size();
		int i;
		for ( i=0;i<seq_var_size;i++ )
		{
			string para=seq_var_name[i];
			str_seq_var_name += para;
			str_seq_var_name += "\t";
		}
	}
	return str_seq_var_name;
}

string scheduler::gen_seq_val_str(const shared_ptr<para_base> &ppara)
{
	string str_seq_var_val("");
	if ( ppara->has_seq_value() )
	{
		// find all sequential parameters to append
		const vector<string>& seq_var_name=ppara->get_seq_var_name_vec();
		int seq_var_size=seq_var_name.size();
		int i;
		for ( i=0;i<seq_var_size;i++ )
		{
			string para=seq_var_name[i];
			any para_val=ppara->get_cur_seq_val(para);
			stringstream str_val;
			if ( !is_int(para_val) )
				str_val<<any_cast<double>(para_val);
			else
				str_val<<any_cast<int>(para_val);
			str_seq_var_val += str_val.str();
			if ( i!=seq_var_size-1 )
				str_seq_var_val += "\t";
		}
	}
	return str_seq_var_val;
}

string scheduler::gen_seq_var_num_str(const shared_ptr<para_base> &ppara)
{
	string str_seq_var_num("");
	if ( ppara->has_seq_value() )
	{
		// find all sequential parameters to append
		const vector<string>& seq_var_name=ppara->get_seq_var_name_vec();
		int seq_var_size=seq_var_name.size();
		str_seq_var_num += "seq_para_var_num=";
		str_seq_var_num += lexical_cast<string>(seq_var_size);
	}
	return str_seq_var_num;
}

string scheduler::gen_out_prefix_common(int algo_type,const shared_ptr<func_base> &pFunc,const shared_ptr<para_base> &ppara)
{
	string str_out_prefix;
	m_algo_name=algo_name[algo_type-1];
	str_out_prefix=m_algo_name;
	str_out_prefix+="_";
	m_func_eval_name=gen_func_eval_name(pFunc);
	str_out_prefix+=m_func_eval_name;
	str_out_prefix+="_run=";
	std::stringstream buf_run;
	buf_run<<ppara->get_max_run();
	str_out_prefix+=buf_run.str();// type
	str_out_prefix+="_dim=";
	std::stringstream buf_dim;
	buf_dim<<ppara->get_dim();;
	str_out_prefix+=buf_dim.str();

	return str_out_prefix;
}

string scheduler::gen_out_prefix_para(const shared_ptr<para_base> &ppara,string delimiter)
{
	if ( ppara->has_seq_value() )
	{
		// find all sequential parameters to append
		const vector<string>& seq_var_name=ppara->get_seq_var_name_vec();
		int seq_var_size=seq_var_name.size();
		int i;
		m_para_out_prefix.clear();
		for ( i=0;i<seq_var_size;i++ )
		{
			string para=seq_var_name[i];;
			m_para_out_prefix += para;
			m_para_out_prefix += "=";
			any para_val=ppara->get_cur_seq_val(para);
			stringstream str_val;
			if ( !is_int(para_val) )
			{
				str_val.precision(4);
				str_val<<any_cast<double>(para_val);
			}
			else
				str_val<<any_cast<int>(para_val);
			m_para_out_prefix += (str_val.str());
			if ( i!=seq_var_size-1 )
				m_para_out_prefix += delimiter;
		}
	}
	return m_para_out_prefix;
}

string scheduler::gen_out_para_val_str(const shared_ptr<para_base> &ppara)
{
	if ( ppara->has_seq_value() )
	{
		// find all sequential parameters to append
		const vector<string>& seq_var_name=ppara->get_seq_var_name_vec();
		int seq_var_size=seq_var_name.size();
		int i;
		m_para_val_prefix.clear();
		for ( i=0;i<seq_var_size;i++ )
		{
			string para=seq_var_name[i];;
			any para_val=ppara->get_cur_seq_val(para);
			stringstream str_val;
			if ( !is_int(para_val) )
			{
				str_val.precision(6);
				str_val<<any_cast<double>(para_val);
			}
			else
				str_val<<any_cast<int>(para_val);
			m_para_val_prefix += (str_val.str());
			m_para_val_prefix += "\t";
		}
	}
	return m_para_val_prefix;
}

void scheduler::output_running_prompt(int counter,const shared_ptr<para_base> &ppara)
{
	string str_upper;
	int (*pf)(int)=toupper;
	transform(m_algo_name.begin(),m_algo_name.end(),back_inserter(str_upper),pf);// transform to upper case
	str_upper=str_upper.substr(0,str_upper.size());
	bool has_seq_value=ppara->has_seq_value();
	string str_para_desc;
	// make parameters value combination more human readable
	if ( has_seq_value )
		str_para_desc=gen_out_prefix_para(ppara,",");
	io_mutex.lock();
	cout<<counter
		<<"."
		<<" "
		<<str_upper
		<<" "
		<<"Algorithm is calculating"
		<<" "
		<<"function"
		<<" "
		<<m_func_eval_name
		<<"("
		<<ppara->get_dim()
		<<" "
		<<"dimensions"
		<<")"
		<<" "
		<<"in"
		<<" "
		<<ppara->get_max_run()
		<<" "
		<<"run";
	if ( has_seq_value )
	{
		cout<<" "
			<<"("
			<<str_para_desc
			<<")";
	}
	cout<<"\n";
	io_mutex.unlock();
}// end function output_running_prompt


