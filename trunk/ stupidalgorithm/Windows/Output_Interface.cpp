#include "Output_Interface.h"
#include "CmdLine_Para_Type.h"

using std::ios;
using std::ostream;
using std::string;
using std::vector;
using std::copy;
using std::ostream_iterator;

using boost::shared_ptr;

void Output_Interface::Print_Run_Title(ostream &os)
{
	os<<"\n"
		<<"gen"
		<<"\t"
		<<"best_so_far"
		<<"\t"
		<<"avg"
		<<"\t"
		<<"pos_diversity"
		<<"\t"
		<<"std"
		<<"\t"
		<<"stag_gen"
		<<"\t"
		<<"eval_count";
}

void Output_Interface::Print_Sep_Line(ostream &os)
{
	int i;
	os<<"\n";
	for ( i=0;i<64;i++ )
		os<<"*";
}

void Output_Interface::Print_Gen_Stat(ostream &os,int cur_gen,const Alg_Stat &alg_stat)
{
	//// set float-point precision for output
	//os.precision(8);
	//os.setf(ios::scientific,ios::floatfield);
	os<<"\n"
		<<cur_gen
		<<" "
		<<alg_stat.gbest.obj[0];
	os.precision(8);
	// os.setf(ios::fixed,ios::floatfield);
	os<<"\t"
		<<alg_stat.avg
		<<"\t"
		<<alg_stat.pos_diver;

	//os.setf(ios::scientific,ios::floatfield);
	os<<"\t"
		<<alg_stat.std
		<<"\t"
		<<alg_stat.num_stag
		<<"\t"
		<<alg_stat.all_eval_num;
}

void Output_Interface::Print_Gen_Avg_Vals(ostream &os,const std::vector<double> &gen_vals)
{
	int i;
	int max_gen=gen_vals.size();
	for ( i=0;i<max_gen;i++ )
		print_gen_val(os,i+1,gen_vals[i]);
}

void Output_Interface::print_gen_vec(ostream &os,int cur_gen,const vector<double> &vec)
{
	os<<cur_gen<<"\t";
	int num_dims=vec.size();
	int j;
	for ( j=0;j<num_dims;j++ )
		os<<vec[j]<<"\t";
	os<<"\n";
}

void Output_Interface::Print_Gen_Avg_Vec(ostream &os,const std::vector<std::vector<double> > &gen_vec)
{
	int i;
	int max_gen=gen_vec.size();
	for ( i=0;i<max_gen;i++ )
		print_gen_vec(os,i+1,gen_vec[i]);
}

void Output_Interface::Print_Best_X(ostream &os,const Individual &bst_ind)
{
	int num_dims=bst_ind.x.size();
	os<<"\n"
		<<"best X value:"
		<<"\n";
	int i;
	for ( i=0;i<num_dims;i++ )
	{
		os<<bst_ind.x[i];
		if ( i!=num_dims-1 )
		{
			if ( 0==i%2 )
				os<<" ";
			else
				os<<"\n";
		}
	}
}

void Output_Interface::Print_Run_Stat_Title(ostream &os)
{
	os<<"\n"
		<<"\n"
		<<"best_of_all_gbest"
		<<"\t"
		<<"worst_of_all_gbest"
		<<"\t"
		<<"mean_of_gbest"
		<<"\t"
		<<"std_of_gbest"
		<<"\t"
		<<"convergence_number"
		<<"\t"
		<<"convergence_ratio"
		<<"\t"
		<<"ob_num"
		<<"\t"
		<<"average_stag_gen_on_stop";
}

void Output_Interface::Print_Run_Stat(ostream &os,const Alg_Stat &alg_stat,int max_run)
{
	// set float-point precision for output
	/*os.precision(8);
	os.setf(ios::scientific,ios::floatfield);*/

	Print_Run_Stat_Title(os);
	os<<"\n"
		<<alg_stat.bst_ind.obj[0]
		<<"\t"
		<<alg_stat.wst_ind.obj[0]
		<<"\t"
		<<alg_stat.run_avg_bst
		<<"\t"
		<<alg_stat.run_bst_std
		<<"\t"
		<<alg_stat.conv_count
		<<"/"
		<<max_run
		<<"\t";

	// percent number format
	os.setf(ios::fixed,ios::floatfield);
	os.precision(2);
	os<<(alg_stat.conv_ratio*100.0)
		<<"%"
	<<"\t"
	<<alg_stat.all_ob_num
	<<"\t"
	<<alg_stat.run_avg_stag_gen;
}

void Output_Interface::Print_Avg_Time(ostream &os,double sec)
{
	long sec_elapsed=static_cast<long>(sec);
	os<<"\n"
		<<"average_sec_per_run:"
		<<" "
		<<sec
		<<" seconds";
	long hours=static_cast<long>(sec/3600);
	long minutes=sec_elapsed%3600/60;
	long seconds=sec_elapsed%3600%60;
	if ( hours )
	{
		os<<" ("
			<<hours
			<<" hour(s)"
			<<" and "
			<<minutes
			<<" "
			<<"minute(s)"
			<<" and "
			<<seconds
			<<" second(s)"
			<<")";
	}
	if ( 0==hours && minutes )
	{
		os
			<<" ( "
			<<minutes
			<<" "
			<<"minute(s)"
			<<" and "
			<<seconds
			<<" second(s)"
			<<" )";
	}
}// end function Print_Avg_Time

void Output_Interface::write_initial_pop(ostream &os,const Population &pop)
{
	int pop_size=pop.size();
	if ( pop_size>0 )
	{
		int num_dims=pop[0].x.size();
		int i,j;
		for ( i=0;i<pop_size;i++ )
		{
			for ( j=0;j<num_dims;j++ )
			{
				os<<pop[i].x[j];
				os<<" ";
			}
			os<<"\n";
		}
	}
}// end function print_initial_pop

void Output_Interface::Print_All_Gbest_Vals(ostream &os,const vector<double> &gbest_vec)
{
	copy(gbest_vec.begin(),gbest_vec.end(),ostream_iterator<double>(os," "));
}

void Print_Algo_Stat_Title(ostream &os,bool has_max_gen)
{
	os<<"func_code"
		<<"\t"
		<<"algo_code"
		<<"\t"
		<<"best"
		<<"\t"
		<<"worst"
		<<"\t"
		<<"avg_val"
		<<"\t"
		<<"std"
		<<"\t"
		<<"conv_count"
		<<"\t"
		<<"conv_ratio"
		<<"\t"
		<<"avg_conv_eval_num"
		<<"\t"
		<<"avg_time"
		<<"\t"
		<<"eval_count"
		<<"\t"
		<<"ob_count"
		<<"\t"
		<<"avg_stag_gen_on_stop";
	if ( has_max_gen )
		os<<"\t"
		<<"avg_gen";
}

void Print_Para_Stat_Title(ostream &os,string para_var_list,string para_num_str,bool has_max_gen)
{
	os<<para_var_list
		<<"\t"
		<<"best"
		<<"\t"
		<<"worst"
		<<"\t"
		<<"avg_val"
		<<"\t"
		<<"std"
		<<"\t"
		<<"conv_count"
		<<"\t"
		<<"conv_ratio"
		<<"\t"
		<<"avg_conv_eval_num"
		<<"\t"
		<<"avg_time"
		<<"\t"
		<<"eval_count"
		<<"\t"
		<<"ob_count"
		<<"\t"
		<<"avg_stag_gen_on_stop";
	if ( has_max_gen )
		os<<"\t"
		<<"avg_gen";
	os<<"\n"
		<<para_num_str;
}

void Print_Common_Filelist_Title(ostream &os)
{
	os<<"func_type"<<"\t"
		<<"algo_type"<<"\t"
		<<"all_bst_val_path"<<"\t"
		<<"avg_bst_path"<<"\t"
		<<"avg_div_path"<<"\t"
		<<"avg_rad_path"
		<<"\n";
}
