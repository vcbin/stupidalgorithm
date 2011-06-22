#include <fstream>
#include <sstream>
#include <string>
#include <iostream>

#include "para_base.h"

using namespace boost::program_options;
namespace po = boost::program_options;

using std::string;
using std::vector;
using std::ifstream;
using std::stringstream;
using std::exception;
using std::logic_error;
using std::cout;

using boost::split;
using boost::is_any_of;
using boost::lexical_cast;
using boost::any;

int para_base::read_para(string conf_path)
{
	string path=conf_path.empty()?m_conf_path:conf_path;
	try
	{
		ifstream ifs(path.c_str());
		if (!ifs)
			return -1;
		store(parse_config_file(ifs,general), m_ext_vm);
		notify(m_ext_vm);

		if ( m_ext_vm.count("func_type") )
			func_type=m_ext_vm["func_type"].as<int>();

		parse_seq_var(m_ext_vm,"vtr",vtr);

		if ( m_ext_vm.count("run") )
			run=m_ext_vm["run"].as<int>();

		if ( m_ext_vm.count("pop_size") )
			pop_size=m_ext_vm["pop_size"].as<int>();

		if ( m_ext_vm.count("stop_type") )
			stop_type=m_ext_vm["stop_type"].as<int>();

		if ( m_ext_vm.count("stop_threshold") )
			stop_val=m_ext_vm["stop_threshold"].as<double>();

		if ( stoptype_gen==stop_type )
			max_gen=static_cast<int>(stop_val);
		else if ( stoptype_eval==stop_type )
			max_gen=static_cast<int>(stop_val/pop_size);

		if ( m_ext_vm.count("ini_type") )
			ini_type=m_ext_vm["ini_type"].as<int>();

		if ( m_ext_vm.count("ini_pop_file") )
			ini_pop_path=m_ext_vm["ini_pop_file"].as<string>();

		if ( m_ext_vm.count("ini_range_file") )
			ini_rng_path=m_ext_vm["ini_range_file"].as<string>();

		double val_low_bnd,val_up_bnd;
		if ( m_ext_vm.count("val_lower_bound") )
			val_low_bnd=m_ext_vm["val_lower_bound"].as<double>();

		if ( m_ext_vm.count("val_upper_bound") )
			val_up_bnd=m_ext_vm["val_upper_bound"].as<double>();

		parse_seq_var(m_ext_vm,"dimension",dims);

		double ini_low_bnd,ini_up_bnd;
		if ( m_ext_vm.count("ini_lower_bound") )
			ini_low_bnd=m_ext_vm["ini_lower_bound"].as<double>();

		if ( m_ext_vm.count("ini_upper_bound") )
			ini_up_bnd=m_ext_vm["ini_upper_bound"].as<double>();

		bool ini_rng_load=false;
		int dim_size=get_dim();
		if ( initype_rng_file==ini_type )
		{
			ini_bounds.resize(dim_size);
			load_ini_range();
			ini_rng_load=true;
		}
		else if ( initype_rnd==ini_type || initype_ortho==ini_type )
		{
			ini_bounds.resize(dim_size);
			int i;
			for ( i=0;i<dim_size;i++ )
			{
				ini_bounds[i].min_val=ini_low_bnd;
				ini_bounds[i].max_val=ini_up_bnd;
			}
		}

		if ( m_ext_vm.count("bnd_value_type") )
			bnd_val_type=m_ext_vm["bnd_value_type"].as<int>();

		if ( m_ext_vm.count("bnd_range_file") )
			bnd_rng_path=m_ext_vm["bnd_range_file"].as<string>();

		parse_seq_var(m_ext_vm,"q_number","q_number");

		val_bounds.resize(dim_size);
		if ( bndvaltype_file==bnd_val_type )
		{
			if ( ini_rng_load && bnd_rng_path==ini_rng_path )
				val_bounds=ini_bounds;
			else
				load_bnd_val_range();
		}
		else
		{
			int i;
			for ( i=0;i<dim_size;i++ )
			{
				val_bounds[i].min_val=val_low_bnd;
				val_bounds[i].max_val=val_up_bnd;
			}
		}

		if ( m_ext_vm.count("num_obj") )
			num_obj=m_ext_vm["num_obj"].as<int>();

		parse_seq_var(m_ext_vm,"trunc_type",trunc_type);

		if ( m_ext_vm.count("plot") )
			plot=m_ext_vm["plot"].as<bool>();

		if ( m_ext_vm.count("plot_script") )
			plot_script=m_ext_vm["plot_script"].as<string>();
	}
	catch (std::exception &err)
	{
		cout<<"\n"
			<<err.what()
			<<" "
			<<"in "
			<<path
			<<"\n";
		return -1;
	}
	return 0;
}// end function read_para


void para_base::load_range_file(string rng_path)
{
	ifstream rng_data(rng_path.c_str());
	int i=0;// entry load number
	int j;
	string inbuf;
	int dim_size=get_dim();
	int item_num=dim_size/2;
	int col_num=4;

	vector<double> rng_tmp(col_num);
	bool inconsist_flag=false;
	val_bounds.resize(dim_size);
	if ( rng_data.fail() )
	{
		stringstream str_err;
		str_err<<"\nopen file "
			<<ini_rng_path
			<<" ERROR!\n"
			<<"please check file path!\n";
		throw logic_error(str_err.str());
	}
	while ( !rng_data.eof() )
	{
		getline(rng_data,inbuf);
		if ( false==is_dataLine(inbuf) )
			continue;
		if ( i==dim_size ) 
			inconsist_flag=true;
		if ( !inconsist_flag )
		{
			stringstream line(inbuf);
			for ( j=0;j<col_num;j++ )
			{
				line>>rng_tmp[j];
			}
			vector<double> x(2),y(2);
			x[0]=rng_tmp[0];
			x[1]=rng_tmp[1];
			y[0]=rng_tmp[2];
			y[1]=rng_tmp[3];
			ini_bounds[i].min_val=x[0];
			ini_bounds[i].max_val=x[1];
			ini_bounds[i+1].min_val=y[0];
			ini_bounds[i+1].max_val=y[1];
			if ( rng_data.fail() ) // read error
			{
				stringstream str_err("ERROR:Read data from file");
				str_err<<" "
					<<rng_path.c_str()
					<<" "
					<<"failed."
					<<"Check your data file for data validity."
					<<"\n";
				throw logic_error(str_err.str());
			}
		}// if consistent
		i=i+2;
	}// while !eof
	if ( inconsist_flag ) // entry count is inconsistent with pop_size in config file,data file is invalid.
	{
		stringstream str_err("ERROR:Invalid Initialization File data!");
		str_err<<"Initial Value Ranges file "
			<<rng_path.c_str()
			<<"'s"
			<<" node size="
			<<item_num
			<<","
			<<"while item number in initialization value ranges file is"
			<<" "
			<<i
			<<".";
		throw logic_error(str_err.str());
	}
}// end function load_range_file


// helper function
bool is_floatpoint(string lit_val) { return ( lit_val.find('.')!=string::npos || lit_val.find('e')!=string::npos || lit_val.find('E')!=string::npos ); }

int para_base::parse_seq_var(variables_map &m_ext_vm,string opt_name,string literal_var)
{
	try
	{
		if ( m_ext_vm.count(opt_name) )
			literal_var=m_ext_vm[opt_name].as<string>();
		// parse real value
		// remove '[]' or '()'
		if ( literal_var.empty() )
		{
			stringstream err;
			err<<"\n"
				<<"entry read error:"
				<<opt_name
				<<"\n";
			throw logic_error(err.str());
		}
		if ( '['==literal_var[0] || '('==literal_var[0] )
			literal_var=literal_var.substr(1,literal_var.size()-2);
		bool has_semicolon=(literal_var.find(':')!=string::npos);
		bool has_comma=(literal_var.find(',')!=string::npos);
		if (  has_semicolon || has_comma )
		{
			m_seq_var_map.insert( make_pair(opt_name,vector<any>()) );
			vector<string> strs; 
			split(strs, literal_var,is_any_of(":,"));
			int size=strs.size();
			vector<double> parsed_val;
			bool is_int=true;
			if ( has_semicolon )
			{
				if ( size!=3 )
				{
					stringstream str_err;
					str_err<<"\n"
						<<"ERROR in parameters"
						<<" "
						<<opt_name
						<<"."
						<<"\n";
					throw logic_error(str_err.str());
				}
				int i;
				// check parameter literal value format to help determine its type
				for ( i=0;i<size;i++ )
				{
					if ( is_floatpoint(strs[i]) )
					{
						is_int=false;
						break;
					}
				}
				for ( i=0;i<size;i++ )
					parsed_val.push_back(lexical_cast<double>(strs[i]));
				// convert to real parameter value
				double span=parsed_val[2]-parsed_val[0];
				double step_size=parsed_val[1];
				int len=static_cast<int>(span/step_size+1);
				for ( i=0;i<len;i++ )
				{
					if ( !is_int )
					{
						double cur_val=parsed_val[0]+i*step_size;
						m_seq_var_map[opt_name].push_back(cur_val);
					}
					else
					{
						int cur_val=static_cast<int>(parsed_val[0]+i*step_size);
						m_seq_var_map[opt_name].push_back(cur_val);
					}
				}
			}
			if ( has_comma )
			{
				int i;
				for ( i=0;i<size;i++ )
				{
					if ( is_floatpoint(strs[i]) )
					{
						is_int=false;
						break;
					}
				}
				for ( i=0;i<size;i++ )
					parsed_val.push_back(lexical_cast<double>(strs[i]));
				for ( i=0;i<size;i++ )
				{
					if ( is_int )
						m_seq_var_map[opt_name].push_back(static_cast<int>(parsed_val[i]));
					else
						m_seq_var_map[opt_name].push_back(parsed_val[i]);
				}
			}
			set_cur_seq_val_idx(opt_name,0);// initialize to first value
			m_seq_var_name.push_back(opt_name);
		}// has_semicolon || has_comma
		else // single value
		{
			if ( is_floatpoint(literal_var) )
				m_seq_var_map[opt_name].push_back( lexical_cast<double>(literal_var) );
			else
				m_seq_var_map[opt_name].push_back( lexical_cast<int>(literal_var) );
			set_cur_seq_val_idx(opt_name,0);
		}
	}// try
	catch (exception &err)
	{
		cout<<"\n"
			<<err.what()
			<<","
			<<"option name:"
			<<opt_name
			<<" "
			<<"in para_base::parse_seq_var"
			<<"\n";
		throw;
	}
	return 0;
}// end function parse_seq_var
