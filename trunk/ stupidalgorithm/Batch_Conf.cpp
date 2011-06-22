#include <fstream> // ifstream
#include <iostream>
#include <sstream>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>

#include "batch_conf.h"
#include "Algo_DataStruct.h"

using std::ifstream;
using std::string;
using std::vector;
using std::cerr;

using boost::lexical_cast;
using boost::bad_lexical_cast;
using boost::split;
using boost::is_any_of;


int batch_conf::load_conf(string path)
{
	if ( path=="" && m_strpath=="" )
		return -1;

	string conf_path=m_strpath.empty()? path:m_strpath;
	ifstream para_file(conf_path.c_str());
	if (!para_file)
	{
		std::cout<<"Open Parameter File "
			<<conf_path
			<<" "
			<<"Error!"
			<<std::endl;
		return -1;
	}
	string inbuf;

	string prev_path;
	int line_counter=0;
	while (!para_file.eof())
	{
		getline(para_file,inbuf);
		line_counter++;
		if ( false==is_dataLine(inbuf) ) 
			continue;

		vector<string> split_token;
		split( split_token, inbuf, is_any_of(" \t") );

		int algo;
		string conf_path;
		try
		{
			algo=lexical_cast<int>(split_token[0].c_str());
		}
		catch(bad_lexical_cast &blc_err)
		{
			cerr<<blc_err.what()
			<<" in line:"
			<<line_counter
			<<" of file:"
			<<path
			<<"\n";
			return -1;
		}
		conf_path=split_token[1];
		if ( conf_path=="-" )
			conf_path=prev_path;
		struct task_info task_tmp={algo,conf_path};
		vec_task.push_back(task_tmp);

		prev_path=conf_path;
	}// while !eof


	return 0;
}// end function load_conf
