#include <iostream>
#include <boost/progress.hpp>

#include "Scheduler.h"

using namespace boost::program_options;

using boost::progress_timer;
using std::cout;
using std::string;
using std::exception;


int main(int argc, char *argv[])
{
	string batch_conf_path;
	int max_thread;
	const string ver_str="version number:2.0.10\nauthor:Jerry Liu(john.rockmania@gmail.com),Long Li\n";
	try 
	{
		// Declare two groups of options.
		options_description opt_desc("");
		opt_desc.add_options()
			("algo_desc,d","algorithm number description\n"
			"Values:\n"
			" 1: standard PSO algorithm\n"
			" 2: mPSO algorithm\n"
			" 3: arPSO algorithm\n"
			" 4: dPSO algorithm\n"
			" 5: my dPSO algorithm\n"
			" 6: PSObc algorithm\n"
			" 7: basic DE algorithm\n"
			" 8: barebones DE algorithm\n"
			" 9: self-adaptive DE algorithm\n"
			" 10: my self-adaptive DE algorithm\n"
			" 11: self-adaptive pareto DE algorithm\n"
			" 12: self-adaptive jDE algorithm\n"
			" 13: fast Evolutionary Programming algorithm\n"
			" 14: improved fast Evolutionary Programming algorithm\n"
			" 15: diversity guided Evolutionary Algorithm\n"
			" 16: Enhanced Differential Evolution Algorithm\n"
			" 17: DE/EDA Algorithm\n"
			" 18: Eda-directed Mutation Differential Evolution Algorithm\n"
			" 19: basic Estimation of Distribution Algorithm\n"
			" 20: bi-Normal Self-adaptive DE algorithm\n"
			" 21: Multi-Objective Self-Adaptive Differential Evolution Algorithm\n"
			" 22: Multi-Objective Differential Evolution Algorithm\n"
			);

		options_description cmdline_opts("General options");
		cmdline_opts.add_options()
			("algo_desc,d","display algorithm code number\n")
			("config,c",value<string>(&batch_conf_path)->default_value("batch_run.txt"), "name of batch run configuration file")
			("max_thrd,t",value<int>(&max_thread)->default_value(1), "allowed maximum thread number")
			("help,h", "produce help message")
			("version,V", "output the version number")
			;

		variables_map opt_vm;
		store(parse_command_line(argc, argv, cmdline_opts), opt_vm);

		if ( opt_vm.count("algo_desc") ) 
		{
			cout << opt_desc;
			return 0;
		}
		if ( opt_vm.count("help") ) 
		{
			cout << cmdline_opts;
			return 0;
		}
		if (opt_vm.count("version")) 
		{
			cout<<ver_str;
			return 0;
		}
		if (opt_vm.count("max_thrd"))
			max_thread=opt_vm["max_thrd"].as<int>();

		if (opt_vm.count("config"))
			batch_conf_path=opt_vm["config"].as<string>();
		else
			return 0;

		progress_timer total_timer;
		{
			scheduler sch(batch_conf_path);
			sch.batch_run(max_thread);
			cout<<"total time:"<<" ";
		}

	}// try
	catch(exception& e) 
	{
		cout << e.what() << "\n";
		return EXIT_FAILURE;
	}
#ifdef WIN32
	system("PAUSE");// pause command only works on windows OS
#endif
	return 0;
}
