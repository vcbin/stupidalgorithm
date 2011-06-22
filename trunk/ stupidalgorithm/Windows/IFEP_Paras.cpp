#include "IFEP_Paras.h"

using namespace boost;
using namespace boost::program_options;
namespace po = boost::program_options;

using std::ifstream;
using std::string;
using std::cout;

namespace ep
{
	int IFEP_Paras::ReadPara(string conf_path)
	{
		try
		{
			variables_map vm;
			string path=conf_path.empty()?m_path:conf_path;
			ifstream ifs(path.c_str());
			if (!ifs)
				return -1;
			store(parse_config_file(ifs,general), vm);
			notify(vm);

			double up_bnd,low_bnd;
			if ( vm.count("val_lower_bound") )
				low_bnd=vm["val_lower_bound"].as<double>();

			if ( vm.count("val_upper_bound") )
				up_bnd=vm["val_upper_bound"].as<double>();

			if ( vm.count("dimension") )
			{
				dims=vm["dimension"].as<int>();
				val_bounds.resize(dims);
				int i;
				for ( i=0;i<dims;i++ )
				{
					val_bounds[i].min_val=low_bnd;
					val_bounds[i].max_val=up_bnd;
				}
			}

			if ( vm.count("ini_lower_bound") )
				low_bnd=vm["ini_lower_bound"].as<double>();

			if ( vm.count("ini_upper_bound") )
				up_bnd=vm["ini_upper_bound"].as<double>();

			ini_bounds.resize(dims);
			int i;
			for ( i=0;i<dims;i++ )
			{
				ini_bounds[i].min_val=low_bnd;
				ini_bounds[i].max_val=up_bnd;
			}

			if ( vm.count("stop_type") )
				stop_type=vm["stop_type"].as<int>();

			if ( vm.count("stop_threshold") )
				stop_val=vm["stop_threshold"].as<double>();

			if ( vm.count("pop_size") )
				pop_size=vm["pop_size"].as<int>();

			if ( 1==stop_type )
				max_gen=static_cast<int>(stop_val);
			else if ( 2==stop_type )
				max_gen=static_cast<int>(stop_val)/pop_size;
			else
				max_gen=0;

			/*if ( vm.count("algo") )
			alg_type=vm["algo"].as<int>();*/
			if ( vm.count("func_type") )
				func_type=vm["func_type"].as<int>();

			if ( vm.count("run") )
				run=vm["run"].as<int>();
			// end of common parameters for all algorithms

			if ( vm.count("tour_size") )
				tour_size=vm["tour_size"].as<int>();
		}
		catch ( const po::multiple_occurrences& e ) 
		{
			std::cerr << e.what() 
				<< " from option: " 
				<< e.get_option_name() 
				<< std::endl;
			return -1;
		}
		catch (std::exception &err)
		{
			cout<<"\n"
				<<err.what()
				<<"\n";
			return -1;
		}
		return 0;
	}// end function LoadPara
}// namespace ep