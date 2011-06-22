#include "GA_Paras.h"

using namespace boost;
using namespace boost::program_options;
namespace po = boost::program_options;

using std::ifstream;
using std::string;
using std::cout;

namespace ga
{
	int GA_Paras::ReadPara(string conf_path)
	{
		string path=conf_path.empty()?m_conf_path:conf_path;
		try
		{
			Para_Base::ReadPara(conf_path);
			// end of common parameters for all algorithms

			// parameters for ga
			parse_seq_var(m_ext_vm,"pr",pr);

			parse_seq_var(m_ext_vm,"pm",pm);

			parse_seq_var(m_ext_vm,"d_l",d_l);

			parse_seq_var(m_ext_vm,"d_h",d_h);
		}
		catch ( const po::multiple_occurrences& e ) 
		{
			std::cerr << e.what() 
				<< " from option: " 
				<< e.get_option_name() 
				<<" "
				<<"in "
				<<path
				<< std::endl;
			return -1;
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
	}// end function LoadPara
}// namespace de
