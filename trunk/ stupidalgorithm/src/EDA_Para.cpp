#include "EDA_Para.h"

using namespace boost;
using namespace boost::program_options;
namespace po = boost::program_options;

using std::ifstream;
using std::string;
using std::cout;
using std::endl;
using std::cerr;

namespace eda
{
	int eda_para::read_para(string conf_path)
	{
		string path=conf_path.empty()?m_conf_path:conf_path;
		try
		{
			para_base::read_para(conf_path);
			// end of common parameters for all algorithms

			parse_seq_var(m_ext_vm,"m_ratio",m_ratio);

			// generate sequential variable size array
			size_t i;
			for ( i=0;i<m_seq_var_name.size();i++ )
				m_seq_dim_size.push_back(m_seq_var_map[m_seq_var_name[i]].size());
		}
		catch ( const po::multiple_occurrences& e ) 
		{
			cerr << e.what() 
				<< " from option: " 
				<< e.get_option_name() 
				<<" "
				<<"in "
				<<path
				<< endl;
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

}// namespace eda