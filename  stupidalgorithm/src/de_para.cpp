#include "de_para.h"

using namespace boost;
using namespace boost::program_options;
namespace po = boost::program_options;

using std::ifstream;
using std::string;
using std::cout;
using std::endl;

namespace de
{
	int de_para::read_para(string conf_path)
	{
		string path=conf_path.empty()?m_conf_path:conf_path;
		try
		{
			para_base::read_para(conf_path);
			// end of common parameters for all algorithms
			parse_seq_var(m_ext_vm,"pr_val",pr_val);

			parse_seq_var(m_ext_vm,"f_val",f_val);

			parse_seq_var(m_ext_vm,"pr_mean",pr_mean);

			parse_seq_var(m_ext_vm,"pr_sigma",pr_sigma);

			parse_seq_var(m_ext_vm,"f_mean",f_mean);

			parse_seq_var(m_ext_vm,"f_sigma",f_sigma);

			parse_seq_var(m_ext_vm,"ini_f_type",ini_f_type);

			parse_seq_var(m_ext_vm,"ini_f_uni_low_bound",ini_f_uni_low_bnd);

			parse_seq_var(m_ext_vm,"ini_f_uni_up_bound",ini_f_uni_up_bnd);

			parse_seq_var(m_ext_vm,"ini_f_mean",ini_f_mean);

			parse_seq_var(m_ext_vm,"ini_f_sigma",ini_f_sigma);

			parse_seq_var(m_ext_vm,"learn_p",learn_p);

			parse_seq_var(m_ext_vm,"f_per_dim",f_per_dim);

			parse_seq_var(m_ext_vm,"f_lower_bound",f_low_bnd);

			parse_seq_var(m_ext_vm,"f_upper_bound",f_up_bnd);

			parse_seq_var(m_ext_vm,"tau_1",tau_1);

			parse_seq_var(m_ext_vm,"tau_2",tau_2);

			parse_seq_var(m_ext_vm,"pr_stra",pr_stra);

			parse_seq_var(m_ext_vm,"delta",delta);

			parse_seq_var(m_ext_vm,"m_ratio",m_ratio);

			parse_seq_var(m_ext_vm,"max_archive",max_archive);

			// generate sequential variable size array
			size_t i;
			for ( i=0;i<m_seq_var_name.size();i++ )
				m_seq_dim_size.push_back(m_seq_var_map[m_seq_var_name[i]].size());
		}
		catch ( const po::multiple_occurrences& e ) 
		{
			std::cerr << e.what() 
				<< " from option: " 
				//<< e.get_option_name() 
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
}// namespace de
