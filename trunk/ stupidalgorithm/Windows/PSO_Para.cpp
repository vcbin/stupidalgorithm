#include <fstream>
#include <iostream>

#include "pso_para.h"

using namespace boost;
using namespace boost::program_options;
namespace po = boost::program_options;

using std::ifstream;
using std::string;
using std::cout;
using std::stringstream;
using std::map;
using std::make_pair;
using std::logic_error;
using std::vector;

using boost::split;
using boost::is_any_of;

namespace pso
{
	int pso_para::read_para(string conf_path)
	{
		string path=conf_path.empty()?m_conf_path:conf_path;
		try
		{
			para_base::read_para(conf_path);
			// end of common parameters for all algorithms

			parse_seq_var(m_ext_vm,"Vmax_type",vmax_type);

			parse_seq_var(m_ext_vm,"Vmax_cof",vmax_cof);

			parse_seq_var(m_ext_vm,"ob_type",ob_type);

			parse_seq_var(m_ext_vm,"omega_max",w_max);

			parse_seq_var(m_ext_vm,"omega_min",w_min);

			parse_seq_var(m_ext_vm,"d_l",d_l);

			parse_seq_var(m_ext_vm,"d_h",d_h);

			parse_seq_var(m_ext_vm,"cr1",cr1);

			parse_seq_var(m_ext_vm,"cr2",cr2);

			if ( stoptype_gen==stop_type )
				max_gen=static_cast<int>(stop_val);
			else if ( stoptype_eval==stop_type )
				max_gen=static_cast<int>(stop_val)/pop_size;
			else
				max_gen=0;

			// generate sequential variable size array
			size_t i;
			for ( i=0;i<m_seq_var_name.size();i++ )
				m_seq_dim_size.push_back(m_seq_var_map[m_seq_var_name[i]].size());
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

}// end namespace pso
