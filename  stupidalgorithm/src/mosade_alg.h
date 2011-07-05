#ifndef STUPIDALGO_MOP_MOSADE
#define STUPIDALGO_MOP_MOSADE

#include "de_alg.h"
#include "mop_common.h"

namespace de
{
	namespace mosade
	{
		class mosade_alg:public de_alg,public mop::mop_common
		{
		public:
			mosade_alg(std::string conf_path):
			  de_alg(conf_path) {}
		protected:
			int run();
			void crowd_tour(const mop::elite_archive &ext_archive,
							individual &trial,individual &ind,int index,int trunc_type,
							int neighbor_num=2);
			void bound_check(double &tri_x,double orig_x,int dim);// mosade SPECIFIC boundaries checking and reparing routine from PDE
			void record_gen_vals(alg_stat &alg_stat,int cur_run);
			void write_stat_vals();
		};//end class mosade_alg

	}// end namespace mosade
}// end namespace de

#endif