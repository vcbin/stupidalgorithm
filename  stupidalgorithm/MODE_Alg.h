#ifndef JERRYLIU_MOP_MODE
#define JERRYLIU_MOP_MODE

#include "DE_Alg.h"
#include "MOP_Common.h"

namespace de
{
	namespace mode
	{
		class mode_alg:public de_alg,public mop::mop_common
		{
		public:
			mode_alg(std::string conf_path):
			  de_alg(conf_path) {}
		protected:
			int run();
			bool is_learn_gen(int cur_gen,int learn_p);

			void bound_check(double &tri_x,double orig_x,int dim);// mosade SPECIFIC boundaries checking and reparing routine from PDE
			void record_gen_vals(alg_stat &alg_stat,int cur_run);
			void write_stat_vals();	

			/*bi_norm_var m_bi_norm_var;
			d_array m_succ_f;
			d_array m_succ_cr;*/
		};//end class mosade_alg

	}// end namespace mosade
}// end namespace de

#endif