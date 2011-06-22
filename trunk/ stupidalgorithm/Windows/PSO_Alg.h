#ifndef STUPIDALGO_PSO_STD_ALGORITHM
#define STUPIDALGO_PSO_STD_ALGORITHM

// #include <list>
#include <boost/progress.hpp>
#include <boost/timer.hpp>
#include <boost/shared_ptr.hpp>

#include <fstream>

#include "Rand_Val.h" // declaration of necessary boost::random header

#include "pso_para.h"
#include "com_alg.h"
#include "Boundary Handler.h"

using namespace boundary_condition;

namespace pso
{
	void alloc_bound_handle(int ob_type,boost::shared_ptr<bound_handle> &pBH);

	class pso_alg:public com_alg<pso_para>,protected bound_handle
	{
	public:
		// ctor
		pso_alg(std::string conf_path):
				com_alg(conf_path) {}
		pso_alg(boost::shared_ptr<pso_para> ppara):// set algorithm parameters pointer directly
				com_alg(ppara) {}
		  // virtual dtor
		  virtual ~pso_alg() { } // for polymorphism

		  void initialize();

		  int run();
		  inline void  set_gen_vel_div_file_name( std::string str_filename ){ m_avg_vel_div_path=str_filename; }
		  void print_run_title(std::ostream &os);
		  void print_gen_stat(std::ostream &os,int cur_gen,const alg_stat &alg_stat);
	protected:
		void set_orig_pop(population &pop);
		void update_pop(population &pop);
		void reinit_x_range();
		void set_dynamic_vmax();
		void set_static_vmax();
		virtual void set_ini_max_velocity();
		void update_omega(int cur_gen);
		void update_diversity(const population &pop);
		virtual double calc_omega(int cur_gen);
		virtual void update_speed(population &pop);

		void reset_run_stat();

		void update_x_rng(const population &pop);

		boost::shared_ptr<bound_handle> m_pBH;
		enum {dynamic_rng=1,static_rng};

		d_array m_vmax;// maximum velocity of every dimenson
		d_array m_v_dim_diver;// velocity diversity of every dimenson
		double m_vel_diver;// population velocity diversity
		d_mat m_gen_vel_div_val;// aux var,generational velocity diversity of all run
		d_array m_gen_avg_vel_div_val;// average velocity diversity of all run in every generation
		void write_avg_vel_div_per_gen();
		void write_stat_vals();
		std::string m_avg_vel_div_path;
		
		void record_vel_div(int cur_run);
		void record_gen_vals(alg_stat &alg_stat,int cur_run);
		void calc_gen_avg_vals();
		val_range m_x_rng;// DYNAMIC x range of every dimenson
		std::vector<std::vector<bool> > m_outbound;// variable for boundary condition check

		double m_omega;// inertia coefficient value
	};// end class pso_alg declaration

}// namespace pso

#endif