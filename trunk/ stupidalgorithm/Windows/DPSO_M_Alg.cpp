#include "DPSO_M_Alg.h"

namespace pso
{
	namespace dpso_m
	{
		// momentum=0.5*(1/fitness)*$sum(V_ij*V_ij)
		void dpso_m_alg::update_momentum(population &pop)
		{
			size_t pop_size=pop.size();
			size_t num_dims=m_ppara->get_dim();
			unsigned i,j;
			m_moment.assign(m_moment.size(),0.0);// reinitialize to 0
			for ( i=0;i<pop_size;i++ )
			{
				double diff;
				for ( j=0;j<num_dims;j++ )
				{
					diff=m_alg_stat.delta_x[i][j];
					m_moment[i] += diff*diff;
				}
				m_moment[i] *= (1.0/pop[i].obj[0]);
				m_moment[i] *= 0.5;
			}
		}

	}// namespace dpso_m
}// namespace pso