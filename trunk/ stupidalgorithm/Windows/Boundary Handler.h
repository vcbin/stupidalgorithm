#ifndef STUPIDALGO_BOUNDARY_HANDLERS
#define STUPIDALGO_BOUNDARY_HANDLERS

namespace boundary_condition
{
	// OB handle class
	struct bound_handle
	{
		virtual void handle(double &velocity) {}
	};
	// reflecting boundary condition:change direction
	struct bh_reflect:public bound_handle
	{
		void handle(double &velocity) { velocity *= -1.0; }
	};
	// absorbing boundary condition, velocity=0.0;
	struct bh_absorb:public bound_handle
	{
		void handle(double &velocity) { velocity = 0.0; }
	};
	struct bh_damp:public bound_handle
	{
		void handle(double &velocity);
	};
}// end namespace boundary_condition

#endif