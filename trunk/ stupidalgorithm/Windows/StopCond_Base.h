#ifndef JERRYLIU_TERMINATION_CRITERION_BASE
#define JERRYLIU_TERMINATION_CRITERION_BASE

template<typename Real_Alg>
class Stop_Cond // virtual base class for stop criterion
{
public:
	Stop_Cond(Real_Alg *pself)
		:self(pself) {}
	virtual operator bool() {}
protected:
	Real_Alg *self;// pimp idiom:self pointer to container class for nested class
};

template<typename Real_Alg>
class stop_gen:public Stop_Cond<Real_Alg>
{
public:
	stop_gen(Real_Alg *pself)
		:Stop_Cond<Real_Alg>(pself) {}
	operator bool() { return self->m_cur_gen < self->m_pPara->GetMaxGen(); }

};

template<typename Real_Alg>
class stop_eval:public Stop_Cond<Real_Alg>
{
public:
	stop_eval(Real_Alg *pself)
		:Stop_Cond<Real_Alg>(pself) {}
	operator bool() { return self->m_alg_stat.eval_num < self->m_pPara->GetStopVal(); }
};

template<typename Real_Alg>
class stop_stag:public Stop_Cond<Real_Alg>
{
public:
	stop_stag(Real_Alg *pself)
		:self(pself) {}
	operator bool() 	{ return self->m_alg_stat.num_stag <= self->m_pPara->GetStopVal(); }
};

#endif