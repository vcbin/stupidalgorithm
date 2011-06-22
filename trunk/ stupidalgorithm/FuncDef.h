#ifndef JERRYLIU_FUNC_DEFS
#define JERRYLIU_FUNC_DEFS

#include <functional>
#include <boost/array.hpp>

#include "Algo_DataStruct.h"

namespace benchmark
{
	class func_base:public std::unary_function<individual&,void>
	{
	public:
		func_base(size_t dims)
			:m_num_dims(dims) {}
		virtual void operator() (individual &ind){}
		virtual bool is_convergent(double fit,double vtr)=0;
		virtual ~func_base() {}
	protected:
		int m_num_dims;
	};

	class func_Sphere :public func_base
	{
	public:
		func_Sphere(size_t dims)
			:func_base(dims) {}
		void operator() (individual &ind);
		bool is_convergent(double fit,double vtr) { return (fit<=vtr || fit==optimum); }
	private:
		static const double optimum;
	};

	class func_Griewank :public func_base
	{
	public:
		func_Griewank(size_t dims)
			:func_base(dims) {}
		void operator() (individual &ind);
		bool is_convergent(double fit,double vtr) { return (fit<=vtr || fit==optimum); }
	private:
		static const double optimum;
	};

	class func_AckeyF1:public func_base
	{
	public:
		func_AckeyF1(size_t dims)
			:func_base(dims) {}
		void operator() (individual &ind);
		bool is_convergent(double fit,double vtr) { return (fit<=vtr || fit==optimum); }
	private:
		static const double optimum;
	};

	class func_Rastrigin :public func_base
	{
	public:
		func_Rastrigin(size_t dims)
			:func_base(dims) {}
		void operator() (individual &ind);
		bool is_convergent(double fit,double vtr) { return (fit<=vtr || fit==optimum); }
	private:
		static const double optimum;
	};

	class func_F5 :public func_base
	{
	public:
		func_F5(size_t dims)
			:func_base(dims) {}
		void operator() (individual &ind);
		bool is_convergent(double fit,double vtr) { return (fit<=vtr || fit==optimum); }
	private:
		static const double optimum;
	};

	class func_Rosenbrock :public func_base
	{
	public:
		func_Rosenbrock(size_t dims)
			:func_base(dims) {}
		void operator() (individual &ind);
		bool is_convergent(double fit,double vtr) { return (fit<=vtr || fit==optimum); }
	private:
		static const double optimum;
	};

	class func_Step:public func_base
	{
	public:
		func_Step(size_t dims)
			:func_base(dims) {}
		void operator() (individual &ind);
		bool is_convergent(double fit,double vtr) { return (fit<=vtr || fit==optimum); }
	private:
		static const double optimum;
	};

	class func_Quartic_with_Noise:public func_base
	{
	public:
		func_Quartic_with_Noise(size_t dims)
			:func_base(dims) {}
		void operator() (individual &ind);
		bool is_convergent(double fit,double vtr) { return (fit<=vtr || fit==optimum); }
	private:
		static const double optimum;
	};

	class func_WS_Location:public func_base
	{
	public:
		func_WS_Location(size_t dims,std::string opt_coor_path="original_coor.txt");
		void operator() (individual &ind);
		bool is_convergent(double fit,double vtr) { return (fit<=vtr || fit==optimum); }
	private:
		static const double optimum;
		d_mat m_opt_coor;
	};

	class func_F2:public func_base
	{
	public:
		func_F2(size_t dims)
			:func_base(dims) {}
		void operator() (individual &ind);
		bool is_convergent(double fit,double vtr) { return (fit<=vtr || fit==optimum); }
	private:
		static const double optimum;
	};

	class func_F8:public func_base
	{
	public:
		func_F8(size_t dims)
			:func_base(dims) {}
		void operator() (individual &ind);
		bool is_convergent(double fit,double vtr) { return (fit<=vtr || fit==optimum); }
	private:
		static const double optimum;
	};

	class func_CamelBack:public func_base
	{
	public:
		func_CamelBack(size_t dims)
			:func_base(dims) {}
		void operator() (individual &ind);
		bool is_convergent(double fit,double vtr) { return (fit<=vtr || fit==optimum); }
	private:
		static const double optimum;
	};

	inline double u(double x,double a,double k,double m);

	class func_F12:public func_base
	{
	public:
		func_F12(size_t dims)
			:func_base(dims) {}
		void operator() (individual &ind);
		bool is_convergent(double fit,double vtr) { return (fit<=vtr || fit==optimum); }
	private:
		static const double optimum;
	}; 

	class func_F13:public func_base
	{
	public:
		func_F13(size_t dims)
			:func_base(dims) {}
		void operator() (individual &ind);
		bool is_convergent(double fit,double vtr) { return (fit<=vtr || fit==optimum); }
	private:
		static const double optimum;
	}; 

	class func_F3:public func_base
	{
	public:
		func_F3(size_t dims)
			:func_base(dims) {}
		void operator() (individual &ind);
		bool is_convergent(double fit,double vtr) { return (fit<=vtr || fit==optimum); }
	private:
		static const double optimum;
	};

	class func_F4:public func_base
	{
	public:
		func_F4(size_t dims)
			:func_base(dims) {}
		void operator() (individual &ind);
		bool is_convergent(double fit,double vtr) { return (fit<=vtr || fit==optimum); }
	private:
		static const double optimum;
	}; 

	class func_F14:public func_base
	{
	public:
		func_F14(size_t dims)
			:func_base(dims) {}
		void operator() (individual &ind);
		bool is_convergent(double fit,double vtr) { return (fit<=vtr || fit==optimum); }
	private:
		static const double optimum;
		static const double a[2][25];
	}; 

	class func_F15:public func_base
	{
	public:
		func_F15(size_t dims)
			:func_base(dims) {}
		void operator() (individual &ind);
		bool is_convergent(double fit,double vtr) { return (fit<=vtr || fit==optimum); }
	private:
		static const double optimum;
		static const boost::array<std::pair<double,double>,11> tab;
	}; 

	class func_Mod_Shekel:public func_base
	{
	public:
		func_Mod_Shekel(size_t dims)
			:func_base(dims) {}
		void operator() (individual &ind);
		bool is_convergent(double fit,double vtr) { return (fit<=vtr || fit==optimum); }
	protected:
		static const double m_a[30][10];
		static const double m_c[];
	private:
		static const double optimum;
	};

	class func_Mod_Langerman:public func_Mod_Shekel
	{
	public:
		func_Mod_Langerman(size_t dims)
			:func_Mod_Shekel(dims) {}
		void operator() (individual &ind);
		static const double optimum;
	};

	class func_Mod_Michalewicz:public func_base
	{
	public:
		func_Mod_Michalewicz(size_t dims)
			:func_base(dims) {}
		void operator() (individual &ind);
		bool is_convergent(double fit,double vtr) { return (fit<=vtr || fit==optimum); }
	private:
		static const double optimum;
		static const int m_m=10;
	};

	class func_Chebyshev:public func_base
	{
	public:
		func_Chebyshev(size_t dims)
			:func_base(dims) {}
		void operator() (individual &ind);
		bool is_convergent(double fit,double vtr) { return (fit<=vtr || fit==optimum); }
	private:
		static const double optimum;
	};

	class func_MOP:public std::unary_function<individual&,void>
	{
	public:
		func_MOP(size_t dims,size_t objs)
			:m_num_dims(dims){}
		virtual void operator() (individual &ind){}
		virtual ~func_MOP() {}
	protected:
		int m_num_dims;
	};

	class func_ZDT3:public func_base
	{
	public:
		func_ZDT3(size_t dims)
			:func_base(dims) {}
		void operator() (individual &ind);
		bool is_convergent(double fit,double vtr) {return false;}
	};
}// end namespace benchmark
#endif