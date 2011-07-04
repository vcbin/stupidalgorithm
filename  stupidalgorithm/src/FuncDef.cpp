#define _USE_MATH_DEFINES // for math constant eg. pi,e... 
#include <cmath>

#include "FuncDef.h"
#include "Rand_Val.h"
#include "Algo_DataStruct.h"

#include <fstream>
#include <sstream>
#include <string>

#include <boost/array.hpp>

using std::vector;
using std::pair;
using std::make_pair;
using boost::array;

using boost::uniform_01;
using boost::normal_distribution;
//using boost::cauchy_distribution;
using boost::mt19937;
using boost::variate_generator;

using std::string;
using std::vector;
using std::ifstream;
using std::stringstream;
using std::logic_error;

extern boost::mt19937 gen;

namespace benchmark
{
	// global optimum of test functions
	const double func_Sphere::optimum=0.0;
	const double func_F5::optimum=-1.0;
	const double func_Griewank::optimum=0.0;
	const double func_AckeyF1::optimum=0.0;
	const double func_Rastrigin::optimum=0.0;
	const double func_Rosenbrock::optimum=0.0;
	const double func_Step::optimum=0.0;
	const double func_Quartic_with_Noise::optimum=0.0;
	const double func_WS_Location::optimum=0.0;
	const double func_F2::optimum=0.0;
	const double func_F8::optimum=-12569.5;
	const double func_CamelBack::optimum=-1.0316285;
	const double func_F12::optimum=0;
	const double func_F13::optimum=0;
	const double func_F3::optimum=0;
	const double func_F4::optimum=0;
	const double func_F14::optimum=1;
	const double func_F15::optimum=0.0003075;

	// No single optimum on different dimension
	const double func_Mod_Shekel::optimum=MIN;
	const double func_Mod_Langerman::optimum=MIN;
	const double func_Mod_Michalewicz::optimum=MIN;
	const double func_Chebyshev::optimum=MIN;

	const double func_Mod_Shekel::m_a[30][10] = 
	{
		{9.681, 0.667, 4.783, 9.095, 3.517, 9.325, 6.544, 0.211, 5.122, 2.020},
		{9.400, 2.041, 3.788, 7.931, 2.882, 2.672, 3.568, 1.284, 7.033, 7.374},
		{8.025, 9.152, 5.114, 7.621, 4.564, 4.711, 2.996, 6.126, 0.734, 4.982},
		{2.196, 0.415, 5.649, 6.979, 9.510, 9.166, 6.304, 6.054, 9.377, 1.426},
		{8.074, 8.777, 3.467, 1.863, 6.708, 6.349, 4.534, 0.276, 7.633, 1.567},
		{7.650, 5.658, 0.720, 2.764, 3.278, 5.283, 7.474, 6.274, 1.409, 8.208},
		{1.256, 3.605, 8.623, 6.905, 0.584, 8.133, 6.071, 6.888, 4.187, 5.448},
		{8.314, 2.261, 4.224, 1.781, 4.124, 0.932, 8.129, 8.658, 1.208, 5.762},
		{0.226, 8.858, 1.420, 0.945, 1.622, 4.698, 6.228, 9.096, 0.972, 7.637},
		{7.305, 2.228, 1.242, 5.928, 9.133, 1.826, 4.060, 5.204, 8.713, 8.247},
		{0.652, 7.027, 0.508, 4.876, 8.807, 4.632, 5.808, 6.937, 3.291, 7.016},
		{2.699, 3.516, 5.874, 4.119, 4.461, 7.496, 8.817, 0.690, 6.593, 9.789},
		{8.327, 3.897, 2.017, 9.570, 9.825, 1.150, 1.395, 3.885, 6.354, 0.109},
		{2.132, 7.006, 7.136, 2.641, 1.882, 5.943, 7.273, 7.691, 2.880, 0.564},
		{4.707, 5.579, 4.080, 0.581, 9.698, 8.542, 8.077, 8.515, 9.231, 4.670},
		{8.304, 7.559, 8.567, 0.322, 7.128, 8.392, 1.472, 8.524, 2.277, 7.826},
		{8.632, 4.409, 4.832, 5.768, 7.050, 6.715, 1.711, 4.323, 4.405, 4.591},
		{4.887, 9.112, 0.170, 8.967, 9.693, 9.867, 7.508, 7.770, 8.382, 6.740},
		{2.440, 6.686, 4.299, 1.007, 7.008, 1.427, 9.398, 8.480, 9.950, 1.675},
		{6.306, 8.583, 6.084, 1.138, 4.350, 3.134, 7.853, 6.061, 7.457, 2.258},
		{0.652, 2.343, 1.370, 0.821, 1.310, 1.063, 0.689, 8.819, 8.833, 9.070},
		{5.558, 1.272, 5.756, 9.857, 2.279, 2.764, 1.284, 1.677, 1.244, 1.234},
		{3.352, 7.549, 9.817, 9.437, 8.687, 4.167, 2.570, 6.540, 0.228, 0.027},
		{8.798, 0.880, 2.370, 0.168, 1.701, 3.680, 1.231, 2.390, 2.499, 0.064},
		{1.460, 8.057, 1.336, 7.217, 7.914, 3.615, 9.981, 9.198, 5.292, 1.224},
		{0.432, 8.645, 8.774, 0.249, 8.081, 7.461, 4.416, 0.652, 4.002, 4.644},
		{0.679, 2.800, 5.523, 3.049, 2.968, 7.225, 6.730, 4.199, 9.614, 9.229},
		{4.263, 1.074, 7.286, 5.599, 8.291, 5.200, 9.214, 8.272, 4.398, 4.506},
		{9.496, 4.830, 3.150, 8.270, 5.079, 1.231, 5.731, 9.494, 1.883, 9.732},
		{4.138, 2.562, 2.532, 9.661, 5.611, 5.500, 6.886, 2.341, 9.699, 6.500}
	};

	const double func_Mod_Shekel::m_c[] = 
	{
		0.806,
		0.517,
		0.1,
		0.908,
		0.965,
		0.669,
		0.524,
		0.902,
		0.531,
		0.876,
		0.462,
		0.491,
		0.463,
		0.714,
		0.352,
		0.869,
		0.813,
		0.811,
		0.828,
		0.964,
		0.789,
		0.360,
		0.369,
		0.992,
		0.332,
		0.817,
		0.632,
		0.883,
		0.608,
		0.326
	};

	void func_Sphere::operator() (individual &ind)
	{
		int i;
		ind.obj[0]=0.0;
		for ( i=0;i<m_num_dims;i++ )
			ind.obj[0] += (ind.x[i]*ind.x[i]);
	}

	void func_F5::operator() (individual &ind)
	{
		int i;
		ind.obj[0]=0.0;
		double quad_sum=0.0;
		for ( i=0;i<m_num_dims;i++ )
		{
			quad_sum += (ind.x[i]*ind.x[i]);
		}
		double sin_tmp=sin( sqrt(quad_sum) );
		sin_tmp=(sin_tmp*sin_tmp);
		double numerator,denominator;
		numerator=(sin_tmp-0.5);
		double den_tmp=1+0.001*quad_sum;
		denominator=den_tmp*den_tmp;
		ind.obj[0]=  numerator/denominator - 0.5;
	}

	void func_Griewank::operator() (individual &ind)
	{
		int i;
		static const double cof=1.0/4000.0;
		ind.obj[0]=0.0;
		double tmp=1.0;
		for ( i=0;i<m_num_dims;i++ )
		{
			ind.obj[0] += (ind.x[i]*ind.x[i]);
			tmp *= cos( ind.x[i]/sqrt((i+1)*1.0) );
		}
		ind.obj[0] *= cof;
		ind.obj[0] -= tmp;
		ind.obj[0] += 1.0;
	}

	void func_AckeyF1::operator() (individual &ind)
	{
		double quad_sum=0.0;
		double cos_sum=0.0;
		int i;
		ind.obj[0]=0.0;
		for ( i=0;i<m_num_dims;i++ )
		{
			quad_sum += (ind.x[i]*ind.x[i]);
			cos_sum += cos(2*M_PI*ind.x[i]);
		}
		quad_sum /= m_num_dims;
		cos_sum /= m_num_dims;
		quad_sum=sqrt(quad_sum);
		ind.obj[0] = M_E + 20 - 20*exp(-0.2*quad_sum) - exp(cos_sum);
	}

	void func_Rastrigin::operator() (individual &ind)
	{
		int i;
		ind.obj[0]=0.0;
		for ( i=0;i<m_num_dims;i++ )
			ind.obj[0] += ( ind.x[i]*ind.x[i] + 10 - 10*cos(2*M_PI*ind.x[i]) );
	}

	void func_Rosenbrock::operator() (individual &ind)
	{
		int i;
		ind.obj[0]=0.0;
		for ( i=0;i<m_num_dims-1;i++ )
		{
			double tmp1,tmp2;
			tmp1=ind.x[i+1]-(ind.x[i]*ind.x[i]);
			tmp2=ind.x[i]-1.0;
			ind.obj[0] += (100*tmp1*tmp1 + tmp2*tmp2);
		}
	}

	void func_Step::operator() (individual &ind)
	{
		int i;
		ind.obj[0]=0.0;
		for ( i=0;i<m_num_dims;i++ )
		{
			double tmp=floor(ind.x[i]+0.5);
			ind.obj[0] += (tmp*tmp);
		}
	}

	void func_Quartic_with_Noise::operator() (individual &ind)
	{
		int i;
		ind.obj[0]=0.0;
		for ( i=0;i<m_num_dims;i++ )
		{
			double tmp=ind.x[i]*ind.x[i];
			tmp=tmp*tmp;
			ind.obj[0] += ((i+1)*tmp);
		}
		uniform_01<> uni_01;
		variate_generator<mt19937&, uniform_01<> > rnd_01(gen, uni_01);
		ind.obj[0] += rnd_01();
	}

	func_WS_Location::func_WS_Location(size_t dims,string opt_coor_path)
		:func_base(dims)
	{
		ifstream opt_coor_data(opt_coor_path.c_str());
		int i=0;// entry read number
		int j;
		string inbuf;
		int col_num=2;
		d_array x_tmp(col_num);
		if ( opt_coor_data.fail() )
		{
			stringstream str_err;
			str_err<<"\nopen file "
				<<opt_coor_path
				<<" ERROR!\n"
				<<"please check file path!\n";
			throw logic_error(str_err.str());
		}
		while ( !opt_coor_data.eof() )
		{
			getline(opt_coor_data,inbuf);
			if ( false==is_dataLine(inbuf) )
				continue;

			stringstream line(inbuf);
			for ( j=0;j<col_num;j++ )
			{
				line>>x_tmp[j];
			}
			m_opt_coor.resize(m_opt_coor.size()+1,d_array(col_num));
			m_opt_coor[i]=x_tmp;
			if ( opt_coor_data.fail() ) // read error
			{
				stringstream str_err("ERROR:Read data from file");
				str_err<<" "
					<<opt_coor_path.c_str()
					<<" "
					<<"failed."
					<<"Check your data file for data validity."
					<<"\n";
				throw logic_error(str_err.str());
			}
			i++;
		}// while !eof
	}

	void func_WS_Location::operator() (individual &ind)
	{
		int opt_coor_num=m_opt_coor.size();
		int node_num=m_num_dims/2;
		if ( opt_coor_num!=node_num )
		{
			stringstream str_err("ERROR:data inconsistant!");
			str_err<<" "
				<<" "
				<<"failed."
				<<"Check your data file for data validity."
				<<"\n";
			throw logic_error(str_err.str());
		}
		d_array coor_tmp(2),opt_tmp(2);
		int i;
		int opt_coor_idx=0;
		ind.obj[0]=0.0;

		for ( i=0;i<m_num_dims-1;i+=2 )
		{
			coor_tmp[0]=ind.x[i];// x coordinate
			coor_tmp[1]=ind.x[i+1];// y coordinate
			ind.obj[0] += (vec_distance(coor_tmp,m_opt_coor[opt_coor_idx]));
			opt_coor_idx++;
		}
		normal_distribution<> std_norm_dist;
		variate_generator<mt19937&, normal_distribution<> > std_norm(gen, std_norm_dist);
		ind.obj[0] /= node_num;
		ind.obj[0] += std_norm();
	}

	void func_F2::operator() (individual &ind)
	{
		int i;
		ind.obj[0]=0.0;
		double tmp1,tmp2;
		tmp1=0.0;
		tmp2=1.0;
		for ( i=0;i<m_num_dims;i++ )
		{
			tmp1 += abs(ind.x[i]);
			tmp2 *= abs(ind.x[i]);
		}
		ind.obj[0] = (tmp1+tmp2);
	}

	void func_F8::operator() (individual &ind)
	{
		int i;
		ind.obj[0]=0.0;
		for ( i=0;i<m_num_dims;i++ )
			ind.obj[0] += ( ind.x[i]*sin(sqrt(abs(ind.x[i]))) );
		ind.obj[0] *= -1.0;
	}

	void func_CamelBack::operator() (individual &ind)
	{
		ind.obj[0] = 4*ind.x[0]*ind.x[0]-2.1*pow(ind.x[0],4)+1/3.0*pow(ind.x[0],6)+
			ind.x[0]*ind.x[1]-
			4*ind.x[1]*ind.x[1]+4*pow(ind.x[1],4);
	}

	// u function for F12,F13
	inline double u(double x,double a,double k,double m)
	{
		if ( x>a )
			return k*pow((x-a),m);
		else if ( abs(x)<=a )
			return 0;
		else //  x<-a
			return k*pow(-x-a,m);
	}

	void func_F12::operator() (individual &ind)
	{
		int i;
		double sum_x=0.0,y_i,y_next;
		double sum_u=0.0;
		double cof=1/4.0;
		double y1,yn;
		for ( i=0;i<m_num_dims;i++ )
		{
			y_i=1+cof*(ind.x[i]+1);
			if ( 0==i )
				y1=y_i;
			if ( i==m_num_dims-1 )
				yn=y_i;
			if (i<m_num_dims-1)
			{
				y_next=1+cof*(ind.x[i+1]+1);
				double tmp_y_i=y_i-1;
				double tmp_sin=sin(M_PI*y_next);
				sum_x += ( tmp_y_i*tmp_y_i*(1+10*tmp_sin*tmp_sin) );
			}
			sum_u += u(ind.x[i],10,100,4);
		}
		double tmp_sin=sin(M_PI*y1);
		double tmp_yn=yn-1;
		ind.obj[0] = M_PI/m_num_dims*(10*tmp_sin*tmp_sin+sum_x+tmp_yn*tmp_yn)+
			sum_u;
	}

	void func_F13::operator() (individual &ind)
	{
		int i;
		double sum_x=0.0;
		double sum_u=0.0;
		int max_idx=m_num_dims-1;
		for ( i=0;i<m_num_dims;i++ )
		{
			if (i<max_idx)
			{
				double tmp_x_i=ind.x[i]-1;
				double tmp_sin=sin(3*M_PI*ind.x[i+1]);
				sum_x += ( tmp_x_i*tmp_x_i*(1+tmp_sin*tmp_sin) );
			}
			sum_u += u(ind.x[i],5,100,4);
		}
		double tmp_sin1=sin(M_PI*3*ind.x[0]);
		double tmp_xn=ind.x[max_idx]-1;
		double tmp_sin2=sin(2*M_PI*ind.x[max_idx]);
		ind.obj[0] = 0.1*( tmp_sin1*tmp_sin1+sum_x+tmp_xn*tmp_xn*(1+tmp_sin2*tmp_sin2)
			) + sum_u;
	}

	void func_F3::operator() (individual &ind)
	{
		int i,j;
		ind.obj[0]=0.0;
		for ( i=0;i<m_num_dims;i++ )
		{
			double sum_x=0.0;
			for ( j=0;j<i;j++ )
				sum_x += ind.x[j];
			sum_x=sum_x*sum_x;
			ind.obj[0] += sum_x;
		}
	}

	void func_F4::operator() (individual &ind)
	{
		int i;
		ind.obj[0]=abs(ind.x[0]);
		double tmp;
		for ( i=1;i<m_num_dims;i++ )
		{
			tmp=abs(ind.x[i]);
			if ( tmp>ind.obj[0] )
				ind.obj[0]=tmp;
		}
	}

	const double func_F14::a[2][25] = 
	{
		{-32,-16,0,16,32,-32,-16,0,16,32,-32,-16,0,16,32,-32,-16,0,16,32,-32,-16,0,16,32},
		{-32,-32,-32,-32,-32,-16,-16,-16,-16,-16,0,0,0,0,0,16,16,16,16,16,32,32,32,32,32}
	};

	void func_F14::operator() (individual &ind)
	{
		int i, j;
		double sum,x_sum;

		ind.obj[0]=0.0;
		sum = 0.0;
		for (j = 0; j < 25; j++) 
		{
			x_sum=0.0;
			for (i = 0; i < 2; i++) 
				x_sum += pow(ind.x[i] - a[i][j],6);
			x_sum += (j+1);
			x_sum=1.0/x_sum;
			sum += x_sum;
		}
		ind.obj[0]=(1.0/500+sum);
		ind.obj[0]=(1.0/ind.obj[0]);
	}

	const array<pair<double,double>,11> func_F15::tab=
	{ make_pair(0.1957, 0.25),
		make_pair(0.1947, 0.5), make_pair(0.1735, 1), make_pair(0.16, 2),
		make_pair(0.0844, 4), make_pair(0.0627, 6), make_pair(0.0456, 8),
		make_pair(0.0342, 10), make_pair(0.0323, 12), make_pair(0.0235, 14),
		make_pair(0.0246, 16) };

	void func_F15::operator() (individual &ind)
	{
		int i;
		double bi,bi_sqr,tmp;

		ind.obj[0]=0.0;
		for (i = 0; i < 11; i++) 
		{
			bi=1.0/tab[i].second;
			bi_sqr=bi*bi;
			tmp=( tab[i].first-(ind.x[0]*(bi_sqr+bi*ind.x[1]))/(bi_sqr+bi*ind.x[2]+ind.x[3]) );
			tmp=tmp*tmp;
			ind.obj[0] += tmp;
		}
	}

	void func_Mod_Shekel::operator() (individual &ind)
	{		
		int i, j;
		double sp, h, result = 0.0;

		for (i = 0; i < 30; i++) 
		{
			sp = 0.0;
			for (j = 0; j < m_num_dims; j++) 
			{
				h = ind.x[j] - m_a[i][j];
				sp += h * h;
			}
			result += 1.0 / (sp + m_c[i]);
		}
		ind.obj[0]=-1.0*result;
	}

	void func_Mod_Langerman::operator() (individual &ind)
	{		
		int i,j;
		double sum,d,dist,temp1,temp2,temp20,temp21;
		sum = 0.0;
		for ( i = 0; i < 5; i++ )
		{	
			dist = 0.0;
			for ( j= 0; j<m_num_dims; j++ )
			{
				d =ind.x[j] - m_a[i][j];
				temp1=(d*d);
				dist =dist + temp1;
				//printf("%1f*%1f|",x[j],m_a[i][j]);
			}

			//dist = SqrDst(x, a[i], nn);
			temp20=exp(-dist/M_PI);
			temp21=cos( M_PI * dist ) ;
			temp2=m_c[i] * (temp20*temp21);
			//printf("\n dist=%1f ++ %1f,%1f,%1f \n",dist,temp20,temp21,temp2);

			sum -= temp2;
			//printf("\nSum=*%1f**",-sum);
		}
		//printf("\n E sum=*%1f**",sum);
		ind.obj[0]=(-(double)sum);
	}

	void func_Mod_Michalewicz::operator() (individual &ind)
	{
		int i;
		double sum;
		double y;
		double tmp;

		sum=0.0;
		for ( i=0;i<m_num_dims;i++ )
		{
			if ( i!=m_num_dims-1 )
			{
				if ( 1==(i+1)%2 )
					y=ind.x[i]*cos(M_PI/6)-ind.x[i+1]*sin(M_PI/6);
				else
					y=ind.x[i]*cos(M_PI/6)+ind.x[i+1]*sin(M_PI/6);
			}
			else
				y=ind.x[i];

			tmp=sin((i+1)*y*y/M_PI);
			tmp=pow(tmp,2*m_m);
			sum += sin(y)*tmp;
		}
		ind.obj[0]=-1.0*sum;
	}

	void func_Chebyshev::operator() (individual &ind)  /*  Chebychev Polynomial */
	{                                                    /* Valid for D=9 or D=17 */
		int i,j;
		double px,x=-1,dx;

		if ( m_num_dims == 9 )
			dx=72.66066;
		if ( m_num_dims == 17 )
			dx=10558.145;

		ind.obj[0]=0.0;
		for (i=0;i<=100;i++)
		{
			px=ind.x[0];
			for (j=1;j<m_num_dims;j++) px=x*px+ind.x[j];
			if (px<-1 || px>1) ind.obj[0]+=(1.-px)*(1.-px);
			x+=.02;
		}
		px=ind.x[0];
		for (j=1;j<m_num_dims;j++) px=1.2*px+ind.x[j];
		px=px-dx;
		if (px<0) ind.obj[0]+=px*px;
		px=ind.x[0];
		for (j=1;j<m_num_dims;j++) px=-1.2*px+ind.x[j];
		px=px-dx;
		if (px<0) ind.obj[0]+=px*px;
	}

	void func_ZDT3::operator() (individual &ind)
	{
		double f1, f2, g, h;
		int i;
		f1 = ind.x[0];
		g = 0.0;
		for (i=1; i<30; i++)
			g += ind.x[i];
		g = 9.0*g/29.0;
		g += 1.0;
		h = 1.0 - sqrt(f1/g) - (f1/g)*sin(10.0*M_PI*f1);
		f2 = g*h;
		ind.obj[0] = f1;
		ind.obj[1] = f2;
		return;
	}

}// end namespace benchmark
