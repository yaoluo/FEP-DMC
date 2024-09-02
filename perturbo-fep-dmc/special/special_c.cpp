#include<boost/math/special_functions/erf.hpp>
#include<math.h>
using namespace boost;
using namespace std;
extern "C" void erf_c_(double &x, double&y)
{
	const double invsq2 = 1.0/sqrt(2.0);
	double z;
	z = x*invsq2;
	y = 0.5*math::erf(z) + 0.5; 
	
};

extern "C" void inverf_c_(double &x, double&y)
{
 	const double sq2 = sqrt(2.0);
	double z;	
        z = 2.0*(x-0.5);	
	
	y = sq2*math::erf_inv(z);

};





