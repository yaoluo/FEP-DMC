#include<boost/math/special_functions/erf.hpp>
#include<iostream>
using namespace std;
using namespace boost;

int main(){
        double x[100],y[100];
	double a = 0.9;
        int i;
	x[0]=1.0;
	for(i=1;i<99;i++){
		x[i] = a*x[i-1];
        	y[i] = math::erf_inv(1.0-x[i]);
		cout<<1-x[i]<<", "<<y[i]<<"\n";
	}
	return 0; 
}; 
