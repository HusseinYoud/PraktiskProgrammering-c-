#include<cmath>
#include<numbers> // c++20
#include<limits>

namespace sfuns{
    constexpr double NaN = std::numeric_limits<double>::quiet_NaN();
    constexpr double pi = 3.14159265358979323846; // pre- c++20
// constexpr double pi = std::numbers::pi; // c++20
double fgamma(double x){
	if(x<0){
        return pi/std::sin(pi*x)/fgamma(1-x);
    }
    
    if(x<9){
        return fgamma(x+1)/x;
    }

	double lnfgamma=x*std::log(x+1/(12*x-1/x/10))-x + std::log(2*pi/x)/2;
	return std::exp(lnfgamma);
}

double lngamma(double x){
    if(x<0) return NaN;

    if(x<9) return lngamma(x+1) - std::log(x);

    double lnfgamma= x*std::log(x+1/(12*x-1/x/10))-x + std::log(2*pi/x)/2;
    return lnfgamma;
}


}

