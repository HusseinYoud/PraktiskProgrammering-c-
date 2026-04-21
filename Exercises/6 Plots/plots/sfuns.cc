#include<cmath>
#include<numbers> // c++20
#include"sfuns.h"
#include<limits>
#include<vector>

namespace sfuns{
constexpr double pi = std::numbers::pi;
constexpr double NaN = std::numeric_limits<double>::quiet_NaN();

double sgamma(double x){
if(x<0)return pi/std::sin(pi*x)/sgamma(1-x);
if(x<9)return sgamma(x+1)/x;
double lnsgamma=std::log(2*pi)/2+(x-0.5)*std::log(x)-x
    +(1.0/12)/x-(1.0/360)/(x*x*x)+(1.0/1260)/(x*x*x*x*x);
return std::exp(lnsgamma);
}



double erf(double x){
// single precision error function (Abramowitz and Stegun, from Wikipedia)
if(x<0) return -erf(-x);
std::vector<double> a {0.254829592,-0.284496736,1.421413741,-1.453152027,1.061405429};
double t=1/(1+0.3275911*x);
double sum=t*(a[0]+t*(a[1]+t*(a[2]+t*(a[3]+t*a[4]))));/* the right thing */
return 1-sum*std::exp(-x*x);
} 

double lngamma(double x){
if(x<=0) return NaN;
if(x<9) return lngamma(x+1)-std::log(x);
return x*std::log(x+1/(12*x-1/x/10))-x+std::log(2*pi/x)/2;
}
}