#include<iostream>
#include<cmath>
#include<complex>
#include<limits>
#include "sfuns.h"
#include<iomanip>
void line(char c = '-'){
    std::cout <<std::string(30,c) << "\n";
}

int main(){
std::cout << "1." << "\n";
line();

float       f=1.0f; while((float)      (1.0f+f) != 1.0f){f/=2.0f;} f*=2.0f;
double      d=1.0d; while((double)     (1.0d+d) != 1.0d){d/=2.0d;} d*=2.0d;
long double l=1.0L; while((long double)(1.0L+l) != 1.0L){l/=2.0L;} l*=2.0L;
std::cout << "float eps=" << f << "   The system "
          << std::numeric_limits<float>::epsilon() << "\n" << " For comparison it should be " << std::pow(2,-23) << "\n";

std::cout << "double eps=" << d << "   The system "
          << std::numeric_limits<double>::epsilon() << "\n" << " For comparison it should be " << std::pow(2,-52) << "\n";

std::cout << "long double eps=" << l << "   The system "
          << std::numeric_limits<long double>::epsilon() << "\n";

line();
std::cout << "2." << "\n";
line();
std::cout << std::fixed << std::setprecision(20);

double tiny = d/2;
double a=1+tiny+tiny;
double b=tiny+tiny+1;
std::cout << "a==b ? " << (a==b ? "true":"false") << "\n";
std::cout << "a>1  ? " << (a>1  ? "true":"false") << "\n";
std::cout << "b>1  ? " << (b>1  ? "true":"false") << "\n";
std::cout << "I assume what is happening is that the compiler is rounding the tiny numbers to nearest number. We see that small numbers might hit a floating point underflow therefore a is not equal to b and why adding the two tiny numbers to 1 is different than adding 1 to the two tiny numbers" << "\n";

std::cout << "tiny=" << tiny << "\n";
std::cout << "1+tiny+tiny=" << a << "\n";
std::cout << "tiny+tiny+1=" << b << "\n";

line();
std::cout << "3." << "\n";
line();
double d1 = 0.1+0.1+0.1+0.1+0.1+0.1+0.1+0.1;
double d2 = 8*0.1;

std::cout << "d1==d2? " << (d1==d2 ? "true":"false") << "\n"; 

std::cout << "d1=" << d1 << "\n";
std::cout << "d2=" << d2 << "\n";

std::cout << "For the implementation of the approx function =" << (sfuns::approx(d1,d2) ? "true":"false") <<"\n"; 

return 0;
}
