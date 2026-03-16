#include<cmath>
#include<complex>
#include<iostream>
#include"sfuns.h"
#include<iostream>
#include<iomanip>
using complex=std::complex<double>;
constexpr double  pi = 3.14159265358979324;
constexpr double  ee = 2.71828182845904523;
constexpr complex ii = complex(0,1);
void line(char c = '-'){
    std::cout <<std::string(10,c) << "\n";
}
int main(){
    
    std::cout << std::setprecision(6) << std::fixed;

    std::cout << "1." << "\n";
    line();   
    std::cout << "   log(ii)=" << std::log(ii)  <<"\n";
    std::cout << "   ii^ii=" << std::pow(ii,ii) <<"\n";
    std::cout << "   pi^ii=" << std::pow(pi,ii) <<"\n";
    std::cout << "   ee^ii=" << std::pow(ee,ii) <<"\n";
    std::cout << "   pi^ee=" << std::pow(pi,ee) <<"\n";
    std::cout << "   ee^ii=" << std::pow(ee,ii) <<"\n";
    std::cout << "   ee^pi=" << std::pow(ee,pi) <<"\n";
    std::cout << "   2^(1/5)= " << std::pow(2.0,1.0/5.0) << "\n";
    std::cout << "   sqrt(2)=" << std::sqrt(2.0) << "\n";
    line();
    std::cout << "2." << "\n";
    line();
    
    for(int i=1; i<=10; i++){
        std::cout << "   fgamma(" << i << ")=" << sfuns::fgamma(i) << "\n";
    }
    line();
    std::cout << "3." << "\n";
    line();
    for(int i=1; i<=10; i++){
        std::cout << "   lngamma(" << i << ")=" << sfuns::lngamma(i) << "\n";
    }
    return 0;
}
