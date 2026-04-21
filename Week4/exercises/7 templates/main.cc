#include<iostream>
#include<cmath>
#include "vec.h"
#include<limits>
#include<iomanip>

void line(char c = '-'){
    std::cout <<std::string(10,c) << "\n";
}
int main(){

	vec<double> a(1.1, 2.3, 3.4);
	vec<double> b(4.6, 5.2, 6.5);

	a.print("a: ");
	b.print("b: ");
	line();
    std::cout << "a + b = " << a + b << "\n";
    std::cout << "a - b = " << a - b << "\n";
    std::cout << "dot(a,b) = " << dot(a,b) << "\n";
    std::cout << "cross(a,b) = " << cross(a,b) << "\n";
    std::cout << "norm(a) = " << a.norm() << "\n";

    vec<float> c(1.11f, 2.21f, 3.3f);
    c.print("c = ");

	return 0;
	}

