#include<iostream>
#include<cmath>
#include "vec.h"
#include<limits>
#include<iomanip>

void line(char c = '-'){
    std::cout <<std::string(10,c) << "\n";
}
int main(){

	vec a(1, 2, 3);
	vec b(4, 5, 6);

	a.print("a: ");
	b.print("b: ");
	line();
	std::cout <<"Vectors done" <<"\n";
	line();
	std::cout <<"We'd now like to test if we can change the vector and use all of our operators" <<"\n";
	line();
	a.set(2,1,5);
	a.print("New set a is = ");

	b=a;
	b.print("b=a, b is now = ");

	vec c = b+a;
	c.print("c = b + a, c is now = ");

	vec d = b-a;
	d.print("d = b - a, d is now = ");

	vec e = 5*a;
	e.print("e = 5*a, e is now = ");

	vec f = a/2;
	f.print("f = a/2, f is now = ");


	double an = a.norm();
	std::cout << "norm of a is = " << an << "\n";

	double adotb = dot(a,b);
	std::cout << "dot product of a and b is = " << adotb << "\n";

	vec acrossb = cross(a,b);
	acrossb.print("cross product of a and b is = ");
	
	return 0;
	}

