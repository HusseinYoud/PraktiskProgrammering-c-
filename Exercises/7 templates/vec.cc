#include"vec.h"
#include<algorithm>
#include<cmath>
#include<iostream>
#include<string>

// void vec::print(const std::string& s) const{
// 	std::cout << s << x << " " << y << " " << z << "\n";
// }


// //operations between vectors:
// vec& vec::operator+=(const vec& v){
//     x += v.x; y += v.y; z += v.z;
//     return *this;
// }

// vec& vec::operator-=(const vec& v){
//     x -= v.x; y -= v.y; z -= v.z;
//     return *this;
// }

// vec& vec::operator*=(double s){
//     x *= s; y *= s; z *= s;
//     return *this;
// }

// vec& vec::operator/=(double s){
//     x /= s; y /= s; z /= s;
//     return *this;
// }

// //vector operations def:
// vec operator+(const vec& a, const vec& b){
//     return vec(a.x + b.x, a.y + b.y, a.z + b.z);
// }
// vec operator-(const vec& a, const vec& b){
//     return vec(a.x - b.x, a.y - b.y, a.z - b.z);
// }

// vec operator*(const vec& v, double s){
//     return vec(v.x * s, v.y * s, v.z * s);
// }

// vec operator/(const vec& v, double s){
//     return vec(v.x / s, v.y / s, v.z / s);
// }

// vec operator-(const vec& v){
//     return vec(-v.x, -v.y, -v.z);
// }

// vec operator*(double s, const vec& v){
//     return v * s;
// }

// double vec::norm() const{
//     double normsq = x*x + y*y + z*z;
//     return std::pow(normsq, 0.5);
// }
// double dot(const vec& a, const vec& b){
//     double scalar = a.x * b.x + a.y * b.y + a.z * b.z;
//     return scalar;
// }

// vec cross(const vec& a, const vec& b){
//     double nx = a.y * b.z - a.z * b.y;
//     double ny = a.z * b.x - a.x * b.z;
//     double nz = a.x * b.y - a.y * b.x;
//     return vec(nx, ny, nz);
// }

// bool approx(double a, double b, double acc, double eps){
//     double diff = std::abs(a - b);
//     if (diff <= acc) return true;
//     double max_ab = std::max(std::abs(a), std::abs(b));
//     return diff <= eps*max_ab;
// }

// bool approx(const vec& a, const vec& b, double acc, double eps){
//     double diff = (a - b).norm();
//     if (diff <= acc) return true;
//     double max_ab = std::max(a.norm(), b.norm());
//     return diff <= eps*max_ab;
// } 


// bool vec_approx(const vec& a, const vec& b, double acc, double eps){
// 	if(!approx(a.x,b.x, acc, eps))return false;
// 	if(!approx(a.y,b.y, acc, eps))return false;
// 	if(!approx(a.z,b.z, acc, eps))return false;
// 	return true;
// }

// std::ostream& operator<<(std::ostream& os, const vec& v){
//     os << "{ " << v.x << ", " << v.y << ", " << v.z << " } ";
//     return os;
//     }