#include"vector.h"
#include<algorithm>
#include<cmath>
#include<iostream>
#include<string>
#include<vector>
#include<algorithm>
#include <cassert>
#include <iomanip>

//entire vector.cc file is taken from the homework 2
namespace pp{
    vec::vec(int n) : data(n) {}
    vec::vec(std::initializer_list<double> list) : data(list) {}
    int vec::size() const { return static_cast<int>(data.size()); }

    double& vec::operator[](int i) { return data[i]; }
    const double& vec::operator[](int i) const { return data[i]; }
    

    double vec::dot(const vec& other) const{
        assert(size() == other.size());
        double result = 0.0;
        for (int i = 0; i < size(); ++i){
            result += (*this)[i] * other[i];
        }
        return result;
    }
    
    double vec::norm() const{
        return std:: sqrt(dot(*this));
    }


    // arithmetic
    vec& vec::operator+=(const vec& other){
        assert(size() == other.size());
        for (int i = 0; i < size(); ++i) (*this)[i] += other[i];
        return *this;

    }

    vec& vec::operator-=(const vec& other){
        assert(size() == other.size());
        for (int i = 0; i < size(); ++i) (*this)[i] -= other[i];
        return *this;
    }

    vec& vec::operator*=(double c){
        for (int i = 0; i < size(); ++i) (*this)[i] *= c;
        return *this;
    }

    vec& vec::operator/=(double c) {
        for (int i=0; i<size(); ++i) (*this)[i] /= c;
        return *this;
    }

    //Print statements, setter statement and others:
    void vec::print(const std::string& s) const {
        if (!s.empty()) {
            std::cout << s << '\n';
        }

        for (double x : data) {
            std::cout << std::setw(12) << x << ' ';
        }

        std::cout << '\n';
    }


    

// non-member operators
vec operator-(vec a, const vec& b) { a -= b; return a; }
vec operator-(vec a) { a*=-1.0; return a; }
vec operator+(vec a, const vec& b) { a += b; return a; }
vec operator*(vec a, double c) { a *= c; return a; }
vec operator*(double c, vec  a) { a *= c; return a; }
vec operator/(vec a, double c) { a /= c; return a; }

std::ostream& operator <<(std::ostream& os, const vec& v) {
    os << "[";
    for (int i = 0; i<v.size(); ++i) {
        os << v[i];
        if (i < v.size() - 1) os << ", ";
    }
    os << "]";
    return os;
}

double dot(const vec& a, const vec& b) {
 return a.dot(b);
}

double norm(const vec& v) {
    return v.norm();
}

// approximate equality
bool approx(double x, double y, double acc, double eps){
    double d = std::abs(x - y);
    if (d <= acc) return true;
    return d <= eps * std::max(std::abs(x), std::abs(y));
}
bool approx(const vec& a, const vec& b, double acc, double eps){
    if (a.size() != b.size()) return false;
    for (int i = 0; i < a.size(); ++i){
        if (!approx(a[i], b[i], acc, eps)) return false;
    }
    return true;
}
} // namespace pp
