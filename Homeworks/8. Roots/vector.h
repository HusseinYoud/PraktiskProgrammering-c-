#pragma once
#include <iostream>
#include <string>
#include <vector>
#include <cmath>
namespace pp{

//Matrix.h, Vector.h, vector.cc and matrix.cc are all copied from earlier homeworks. (ode too)

struct vec {
    std::vector<double> data;

    vec() = default;

    explicit vec(int n);
    vec(std::initializer_list<double> list);
    int size() const;
    
    double& operator[](int i);
    const double& operator[](int i) const;
    
    // arithmetic
    vec& operator+=(const vec&);
    vec& operator-=(const vec&);
    vec& operator*=(double);
    vec& operator/=(double);

    double norm() const;
    double dot(const vec& other) const;
    
    //Print statements, setter statement and others:
    void print(const std::string& s="") const;
};

// non-member operators
vec operator-(vec a, const vec& b);
vec operator-(vec a);
vec operator+(vec a, const vec& b);
vec operator*(vec a, double c);
vec operator*(double c, vec a);
vec operator/(vec a, double c);

double dot(const vec& a, const vec& b);
double norm(const vec& v);


// approximate equality
bool approx(double x, double y, double acc = 1e-9, double eps = 1e-9);
bool approx(const vec& a, const vec& b, double acc = 1e-9, double eps = 1e-9);

std::ostream& operator<<(std::ostream& os, const vec& v);
} // namespace pp