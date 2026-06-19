#pragma once
#include <iostream>
#include <vector>
#include <string>

//vector class taken from the vector exercise
struct vec {
    double x,y,z;

    vec(double x, double y, double z) : x(x), y(y), z(z) {}
    vec() : vec(0,0,0) {}
    vec(const vec&) = default;                                 // copy
    vec(vec&&) = default;                                      // move
    ~vec() = default;                                          // destructor

    // assignment
    vec& operator=(const vec&) = default;                     // copy assignment
    vec& operator=(vec&&) = default;                          // move assignment

    // arithmetic
    vec& operator+=(const vec&);
    vec& operator-=(const vec&);
    vec& operator*=(double);
    vec& operator/=(double);

    //Print statements, setter statement and others:
    void print(const std::string& s="") const;
    void set(double a, double b, double c) {x=a; y=b; z=c;}
    double norm() const;
};
// non-member operators
vec operator-(const vec&);
vec operator-(const vec&, const vec&);
vec operator+(const vec&, const vec&);
vec operator*(const vec&, double);
vec operator*(double, const vec&);
vec operator/(const vec&, double);

double dot(const vec&, const vec&);
vec cross(const vec&, const vec&);

// approximate equality
bool approx(double, double, double acc = 1e-6, double eps = 1e-6);
bool approx(const vec&, const vec&, double acc = 1e-6, double eps = 1e-6);
bool vec_approx(const vec& a, const vec& b, double acc = 1e-6, double eps = 1e-6);

std::ostream& operator<<(std::ostream&, const vec&);

struct matrix{
    std::vector<std::vector<double>> M;
    matrix(int rows, int cols) : M(rows, std::vector<double>(cols)) {}
    matrix(const matrix& other) = default;
    matrix(matrix&& other) = default;
    matrix& operator=(const matrix& other) = default;
    matrix& operator=(matrix&& other) = default;
    inline double& operator()(int i, int j) { return M[i][j]; }
    inline double& operator[](int i, int j) { return M[i][j]; }
    
    
}