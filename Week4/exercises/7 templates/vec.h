#pragma once
#include <iostream>
#include <string>
#include<complex>
#include<cmath>


template<typename T>
struct vec {
    T x, y, z;

    // constructors
    vec(T x, T y, T z) : x(x), y(y), z(z) {}         // templated constructor
    vec() : vec(T(0), T(0), T(0)) {}                                    // default
    vec(const vec&) = default;                                 // copy
    vec(vec&&) = default;                                      // move
    ~vec() = default;                                          // destructor

    // assignment
    vec& operator=(const vec&) = default;                     // copy assignment
    vec& operator=(vec&&) = default;                          // move assignment

    // arithmetic
    vec& operator+=(const vec& other){
        x += other.x; y += other.y; z += other.z;
        return *this;
    }
    vec& operator-=(const vec& other){
        x -= other.x; y -= other.y; z -= other.z;
        return *this;
    }
    vec& operator*=(T c){
        x *= c; y *= c; z *= c;
        return *this;
    }
    vec& operator/=(T c){
        x /= c; y /= c; z /= c;
        return *this;
    }

    //Print statements, setter statement and others:
    void print(const std::string& s="") const {
        std::cout << s << x << " " << y << " " << z << "\n";
    }
    void set(T a, T b, T c) {x=a; y=b; z=c;}
    T norm() const {
        return std::sqrt(x*x + y*y + z*z);
    }
};

// non-member operators
template<typename T>
vec<T> operator-(const vec<T>& a){
    return vec<T>(-a.x, -a.y, -a.z);
}
template<typename T>
vec<T> operator-(const vec<T>& a, const vec<T>&b) {
    return vec<T>(a.x - b.x, a.y - b.y, a.z - b.z);
}


template<typename T>
vec<T> operator+(const vec<T>& a, const vec<T>& b) {
    return vec<T>(a.x + b.x, a.y + b.y, a.z + b.z);
}


template<typename T>
vec<T> operator*(const vec<T>& a, T c){
    return vec<T>(a.x * c, a.y * c, a.z * c);
}


template<typename T>
vec<T> operator*(T c, const vec<T>& a){
    return a*c;
}

template<typename T>
vec<T> operator/(const vec<T>& a, T c){
    return vec<T>(a.x / c, a.y / c, a.z / c);
}

template<typename T>
T dot(const vec<T>& a, const vec<T>& b){
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

template<typename T>
vec<T> cross(const vec<T>& a, const vec<T>& b){
    return vec<T>(
        a.y * b.z - a.z * b.y, 
        a.z * b.x - a.x * b.z, 
        a.x * b.y - a.y * b.x
    );
}



// approximate equality
template<typename T>
bool approx(T a, T b, T acc = 1e-6, T eps = 1e-6){
    return std::abs(a - b) <= acc || std::abs(a - b) <= eps * std::max(std::abs(a), std::abs(b));
}


template<typename T>
bool approx(const vec<T>& a, const vec<T>& b, T acc = 1e-6, T eps = 1e-6){
    return approx(a.x, b.x, acc, eps) && approx(a.y, b.y, acc, eps) && approx(a.z, b.z, acc, eps);
}


template<typename T>
bool vec_approx(const vec<T>& a, const vec<T>& b, T acc = 1e-6, T eps = 1e-6){
    return approx(a,b,acc,eps);
}

template<typename T>
std::ostream& operator<<(std::ostream& out, const vec<T>& a){
    out << a.x << " " << a.y << " " << a.z;
    return out;
}