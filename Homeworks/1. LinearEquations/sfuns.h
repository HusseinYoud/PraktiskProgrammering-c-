#pragma once
#include<iostream>
#include<cstdio>
#include<iomanip>
#include<cmath>
#include<string>
#include<vector>
#include<initializer_list>
#include<functional>
#include<ranges>
#include<utility>
#include<cassert>
#include<random>
#include<stdexcept>

namespace pp{
//vector class taken from the vector exercise. Throughout the homework i will switch between vec and vector, but they are the same. 
struct vec {
    std::vector<double> data;

    vec() = default;
    explicit vec(int n) : data(n) {}
    vec(std::initializer_list<double> values) : data(values) {}
    vec(const vec&) = default;                                 // copy
    vec(vec&&) noexcept = default;                             // move
    ~vec() = default;                                          // destructor

    // assignment
    vec& operator=(const vec&) = default;                     // copy assignment
    vec& operator=(vec&&) noexcept = default;                 // move assignment

    int size() const { return static_cast<int>(data.size()); }
    double& operator[](int i) { return data[i]; }
    const double& operator[](int i) const { return data[i]; }

    double dot(const vec& other) const{
        assert(size() == other.size());
        double result = 0.0;
        for (int i = 0; i < size(); ++i)
            result += (*this)[i] * other[i];
        return result;
    }
    
    double norm() const{
        return std:: sqrt(dot(*this));
    }


    // arithmetic
    vec& operator+=(const vec& other){
        assert(size() == other.size());
        for (int i = 0; i < size(); ++i) (*this)[i] += other[i];
        return *this;

    }

    vec& operator-=(const vec& other){
        assert(size() == other.size());
        for (int i = 0; i < size(); ++i) (*this)[i] -= other[i];
        return *this;
    }

    vec& operator*=(double c){
        for (int i = 0; i < size(); ++i) (*this)[i] *= c;
        return *this;
    }

    vec& operator/=(double c) {
        for (int i=0; i<size(); ++i) (*this)[i] /= c;
        return *this;
    }

    //Print statements, setter statement and others:
    void print(const std::string& s = "") const {
        if (!s.empty()) {
            std::cout << s << '\n';
        }

        for (double x : data) {
            std::cout << std::setw(12) << x << ' ';
        }

        std::cout << '\n';
    }

};
    

// non-member operators
inline vec operator-(vec a, const vec& b) { a -= b; return a; }
inline vec operator-(vec a) { a*=-1.0; return a; }
inline vec operator+(vec a, const vec& b) { a += b; return a; }
inline vec operator*(const vec& a, double c) { vec result = a; result *= c; return result; }
inline vec operator*(double c, const vec& a) { return a * c; }
inline vec operator/(const vec& a, double c) { vec result = a; result /= c; return result; }



// approximate equality
inline bool approx(double x, double y, double acc = 1e-9, double eps = 1e-9){
    double d = std::abs(x - y);
    if (d <= acc) return true;
    return d <= eps * std::max(std::abs(x), std::abs(y));
}
inline bool approx(const vec& a, const vec& b, double acc = 1e-9, double eps = 1e-9){
    if (a.size() != b.size()) return false;
    for (int i = 0; i < a.size(); ++i){
        if (!approx(a[i], b[i], acc, eps)) return false;
    }
    return true;
}

struct matrix {
    std::vector<vec> M;

    matrix() = default;
    matrix(int rows, int columns) : M(columns, vec(rows)) {}

    matrix(const matrix&) = default;
    matrix(matrix&&) noexcept = default;
    matrix& operator=(const matrix&) = default;
    matrix& operator=(matrix&&) noexcept = default;

    int size1() const { return M.empty() ? 0 : M[0].size();} //gives the rows
    int size2() const {return M.size();} //Gives the columns


    double& operator()(int i, int j) {return M[j][i];}
    const double& operator()(int i, int j) const {return M[j][i];}


    vec& operator[](int j) {return M[j];}
    const vec& operator[](int j) const {return M[j];}

    void setid() {
        assert(size1() == size2());

        for (int i = 0; i < size1(); ++i)
            for (int j = 0; j < size2(); ++j)
                (*this)(i, j) = (i == j) ? 1.0 : 0.0;
    }

    matrix transpose() const {
        matrix AT(size2(), size1());

        for (int i = 0; i < size1(); ++i)
            for (int j = 0; j < size2(); ++j)
                AT(j, i) = (*this)(i, j);

        return AT;
    }

    matrix T() const {return transpose();}

    matrix& operator+=(const matrix& B) {
        assert(size1() == B.size1());
        assert(size2() == B.size2());

        for (int j = 0; j < size2(); ++j)
            M[j] += B[j];

        return *this;
    }

    matrix& operator-=(const matrix& B) {
        assert(size1() == B.size1());
        assert(size2() == B.size2());

        for (int j = 0; j < size2(); ++j)
            M[j] -= B[j];

        return *this;
    }

    matrix& operator*=(double c) {
        for (vec& column : M)
            column *= c;

        return *this;
    }

    matrix& operator/=(double c) {
        assert(c != 0.0);

        for (vec& column : M)
            column /= c;

        return *this;
    }

    void print(const std::string& s = "") const{
        if (!s.empty()) std::cout << s << std::endl;
        for (int i = 0; i < size1(); ++i){
            for (int j = 0; j < size2(); ++j)
                std::cout << std::setw(12) << (*this)(i, j) << " ";
            std::cout << std::endl;
        }
    }
};
// non-member operators

inline matrix operator+(matrix A, const matrix& B) { A += B; return A; }
inline matrix operator-(matrix A, const matrix& B) { A -= B; return A; }
inline matrix operator*(matrix A, double c)         { A *= c; return A; }
inline matrix operator*(double c, matrix A)         { A *= c; return A; }
inline matrix operator/(matrix A, double c)         { A /= c; return A; }

// matrix on vector and matrix on matrix multiplication
inline vec operator*(const matrix& A, const vec& x){
    assert(A.size2() == x.size());

    vec y(A.size1());

    for (int i = 0; i < A.size1(); ++i) {
        for (int j = 0; j < A.size2(); ++j)
            y[i] += A(i, j) * x[j];
    }
    return y;
}

inline double dot(const vec& a, const vec& b) {
    return a.dot(b);
}
inline matrix operator*(const matrix& A, const matrix& B) {
    assert(A.size2() == B.size1());
    matrix C(A.size1(), B.size2());

    for (int i = 0; i < B.size2(); ++i)
        C[i] = A * B[i];
    return C;
}

inline bool approx(const matrix& A, const matrix& B, double acc = 1e-9, double eps = 1e-9){
    if (A.size1() != B.size1() || A.size2() != B.size2()) return false;

    for (int i = 0; i < A.size1(); ++i)
        for (int j = 0; j < A.size2(); ++j)
            if (!approx(A(i, j), B(i, j), acc, eps)) return false;

    return true;
}

struct qr {
    //We want to do qr decomposition using a modified Gram-Schmidt
    matrix Q;
    matrix R;

    explicit qr(const matrix& A) : Q(A), R(A.size2(), A.size2()) {
        assert(A.size1() >= A.size2());
        for (int i = 0; i < A.size2(); ++i) {
        R(i,i) = Q[i].norm();
        if (approx(R(i,i), 0.0, 1e-9, 1e-9)) {
            // Handle the case where the norm is zero
            throw std::runtime_error("We cannot perform QR decomposition on a singular matrix.");
        }
        Q[i] /= R(i,i);

        for (int j = i + 1; j < A.size2(); ++j) {
            R(i, j) = dot(Q[i], Q[j]);
            Q[j] -= Q[i] * R(i, j);
        }
    }
}
    vec solve(const vec& b) const {
        assert(Q.size1() == b.size());
        assert(R.size1() == R.size2());

        vec x = Q.T() * b;

        for (int i = R.size1() - 1; i >= 0; --i) {
            for (int j = i + 1; j < R.size2(); ++j)
                x[i] -= R(i, j) * x[j];
            x[i] /= R(i, i);
        }
        return x;
    }

    double det() const{
        assert(R.size1() == R.size2());
        double product = 1.0;
        for (int i = 0; i < R.size1(); ++i) product *= R(i, i);
        return product;
    }

    matrix inverse() const {
        assert(Q.size1() == Q.size2());
        int n = Q.size1();
        matrix A(n,n);

        for (int i = 0; i < n; ++i) {
            vec e(n);
            e[i] = 1.0;
            A[i] = solve(e);
        }
        return A;
    }


}; 

inline matrix rando_M(int rows, int cols, double scale = 1.0, unsigned seed = 1) {
    std::mt19937 gen(seed);
    std::uniform_real_distribution<double> dist(-scale, scale);

    matrix MR(rows, cols);
    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j)
            MR(i, j) = dist(gen);
    return MR;

}

inline vec rando_vec(int n, double scale = 1.0, unsigned seed = 1) {
    std::mt19937 gen(seed);
    std::uniform_real_distribution<double> dist(-scale, scale);
    vec v(n);
    for (int i = 0; i < n; ++i)
        v[i] = dist(gen);
    return v;
}

inline bool uptrig(const matrix& A, double tol = 1e-10){
    assert(A.size1() == A.size2());

    for(int i = 0; i < A.size1(); ++i) 
        for(int j = 0; j < i && j < A.size2(); ++j)
            if (std::abs(A(i, j)) > tol) return false;
    return true;
}



}// namespace pp