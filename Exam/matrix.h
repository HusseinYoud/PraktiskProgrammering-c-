#pragma once
#include "vector.h"
#include <string>

namespace pp{
    struct matrix{
        int Mrows = 0;
        int Mcolumns = 0;
        std::vector<vec> M;

        matrix () = default;

        matrix(int Mrows, int Mcoulmns);

        int size1() const;
        int size2() const;

        double& operator()(int i, int j);
        const double& operator()(int i,int j) const;

        vec& operator[](int j);
        const vec& operator[](int j) const;

        void setid();
        matrix transpose() const;
        matrix T() const;

        void print(const std:: string& s = "") const;

        matrix& operator+=(const matrix& B);
        matrix& operator-=(const matrix& B);
        matrix& operator*=(double c);
        matrix& operator/=(double c);

    };
    matrix operator+(matrix A, const matrix& B);
    matrix operator-(matrix A, const matrix& B);
    matrix operator*(matrix A, double c);
    matrix operator*(double c, matrix A);
    matrix operator/(matrix A, double c);


// matrix on vector and matrix on matrix multiplication
vec operator*(const matrix& A, const vec& x);
matrix operator*(const matrix& A, const matrix& B);

//linear algebra functions
matrix Id(int n);
matrix diag(const vec& v);
matrix rando_M(int rows, int cols, double scale = 1.0, unsigned seed = 1);
matrix rando_sym_M(int n, double scale = 1.0, unsigned seed = 1);

vec rando_vec(int n, double scale = 1.0, unsigned seed = 1);

bool approx(const matrix& A, const matrix& B, double acc = 1e-9, double eps = 1e-9);







} //namespace pp