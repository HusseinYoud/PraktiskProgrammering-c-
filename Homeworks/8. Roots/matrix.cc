#include "matrix.h"
#include<iomanip>
#include<iostream>
#include<random>
#include<stdexcept>
#include<cassert>

//Matrix.h, Vector.h, vector.cc and matrix.cc are all copied from earlier homeworks. (ode too)
namespace pp{
    matrix::matrix(int rows, int columns) : Mrows(rows), Mcolumns(columns), M(Mcolumns, vec(Mrows)) {}

    int matrix::size1() const {return Mrows;}
    int matrix::size2() const {return Mcolumns;}

    double& matrix::operator()(int i, int j) {return M[j][i];}
    const double& matrix::operator()(int i, int j) const {return M[j][i];}
    

    vec& matrix::operator[](int j) {return M[j];}
    const vec& matrix::operator[](int j) const {return M[j];}


    void matrix::setid() {
        assert(size1() == size2());

        for (int i = 0; i < size1(); ++i)
            for (int j = 0; j < size2(); ++j)
                (*this)(i, j) = (i == j) ? 1.0 : 0.0;
    }
    matrix matrix::transpose() const {
        matrix AT(size2(), size1());

        for (int i = 0; i < size1(); ++i)
            for (int j = 0; j < size2(); ++j)
                AT(j, i) = (*this)(i, j);

        return AT;
    }
    matrix matrix::T() const {return transpose();}

    matrix& matrix::operator+=(const matrix& B) {
        assert(size1() == B.size1());
        assert(size2() == B.size2());

        for (int j = 0; j < size2(); ++j)
            M[j] += B[j];

        return *this;
    }

    matrix& matrix::operator-=(const matrix& B) {
        assert(size1() == B.size1());
        assert(size2() == B.size2());

        for (int j = 0; j < size2(); ++j)
            M[j] -= B[j];

        return *this;
    }

    matrix& matrix::operator*=(double c) {
        for (vec& column : M)
            column *= c;

        return *this;
    }

    matrix& matrix::operator/=(double c) {
        assert(c != 0.0);

        for (vec& column : M)
            column /= c;
        return *this;
    }
    

    void matrix::print(const std::string& s) const{
        if (!s.empty()) std::cout << s << std::endl;
        for (int i = 0; i < size1(); ++i){
            for (int j = 0; j < size2(); ++j)
                std::cout << std::setw(12) << (*this)(i, j) << " ";
            std::cout << std::endl;
        }
    }
    matrix operator+(matrix A, const matrix& B) { A += B; return A; }
    matrix operator-(matrix A, const matrix& B) { A -= B; return A; }
    matrix operator*(matrix A, double c)         { A *= c; return A; }
    matrix operator*(double c, matrix A)         { A *= c; return A; }
    matrix operator/(matrix A, double c)         { A /= c; return A; }

    vec operator*(const matrix& A, const vec& x){
        assert(A.size2() == x.size());

        vec y(A.size1());

        for (int i = 0; i < A.size1(); ++i) {
            for (int j = 0; j < A.size2(); ++j)
                y[i] += A(i, j) * x[j];
        }
        return y;
    }

    matrix operator*(const matrix& A, const matrix& B) {
        assert(A.size2() == B.size1());
        matrix C(A.size1(), B.size2());

        for (int i = 0; i < B.size2(); ++i)
            C[i] = A * B[i];
        return C;
    }
    matrix Id(int n) {
        matrix I(n, n);
        I.setid();
        return I;
    }

    matrix diag(const vec& v) {
        int n = v.size();
        matrix D(n, n);
        for (int i = 0; i < n; ++i)
            D(i, i) = v[i];
        return D;
    }



    matrix rando_M(int rows, int cols, double scale, unsigned seed) {
        std::mt19937 gen(seed);
        std::uniform_real_distribution<double> dist(-scale, scale);

        matrix MR(rows, cols);
        for (int i = 0; i < rows; ++i)
            for (int j = 0; j < cols; ++j)
                MR(i, j) = dist(gen);
        return MR;

    }

    matrix rando_sym_M(int n, double scale, unsigned seed){
        matrix RM = rando_M(n,n,scale,seed);
        matrix RMsym(n, n);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                RMsym(i, j) = 0.5 * (RM(i, j) + RM(j, i));
            }
        }
        return RMsym;
    }

    vec rando_vec(int n, double scale, unsigned seed) {
        std::mt19937 gen(seed);
        std::uniform_real_distribution<double> dist(-scale, scale);
        vec v(n);
        for (int i = 0; i < n; ++i)
            v[i] = dist(gen);
        return v;
    }

    bool approx(const matrix& A, const matrix& B, double acc, double eps){
    if (A.size1() != B.size1() || A.size2() != B.size2()) return false;

    for (int i = 0; i < A.size1(); ++i)
        for (int j = 0; j < A.size2(); ++j)
            if (!approx(A(i, j), B(i, j), acc, eps)) return false;

    return true;
}
} //namespace pp