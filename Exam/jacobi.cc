#include "jacobi.h"
#include <cassert>
#include <cmath>
#include <stdexcept>

//jacobi.cc is also taken from homework 2 EVD
namespace pp{
    void EVD::timesJ(pp::matrix& A, int p, int q, double theta) {
        assert(p >= 0 && p < A.size2());
        assert(q >= 0 && q < A.size2());
        assert(p != q);

        const double c = std::cos(theta);
        const double s = std::sin(theta);

        for (int i = 0; i < A.size1(); ++i) {
            const double Aip = A(i, p);
            const double Aiq = A(i, q);

            A(i, p) = c * Aip - s * Aiq;
            A(i, q) = s * Aip + c * Aiq;
        }
    }

    void EVD::Jtimes(pp::matrix& A, int p, int q, double theta){
        assert(p >= 0 && p < A.size1());
        assert(q >= 0 && q < A.size1());
        assert(p != q);

        const double c = std::cos(theta);
        const double s = std::sin(theta);

        for (int j = 0; j < A.size2(); ++j){
            const double Apj = A(p, j);
            const double Aqj = A(q, j);

            A(p, j) = c * Apj + s * Aqj;
            A(q, j) = -s * Apj + c * Aqj;
        }


    }

    //A bit of a crazy implementation but the idea is that we create a copy of the input matrix when the constructor starts.
    //Basically it asks has any diagonal element changed after a full rotation, then it goes to changed = true which adds the element to matrix A and matrix V, if changed = true it asks if it needs to make another sweep or not. If yes it makes another sweep if no it has converged.
    EVD::EVD(matrix A) : w(A.size1()), V(A.size1(), A.size2()){

        assert(A.size1() == A.size2());
        const int n = A.size1();
        V.setid();
        const double tolerance = 1e-12;
        const int max_sweeps = 100;

        bool converged = false;
        for (int sweep = 0; sweep < max_sweeps; ++sweep) {
            bool changed = false;

            for (int p = 0; p < n - 1; ++p){
                for (int q = p + 1; q < n; ++q) {
                    const double App = A(p,p);
                    const double Aqq = A(q,q);
                    const double Apq = A(p,q);


                    const double scale = std::abs(App) + std::abs(Aqq) + 1.0;

                    if (std::abs(Apq) <= tolerance * scale) {
                        continue;
                    }
                    const double theta = 0.5 * std::atan2(2.0 * Apq, Aqq - App);

                    timesJ(A, p, q, theta);
                    Jtimes(A, p, q, -theta);
                    timesJ(V, p, q, theta);
                    changed = true;
                }
            }
            if (!changed) {
                converged = true;
                break;
            }
        }
        if (!converged) {
            throw std::runtime_error("Jacobi did not converge");
        }
        for (int i = 0; i < n; ++i) {
            w[i] = A(i, i);
        }
    } 
}//namespace pp