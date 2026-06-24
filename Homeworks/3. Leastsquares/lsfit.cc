#include "lsfit.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>


namespace pp{
    namespace {
        void qr_decompose(const matrix& A, matrix& Q, matrix& R){
            const int m = A.size2();

            Q = A;
            R = matrix(m, m);

            double scale = 0.0;

            for (int j = 0; j < m; ++j){
                scale = std::max(scale, A[j].norm());
            }

            //We modify the gram schmidt process:

            for (int i = 0; i < m; ++i){

                R(i,i) = Q[i].norm();

                Q[i] /= R(i,i);

                for (int j = i + 1; j < m; ++j){
                    R(i,j) = dot(Q[i], Q[j]);
                    Q[j] -= R(i,j) * Q[i];
                }
            }
        }
        //We want to solve the upper-triangular system Rx = b
        vec back_sub (const matrix& R, const vec& b){
            const int m = R.size1();
            vec x(m);

            for (int i = m - 1; i >= 0; --i){
                double sum = 0.0;

                for (int k = i+1; k < m; ++k){
                    sum += R(i,k) * x[k];
                } 
                x[i] = (b[i] - sum) / R(i,i);
            }
            return x;
        }
        //Want to find R^{-1} by findign R x_j = e_j for each column x_j of the inverse
        matrix inv_up_triang(const matrix& R){
            const int m = R.size1();
            
            matrix Rinv(m,m);

            for (int j = 0; j < m; ++j){
                vec e(m);
                e[j] = 1.0;
                Rinv[j] = back_sub(R, e);
                
            }
            return Rinv;
        }
    } // unnamed namespace
    double eval(const Funcs& fs, const vec& c,double z ){
        double sum = 0.0;

        for(int k = 0; k < c.size(); ++k){
            sum += c[k] * fs[k](z);
        }
        return sum;
    }
    //A fit to minimize chi2 

            std::tuple<vec,matrix> lsfit(const Funcs& fs, const vec& x, const vec& y, const vec& dy){
            const int n = x.size();
            const int m = static_cast<int>(fs.size());

            matrix A(n,m);
            vec b(n);

            for (int i = 0; i < n; ++i){
                b[i] = y[i] / dy[i];

                for (int k = 0; k < m; ++k){
                    A(i,k) = fs[k](x[i]) / dy[i];
                }
        }
        matrix Q;
        matrix R;
        qr_decompose(A, Q, R);
        const vec Qtb = Q.T() * b;
        const vec c = back_sub(R, Qtb);
        const matrix Rinv = inv_up_triang(R);
        const matrix cov = Rinv * Rinv.T();

        return {c,cov};
    }
} // namespace pp

