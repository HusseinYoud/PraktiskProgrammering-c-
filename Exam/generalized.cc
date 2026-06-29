#include"generalized.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <stdexcept>

//GEVD inspired by earlier exercise. 

namespace pp{
    namespace {
        matrix solve_lower (const matrix& L, const matrix& R){
            const int n = L.size1();
            const int m = R.size2();
        
        matrix X(n,m);

        for (int j = 0; j < m; ++j){
            for (int i = 0; i < n; ++i){
                double sum = R(i,j);
                for (int k = 0; k < i; ++k) {
                    sum -= L(i,k) * X(k,j);
                }
                X(i,j) = sum / L(i,i);
            }
        }
        return X;
    }
    matrix solve_lower_transpose(const matrix& L, const matrix& R){
        const int n = L.size1();
        const int m = R.size2();
    
    matrix X(n,m);

    for(int j = 0; j < m; ++j){
        for(int i = n - 1;i >= 0; --i){
            double sum = R(i,j);

            for(int k = i + 1; k < n; ++k){
                sum -= L(k,i) * X(k,j);
            }
            X(i,j) = sum / L(i,i);
        }
    }
    return X;
    }
    }  
    GEVD::GEVD(const matrix& A, const matrix& B) : w(A.size1()), V(A.size1(), A.size2()) {
        assert(A.size1() == A.size2());
        assert(B.size1() == B.size2());
        assert(A.size1() == B.size1());


        const int n = A.size1();

        // Solving B = QS Q^T
        EVD B_evd(B);
        const vec& s= B_evd.w;
        const matrix& Q = B_evd.V;


        //We now can calculate the diagonal elements of S^(-1/2)

        vec invsqrt_s(n);

        double large_s = 0.0;

        for (int i = 0; i < n; ++i) {
            large_s = std::max(large_s, std::abs(s[i]));
        }
        if (large_s <= 0.0) {
            throw std::invalid_argument("Matrix not defentive");
        }
        const double pd_tol = 1e-12 * large_s;

        for (int i = 0; i < n; ++i) {
            if (s[i] <= pd_tol) {
                throw std::invalid_argument("Matrix not positive");
            }
            invsqrt_s[i] = 1.0 / std::sqrt(s[i]);
        }

        //We now can calculate C = Q^T*A*Q
        matrix C = Q.T() * A * Q;

        //This leads us to find A_tilde = S^(-1/2)*C*S(-1/2)

        matrix A_tilde(n,n);
        for (int i = 0; i < n; ++i) {
            for (int j = i; j < n; ++j) {
                const double Cij = 0.5*(C(i,j) + C(j,i));
                const double val = invsqrt_s[i] * Cij * invsqrt_s[j];
                A_tilde(i,j) = val;
                A_tilde(j,i) = val;
            }
        }
        //now we solve A_tilde X = X E
        EVD Atilde_evd(A_tilde);

        w = Atilde_evd.w;
        const matrix& X = Atilde_evd.V;

        //We can now calculate Y = S^(-1/2) X by scaling the rows of X.
        matrix Y(n,n);

        for(int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                Y(i,j) = invsqrt_s[i] * X(i,j);
            }
        }
        V = Q * Y;
    }
    matrix cholesky(const matrix& B){
        const int n = B.size1();
        matrix L(n, n);
        double scale = 0.0;

        for (int i = 0; i < n; ++i){
            scale = std::max(scale,std::abs(B(i,i)));
        }
        const double tol = 1e-14 * std::max(1.0,scale);

        for (int j = 0; j < n; ++j){
            double diagonal = B(j,j);
            for (int k = 0; k < j; ++k){
                diagonal -= L(j,k) * L(j,k);
            }
        
        L(j,j) = std::sqrt(diagonal);
        for(int i = j + 1; i < n; ++i){
            double sum = B(i,j);
            for (int k = 0; k < j; ++k){
                sum -= L(i,k) *L(j,k);
            }
            L(i,j) = sum / L(j,j);
         }
        }
        return L;
    }
    GEVD_cholesky::GEVD_cholesky(const matrix&A, const matrix& B) : w(A.size1()), V(A.size1(), A.size2()){
        const int n = A.size1();

        const matrix L = cholesky(B);

        const matrix Z = solve_lower(L,A);

        matrix C = solve_lower(L,Z.T()).T();

        for(int i = 0; i < n; ++i){
            for (int j = i; j < n; ++j){
                const double value = 0.5 * (C(i,j) + C(j,i));
                C(i,j) = value;
                C(j,i) = value;
            }
        }
        EVD transformed_evd(C);
        w = transformed_evd.w;

        const matrix& Y = transformed_evd.V;

        V = solve_lower_transpose(L,Y);
    }
        
} // namespace pp