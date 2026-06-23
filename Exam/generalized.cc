#include"generalized.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <stdexcept>

namespace pp{
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
} // namespace pp