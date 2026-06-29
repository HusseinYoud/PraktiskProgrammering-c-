#include "lsfit.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace pp {

namespace {


// Solve the upper-triangular system R x = b
vec back_sub(
    const matrix& R,
    const vec& b
)
{
    const int m = R.size1();
    vec x(m);

    for (int i = m - 1; i >= 0; --i) {
        double sum = 0.0;

        for (int k = i + 1; k < m; ++k) {
            sum += R(i, k) * x[k];
        }

        x[i] = (b[i] - sum) / R(i, i);
    }

    return x;
}


// Find R^{-1} by solving R x_j = e_j
matrix inv_up_triang(const matrix& R)
{
    const int m = R.size1();

    matrix Rinv(m, m);

    for (int j = 0; j < m; ++j) {
        vec e(m);
        e[j] = 1.0;

        Rinv[j] = back_sub(R, e);
    }

    return Rinv;
}


} // namespace


void qr_decompose(
    const matrix& A,
    matrix& Q,
    matrix& R
)
{
    const int m = A.size2();

    Q = A;
    R = matrix(m, m);

    double scale = 0.0;

    for (int j = 0; j < m; ++j) {
        scale = std::max(scale, A[j].norm());
    }

    const double tolerance =
        std::numeric_limits<double>::epsilon()
        * std::max(1.0, scale);

    for (int i = 0; i < m; ++i) {
        R(i, i) = Q[i].norm();

        if (R(i, i) <= tolerance) {
            throw std::runtime_error(
                "qr_decompose: matrix is singular "
                "or nearly singular"
            );
        }

        Q[i] /= R(i, i);

        for (int j = i + 1; j < m; ++j) {
            R(i, j) = dot(Q[i], Q[j]);
            Q[j] -= R(i, j) * Q[i];
        }
    }
}


vec qr_solve(
    const matrix& A,
    const vec& b
)
{
    if (A.size1() != b.size()) {
        throw std::invalid_argument(
            "qr_solve: incompatible matrix and vector dimensions"
        );
    }

    matrix Q;
    matrix R;

    qr_decompose(A, Q, R);

    const vec Qtb = Q.T() * b;

    return back_sub(R, Qtb);
}


double eval(
    const Funcs& fs,
    const vec& c,
    double z
)
{
    double sum = 0.0;

    for (int k = 0; k < c.size(); ++k) {
        sum += c[k] * fs[k](z);
    }

    return sum;
}


std::tuple<vec, matrix> lsfit(
    const Funcs& fs,
    const vec& x,
    const vec& y,
    const vec& dy
)
{
    const int n = x.size();
    const int m = static_cast<int>(fs.size());

    matrix A(n, m);
    vec b(n);

    for (int i = 0; i < n; ++i) {
        b[i] = y[i] / dy[i];

        for (int k = 0; k < m; ++k) {
            A(i, k) = fs[k](x[i]) / dy[i];
        }
    }

    matrix Q;
    matrix R;

    qr_decompose(A, Q, R);

    const vec Qtb = Q.T() * b;
    const vec c = back_sub(R, Qtb);

    const matrix Rinv = inv_up_triang(R);
    const matrix cov = Rinv * Rinv.T();

    return {c, cov};
}


} // namespace pp