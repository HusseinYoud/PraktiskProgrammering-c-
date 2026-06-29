#include "roots.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace pp {

namespace {

double column_norm(
    const matrix& A,
    int j
)
{
    double sum = 0.0;

    for (int i = 0; i < A.size1(); ++i) {
        sum += A(i, j) * A(i, j);
    }

    return std::sqrt(sum);
}

//entire qr decomp is copied from lsfit.cc + lsfit.h from homework 3
void qr_decomposition(
    const matrix& A,
    matrix& Q,
    matrix& R
)
{
    const int rows = A.size1();
    const int columns = A.size2();

    if (rows != columns) {
        throw std::invalid_argument(
            "qr_decomposition: A must be square"
        );
    }

    Q = A;
    R = matrix(columns, columns);

    for (int i = 0; i < columns; ++i) {

        const double norm = column_norm(Q, i);


        const double tolerance =
            100.0
            * std::numeric_limits<double>::epsilon();

        if (norm <= tolerance) {
            throw std::runtime_error(
                "qr_decomposition: matrix is singular"
            );
        }

        R(i, i) = norm;


        for (int k = 0; k < rows; ++k) {
            Q(k, i) /= R(i, i);
        }


        for (int j = i + 1; j < columns; ++j) {

            double dot_product = 0.0;

            for (int k = 0; k < rows; ++k) {
                dot_product += Q(k, i) * Q(k, j);
            }

            R(i, j) = dot_product;

            for (int k = 0; k < rows; ++k) {
                Q(k, j) -= Q(k, i) * R(i, j);
            }
        }
    }
}



vec back_substitution(
    const matrix& R,
    const vec& b
)
{
    const int n = R.size1();

    if (R.size2() != n || b.size() != n) {
        throw std::invalid_argument(
            "back_substitution: incompatible dimensions"
        );
    }

    vec x(n);

    for (int i = n - 1; i >= 0; --i) {

        double sum = 0.0;

        for (int j = i + 1; j < n; ++j) {
            sum += R(i, j) * x[j];
        }

        const double diagonal = R(i, i);

        const double tolerance =
            100.0
            * std::numeric_limits<double>::epsilon();

        if (std::abs(diagonal) <= tolerance) {
            throw std::runtime_error(
                "back_substitution: singular matrix"
            );
        }

        x[i] = (b[i] - sum) / diagonal;
    }

    return x;
}

//taken from lsfit.cc + lsfit.h from homework 3
vec qr_solve(
    const matrix& A,
    const vec& b
)
{
    if (A.size1() != A.size2()) {
        throw std::invalid_argument(
            "qr_solve: A must be square"
        );
    }

    if (A.size1() != b.size()) {
        throw std::invalid_argument(
            "qr_solve: incompatible dimensions"
        );
    }

    const int n = A.size1();

    matrix Q;
    matrix R;

    qr_decomposition(A, Q, R);


    vec c(n);

    for (int i = 0; i < n; ++i) {

        c[i] = 0.0;

        for (int k = 0; k < n; ++k) {
            c[i] += Q(k, i) * b[k];
        }
    }

    return back_substitution(R, c);
}


vec automatic_dx(
    const vec& x
)
{
    const double sqrt_machine_epsilon =
        std::ldexp(1.0, -26);

    vec dx(x.size());

    for (int i = 0; i < x.size(); ++i) {
        dx[i] =
            std::max(std::abs(x[i]), 1.0)
            * sqrt_machine_epsilon;
    }

    return dx;
}



bool is_empty(
    const vec& x
)
{
    return x.size() == 0;
}

}



matrix numerical_jacobian(
    const root_function& f,
    const vec& x,
    const vec& fx,
    const vec& dx
)
{
    const int n = x.size();

    if (fx.size() != n) {
        throw std::invalid_argument(
            "numerical_jacobian: f(x) must have the same size as x"
        );
    }

    if (dx.size() != n) {
        throw std::invalid_argument(
            "numerical_jacobian: dx must have the same size as x"
        );
    }

    matrix J(n, n);

    for (int j = 0; j < n; ++j) {

        if (dx[j] == 0.0) {
            throw std::invalid_argument(
                "numerical_jacobian: dx components must be non-zero"
            );
        }

  
        vec x_shifted = x;
        x_shifted[j] += dx[j];

        vec f_shifted = f(x_shifted);

        if (f_shifted.size() != n) {
            throw std::invalid_argument(
                "numerical_jacobian: f must return a vector "
                "with the same size as x"
            );
        }


        for (int i = 0; i < n; ++i) {
            J(i, j) =
                (f_shifted[i] - fx[i]) / dx[j];
        }
    }

    return J;
}


vec newton(
    const root_function& f,
    vec x,
    double acc,
    double alpha_min,
    int max_iterations,
    vec dx
)
{
    if (x.size() == 0) {
        throw std::invalid_argument(
            "newton: initial vector must not be empty"
        );
    }

    if (acc <= 0.0) {
        throw std::invalid_argument(
            "newton: acc must be positive"
        );
    }

    if (alpha_min <= 0.0 || alpha_min > 1.0) {
        throw std::invalid_argument(
            "newton: alpha_min must satisfy 0 < alpha_min <= 1"
        );
    }

    if (max_iterations <= 0) {
        throw std::invalid_argument(
            "newton: max_iterations must be positive"
        );
    }

    const int n = x.size();


    const bool use_automatic_dx = is_empty(dx);

    if (!use_automatic_dx && dx.size() != n) {
        throw std::invalid_argument(
            "newton: dx must have the same size as x"
        );
    }

    vec fx = f(x);

    if (fx.size() != n) {
        throw std::invalid_argument(
            "newton: f must return a vector "
            "with the same size as x"
        );
    }

    for (
        int iteration = 0;
        iteration < max_iterations;
        ++iteration
    ) {

        const double fx_norm = fx.norm();

        if (fx_norm < acc) {
            return x;
        }


        vec current_dx =
            use_automatic_dx
            ? automatic_dx(x)
            : dx;

        matrix J =
            numerical_jacobian(
                f,
                x,
                fx,
                current_dx
            );


        vec minus_fx(n);

        for (int i = 0; i < n; ++i) {
            minus_fx[i] = -fx[i];
        }

        vec Dx = qr_solve(J, minus_fx);


        double alpha = 1.0;

        vec z(n);
        vec fz(n);

        while (true) {

            for (int i = 0; i < n; ++i) {
                z[i] = x[i] + alpha * Dx[i];
            }

            fz = f(z);

            if (fz.size() != n) {
                throw std::invalid_argument(
                    "newton: f must return a vector "
                    "with the same size as x"
                );
            }

            const double fz_norm = fz.norm();


            if (
                fz_norm
                < (1.0 - alpha / 2.0) * fx_norm
            ) {
                break;
            }

            if (fz_norm < fx_norm) {
                break;
            }


            if (alpha <= alpha_min) {
                break;
            }

            alpha /= 2.0;

            if (alpha < alpha_min) {
                alpha = alpha_min;
            }
        }


        vec actual_step(n);

        for (int i = 0; i < n; ++i) {
            actual_step[i] = alpha * Dx[i];
        }


        x = z;
        fx = fz;


        if (fx.norm() < acc) {
            return x;
        }

        if (actual_step.norm() < current_dx.norm()) {
            return x;
        }
    }

    throw std::runtime_error(
        "newton: maximum number of iterations reached"
    );
}

}