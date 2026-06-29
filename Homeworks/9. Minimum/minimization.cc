#include "minimization.h"
#include "lsfit.h"
#include <cmath>
#include <limits>
#include <stdexcept>

namespace pp {

namespace {

constexpr double armijo_alpha = 1e-4;
constexpr double lambda_min = 1.0 / 1024.0;

constexpr double initial_regularization = 1e-6;
constexpr int max_regularization_attempts = 12;



double gradient_step(double xi)
{    const double epsilon =
        std::numeric_limits<double>::epsilon();
    return (1.0 + std::abs(xi))
         * std::sqrt(epsilon);
}



double hessian_step(double xi)
{
    const double epsilon =
        std::numeric_limits<double>::epsilon();
    return (1.0 + std::abs(xi))
         * std::pow(epsilon, 0.25);
}


bool is_finite(const vec& x)
{
    for (int i = 0; i < x.size(); ++i) {
        if (!std::isfinite(x[i])) {
            return false;}
    }
    return true;
}

vec gradient_from_value(const objective_function& phi,const vec& x, double phi_x)
{
    const int n = x.size();
    vec g(n);
    vec x_shifted = x;
    for (int i = 0; i < n; ++i) {
        const double dx = gradient_step(x[i]);
        x_shifted[i] = x[i] + dx;
        const double phi_shifted =
            phi(x_shifted);
        g[i] = (phi_shifted - phi_x) / dx;
        x_shifted[i] = x[i];
    }
    return g;
}

vec negative(const vec& x)
{
    vec result(x.size());
    for (int i = 0; i < x.size(); ++i) {
        result[i] = -x[i];
    }
    return result;
}

matrix hessian_from_gradient(const objective_function& phi, const vec& x, const vec& gradient_at_x)
{
    const int n = x.size();

    matrix H(n, n);
    vec x_shifted = x;

    for (int j = 0; j < n; ++j) {
        const double dx = hessian_step(x[j]);

        x_shifted[j] = x[j] + dx;
        const double phi_shifted =
            phi(x_shifted);

        const vec shifted_gradient =
            gradient_from_value(phi, x_shifted, phi_shifted);

        for (int i = 0; i < n; ++i) {
            H(i, j) =
                (
                    shifted_gradient[i]
                    - gradient_at_x[i]
                ) / dx;
        }

        x_shifted[j] = x[j];
    }


    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            const double average =
                0.5 * (H(i, j) + H(j, i));

            H(i, j) = average;
            H(j, i) = average;
        }
    }
    return H;
}



vec make_trial_point(
    const vec& x,
    const vec& dx,
    double lambda
)
{
    vec trial = x;

    for (int i = 0; i < x.size(); ++i) {
        trial[i] += lambda * dx[i];
    }

    return trial;
}
} // namespace of unnamed 


vec gradient(
    const objective_function& phi,
    const vec& x
)
{
    if (x.size() == 0) {
        throw std::invalid_argument(
            "gradient: x must not be empty"
        );
    }

    const double phi_x = phi(x);

    if (!std::isfinite(phi_x)) {
        throw std::runtime_error(
            "gradient: phi(x) is not finite"
        );
    }

    return gradient_from_value(phi, x, phi_x);
}


matrix hessian(
    const objective_function& phi,
    const vec& x
)
{
    if (x.size() == 0) {
        throw std::invalid_argument(
            "hessian: x must not be empty"
        );
    }

    const double phi_x = phi(x);

    if (!std::isfinite(phi_x)) {
        throw std::runtime_error(
            "hessian: phi(x) is not finite"
        );
    }

    const vec g =
        gradient_from_value(phi, x, phi_x);

    return hessian_from_gradient(phi, x, g);
}


minimization_result newton_minimize(
    const objective_function& phi,
    vec x,
    double acc,
    int max_steps
)
{
    if (x.size() == 0) {
        throw std::invalid_argument(
            "newton_minimize: x must not be empty"
        );
    }

    if (acc <= 0.0) {
        throw std::invalid_argument(
            "newton_minimize: acc must be positive"
        );
    }

    if (max_steps < 0) {
        throw std::invalid_argument(
            "newton_minimize: max_steps "
            "must be non-negative"
        );
    }

    const int n = x.size();

    double phi_x = phi(x);

    if (!std::isfinite(phi_x)) {
        throw std::runtime_error(
            "newton_minimize: phi(x) is not finite "
            "at the initial point"
        );
    }

    for (int steps = 0; steps < max_steps; ++steps) {

        const vec g =
            gradient_from_value(phi, x, phi_x);


        if (g.norm() < acc) {
            return {
                x,
                phi_x,
                steps,
                true
            };
        }

        const matrix H =
            hessian_from_gradient(phi, x, g);


        vec dx(n);

        bool direction_found = false;
        double mu = initial_regularization;

        const vec minus_g = negative(g);

        for (
            int attempt = 0;
            attempt < max_regularization_attempts;
            ++attempt
        ) {
            matrix H_regularized = H;

            for (int i = 0; i < n; ++i) {
                H_regularized(i, i) += mu;
            }

            try {
                const vec candidate =
                    qr_solve(H_regularized, minus_g);


                if (
                    is_finite(candidate)
                    && dot(g, candidate) < 0.0
                ) {
                    dx = candidate;
                    direction_found = true;
                    break;
                }
            }
            catch (const std::runtime_error&) {

            }

            mu *= 10.0;
        }


        if (!direction_found) {
            dx = minus_g;
        }

        const double slope = dot(g, dx);


        double lambda = 1.0;

        vec x_trial = x;
        double phi_trial = phi_x;

        bool accepted = false;

        while (lambda >= lambda_min) {
            x_trial =
                make_trial_point(x, dx, lambda);

            phi_trial = phi(x_trial);

            const double armijo_limit =
                phi_x
                + armijo_alpha
                * lambda
                * slope;

            if (
                std::isfinite(phi_trial)
                && phi_trial < armijo_limit
            ) {
                accepted = true;
                break;
            }

            lambda *= 0.5;
        }


        if (!accepted) {
            lambda = lambda_min;

            x_trial =
                make_trial_point(x, dx, lambda);

            phi_trial = phi(x_trial);

            if (!std::isfinite(phi_trial)) {
                return {
                    x,
                    phi_x,
                    steps,
                    false
                };
            }
        }

        x = x_trial;
        phi_x = phi_trial;
    }


    const vec final_gradient =
        gradient_from_value(phi, x, phi_x);

    return {
        x,
        phi_x,
        max_steps,
        final_gradient.norm() < acc
    };
}


} // namespace pp