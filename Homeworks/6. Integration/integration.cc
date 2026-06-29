#include "integration.h"

#include <cmath>
#include <limits>
#include <numbers>
#include <stdexcept>

namespace {

void check_accuracy(
    double acc,
    double eps
)
{
    if (acc < 0.0) {
        throw std::invalid_argument(
            "acc must be non-negative"
        );
    }

    if (eps < 0.0) {
        throw std::invalid_argument(
            "eps must be non-negative"
        );
    }
}



pp::integration_result integrate_recursive(
    const pp::integrand& f,
    double a,
    double b,
    double acc,
    double eps,
    double f2,
    double f3
)
{
    const double h =
        b - a;

    const double f1 =
        f(a + h / 6.0);

    const double f4 =
        f(a + 5.0 * h / 6.0);


    const double Q =
        h
        * (2.0 * f1 + f2 + f3 + 2.0 * f4)
        / 6.0;


    const double q =
        h
        * (f1 + f2 + f3 + f4)
        / 4.0;

    const double error =
        std::abs(Q - q);

    const double tolerance =
        acc + eps * std::abs(Q);

    if (error <= tolerance) {
        return {
            Q,
            error
        };
    }


    const double midpoint =
        (a + b) / 2.0;


    if (midpoint == a || midpoint == b) {
        return {
            Q,
            error
        };
    }

    const double sub_acc =
        acc / std::sqrt(2.0);

    const pp::integration_result left =
        integrate_recursive(
            f,
            a,
            midpoint,
            sub_acc,
            eps,
            f1,
            f2
        );

    const pp::integration_result right =
        integrate_recursive(
            f,
            midpoint,
            b,
            sub_acc,
            eps,
            f3,
            f4
        );


    return {
        left.value + right.value,
        std::hypot(
            left.error,
            right.error
        )
    };
}


pp::integration_result clenshaw_curtis_finite(
    const pp::integrand& f,
    double a,
    double b,
    double acc,
    double eps
)
{
    const double width =
        b - a;

    const double half_width =
        width / 2.0;

    const pp::integrand transformed =
        [&f, a, width, half_width](double theta) {
            const double half_theta =
                theta / 2.0;

            const double c =
                std::cos(half_theta);


            const double x =
                a + width * c * c;

            const double jacobian =
                half_width * std::sin(theta);

            return f(x) * jacobian;
        };

    return pp::integrate_with_error(
        transformed,
        0.0,
        std::numbers::pi,
        acc,
        eps
    );
}

} // namespace


namespace pp {


integration_result integrate_with_error(
    const integrand& f,
    double a,
    double b,
    double acc,
    double eps
)
{
    check_accuracy(
        acc,
        eps
    );

 
    if (!std::isfinite(a) ||
        !std::isfinite(b)) {

        throw std::invalid_argument(
            "pp::integrate_with_error requires finite limits"
        );
    }

    if (a == b) {
        return {
            0.0,
            0.0
        };
    }


    if (a > b) {
        integration_result result =
            integrate_with_error(
                f,
                b,
                a,
                acc,
                eps
            );

        result.value =
            -result.value;

        return result;
    }

    const double h =
        b - a;


    const double f2 =
        f(a + 2.0 * h / 6.0);

    const double f3 =
        f(a + 4.0 * h / 6.0);

    return integrate_recursive(
        f,
        a,
        b,
        acc,
        eps,
        f2,
        f3
    );
}



double integrate(
    const integrand& f,
    double a,
    double b,
    double acc,
    double eps
)
{
    return integrate_with_error(
        f,
        a,
        b,
        acc,
        eps
    ).value;
}



integration_result clenshaw_curtis_with_error(
    const integrand& f,
    double a,
    double b,
    double acc,
    double eps
)
{
    check_accuracy(
        acc,
        eps
    );

    if (std::isnan(a) ||
        std::isnan(b)) {

        throw std::invalid_argument(
            "integration limits cannot be NaN"
        );
    }


    if (a == b) {
        if (std::isfinite(a)) {
            return {
                0.0,
                0.0
            };
        }

        throw std::invalid_argument(
            "equal infinite limits are undefined"
        );
    }


    if (a > b) {
        integration_result result =
            clenshaw_curtis_with_error(
                f,
                b,
                a,
                acc,
                eps
            );

        result.value =
            -result.value;

        return result;
    }

 
    if (std::isfinite(a) &&
        std::isfinite(b)) {

        return clenshaw_curtis_finite(
            f,
            a,
            b,
            acc,
            eps
        );
    }


    if (std::isfinite(a) &&
        std::isinf(b) &&
        b > 0.0) {

        const integrand transformed =
            [&f, a](double t) {
                const double denominator =
                    1.0 - t;

                const double x =
                    a + t / denominator;

                const double jacobian =
                    1.0
                    / (denominator * denominator);

                return f(x) * jacobian;
            };

        return clenshaw_curtis_finite(
            transformed,
            0.0,
            1.0,
            acc,
            eps
        );
    }


    if (std::isinf(a) &&
        a < 0.0 &&
        std::isfinite(b)) {

        const integrand transformed =
            [&f, b](double t) {
                const double denominator =
                    1.0 - t;

                const double x =
                    b - t / denominator;

                const double jacobian =
                    1.0
                    / (denominator * denominator);

                return f(x) * jacobian;
            };

        return clenshaw_curtis_finite(
            transformed,
            0.0,
            1.0,
            acc,
            eps
        );
    }

 
    if (std::isinf(a) &&
        a < 0.0 &&
        std::isinf(b) &&
        b > 0.0) {

        const integrand transformed =
            [&f](double t) {
                const double denominator =
                    1.0 - t * t;

                const double x =
                    t / denominator;

                const double jacobian =
                    (1.0 + t * t)
                    /
                    (denominator * denominator);

                return f(x) * jacobian;
            };

        return clenshaw_curtis_finite(
            transformed,
            -1.0,
            1.0,
            acc,
            eps
        );
    }

    throw std::invalid_argument(
        "unsupported integration limits"
    );
}



double clenshaw_curtis(
    const integrand& f,
    double a,
    double b,
    double acc,
    double eps
)
{
    return clenshaw_curtis_with_error(
        f,
        a,
        b,
        acc,
        eps
    ).value;
}



double erf(
    double z,
    double acc,
    double eps
)
{
    const double factor =
        2.0
        / std::sqrt(std::numbers::pi);

    if (z < 0.0) {
        return -erf(
            -z,
            acc,
            eps
        );
    }

    if (z <= 1.0) {
        const integrand gaussian =
            [](double x) {
                return std::exp(-x * x);
            };

        return factor
            * integrate(
                gaussian,
                0.0,
                z,
                acc,
                eps
            );
    }

    const integrand transformed =
        [z](double t) {
            const double x =
                z + (1.0 - t) / t;

            return std::exp(-x * x)
                / (t * t);
        };

    return 1.0
        - factor
        * integrate(
            transformed,
            0.0,
            1.0,
            acc,
            eps
        );
}

} // namespace pp