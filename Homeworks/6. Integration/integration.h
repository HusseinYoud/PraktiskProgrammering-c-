#pragma once

#include <functional>

namespace pp {

using integrand =
    std::function<double(double)>;

struct integration_result {
    double value;
    double error;
};


double integrate(
    const integrand& f,
    double a,
    double b,
    double acc = 1e-3,
    double eps = 1e-3
);


integration_result integrate_with_error(
    const integrand& f,
    double a,
    double b,
    double acc = 1e-3,
    double eps = 1e-3
);


double clenshaw_curtis(
    const integrand& f,
    double a,
    double b,
    double acc = 1e-3,
    double eps = 1e-3
);


integration_result clenshaw_curtis_with_error(
    const integrand& f,
    double a,
    double b,
    double acc = 1e-3,
    double eps = 1e-3
);


double erf(
    double z,
    double acc = 1e-3,
    double eps = 1e-3
);

} // namespace pp