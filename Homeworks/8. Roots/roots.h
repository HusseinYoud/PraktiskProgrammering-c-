#pragma once

#include "matrix.h"
#include "vector.h"

#include <functional>

namespace pp {

using root_function =
    std::function<vec(const vec&)>;


matrix numerical_jacobian(
    const root_function& f,
    const vec& x,
    const vec& fx,
    const vec& dx
);


vec newton(
    const root_function& f,
    vec x,
    double acc = 1e-12,
    double alpha_min = 1e-3,
    int max_iterations = 100,
    vec dx = {}
);

}