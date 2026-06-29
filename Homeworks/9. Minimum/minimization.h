#pragma once

#include "matrix.h"
#include "vector.h"

#include <functional>

namespace pp {

using objective_function =
    std::function<double(const vec&)>;

struct minimization_result {
    vec x;
    double value;
    int steps;
    bool converged;
};

vec gradient(
    const objective_function& phi,
    const vec& x
);

matrix hessian(
    const objective_function& phi,
    const vec& x
);

minimization_result newton_minimize(
    const objective_function& phi,
    vec x,
    double acc = 1e-6,
    int max_steps = 1000
);

}