#pragma once

#include "matrix.h"

#include <functional>
#include <tuple>
#include <vector>

namespace pp {

using Func = std::function<double(double)>;
using Funcs = std::vector<Func>;


// QR decomposition A = Q R
void qr_decompose(
    const matrix& A,
    matrix& Q,
    matrix& R
);


// Solve A x = b using QR decomposition
vec qr_solve(
    const matrix& A,
    const vec& b
);


std::tuple<vec, matrix> lsfit(
    const Funcs& fs,
    const vec& x,
    const vec& y,
    const vec& dy
);


double eval(
    const Funcs& fs,
    const vec& c,
    double z
);

} // namespace pp