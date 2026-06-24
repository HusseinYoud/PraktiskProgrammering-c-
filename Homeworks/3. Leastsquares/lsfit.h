#pragma once
#include "matrix.h"
#include <tuple>
#include <functional>
#include <vector>

namespace pp {

using Func = std::function<double(double)>;
using Funcs = std::vector<Func>;

std::tuple<vec, matrix> lsfit(const Funcs& fs, const vec& x, const vec& y, const vec& dy);

// struct fit_result {
//     vec c;          // Best-fit coefficients
//     matrix cov;     // Covariance matrix
//     double chi2 = 0.0;
//     int dof = 0;
// };

//fit_result lsfit(const Funcs& fs,const vec& x, const vec& y, const vec& dy);

double eval(const Funcs& fs, const vec& c,double z );

} // namespace pp