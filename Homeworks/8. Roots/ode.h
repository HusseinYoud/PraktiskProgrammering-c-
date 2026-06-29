#pragma once
#include "vector.h"
#include <functional>
#include <limits>
#include <tuple>
#include <vector>

namespace pp {

using ode_function =
    std::function<vec(double, const vec&)>;

std::tuple<vec, vec> rkstep12(
    const ode_function& f,
    double x,
    const vec& y,
    double h
);

std::tuple<std::vector<double>, std::vector<vec>> driver(
    const ode_function& f,
    double a,
    const vec& yinit,
    double b,
    double h = 0.125,
    double acc = 0.01,
    double eps = 0.01,
    double hmax = std::numeric_limits<double>::infinity()
);

}