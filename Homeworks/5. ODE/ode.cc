#include "ode.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace pp {

    std::tuple<vec, vec> rkstep12(
        const ode_function& f,
        double x,
        const vec& y,
        double h
    )
    {
        vec k0 = f(x, y);

        vec k1 = f(
            x + h / 2.0,
            y + k0 * (h / 2.0)
        );

        vec yh = y + k1 * h;
        vec dy = (k1 - k0) * h;

        return {yh, dy};
    }


    std::tuple<std::vector<double>, std::vector<vec>> driver(
        const ode_function& f,
        double a,
        const vec& yinit,
        double b,
        double h,
        double acc,
        double eps
    )
    {

        double x = a;
        vec y = yinit;

        std::vector<double> xlist;
        std::vector<vec> ylist;

        xlist.push_back(x);
        ylist.push_back(y);

        if (a == b) {
            return {xlist, ylist};
        }

        while (x < b) {
            if (x + h > b) {
                h = b - x;
            }

            auto [yh, dy] = rkstep12(f, x, y, h);

            double tol =
                (acc + eps * yh.norm())
                * std::sqrt(h / (b - a));

            double err = dy.norm();

            if (err <= tol) {
                x += h;
                y = yh;

                xlist.push_back(x);
                ylist.push_back(y);
            }

            if (err > 0.0) {
                double factor =
                    0.95 * std::pow(tol / err, 0.25);

                h *= std::min(factor, 2.0);
            }
            else {
                h *= 2.0;
            }

            if (x + h == x) {
                throw std::runtime_error(
                    "driver: step size became too small"
                );
            }
        }

        return {xlist, ylist};
    }

}