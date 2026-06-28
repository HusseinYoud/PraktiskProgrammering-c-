#include "ode.h"

#include <cmath>
#include <iomanip>
#include <iostream>

int main()
{
    pp::ode_function f =
        [](double x, const pp::vec& y) {
            static_cast<void>(x);

            pp::vec dydx(2);

            dydx[0] = y[1];
            dydx[1] = -y[0];

            return dydx;
        };

    double a = 0.0;
    double b = 10.0;

    pp::vec yinit = {0.0, 1.0};

    auto [xlist, ylist] =
        pp::driver(
            f,
            a,
            yinit,
            b,
            0.125,
            1e-5,
            1e-5
        );

    std::cout
        << "# x"
        << " numerical_u"
        << " exact_u"
        << " numerical_du"
        << " exact_du"
        << '\n';

    std::cout << std::setprecision(12);

    for (std::size_t i = 0; i < xlist.size(); ++i) {
        double x = xlist[i];

        std::cout
            << x << ' '
            << ylist[i][0] << ' '
            << std::sin(x) << ' '
            << ylist[i][1] << ' '
            << std::cos(x) << '\n';
    }

    return 0;
}