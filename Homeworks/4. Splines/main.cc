#include "Spliner.h"

#include <cmath>
#include <iomanip>
#include <iostream>

int main()
{
    const double x_min = 0.0;
    const double x_max = 9.0;
    const double dx_data = 0.5;

    const int n_data =
        static_cast<int>((x_max - x_min) / dx_data) + 1;

    pp::vec x(n_data);
    pp::vec y(n_data);

    for (int i = 0; i < n_data; ++i) {
        x[i] = x_min + i * dx_data;
        y[i] = std::cos(x[i]);
    }

    pp::qspline qcos(x, y);

    std::cerr << std::setprecision(12);



    const double z_test = 0.5;

    const double spline_integral =
        pp::linterp_intg(x, y, z_test);

    const double exact_integral =
        std::sin(z_test);

    std::cerr
        << "A. Linear spline (linear interpolation)\n"
        << "----------------------------------------\n";

    std::cerr << "Integral of cos from 0 to 0.5\n";
    std::cerr << "Linear spline: "
              << spline_integral << '\n';

    std::cerr << "Exact cos-integral: "
              << exact_integral << '\n';




    std::cerr
        << "\nB. Quadratic spline\n"
        << "----------------------------------------\n";

    pp::vec x_test{
        1.0, 2.0, 3.0, 4.0, 5.0
    };

    pp::vec y_constant{
        1.0, 1.0, 1.0, 1.0, 1.0
    };

    pp::vec y_linear{
        1.0, 2.0, 3.0, 4.0, 5.0
    };

    pp::vec y_quadratic{
        1.0, 4.0, 9.0, 16.0, 25.0
    };

    pp::qspline constant_spline(x_test, y_constant);
    pp::qspline linear_spline(x_test, y_linear);
    pp::qspline quadratic_spline(x_test, y_quadratic);

    std::cerr << "\ny = 1\n";
    std::cerr << "b = " << constant_spline.b << '\n';
    std::cerr << "c = " << constant_spline.c << '\n';

    std::cerr << "\ny = x\n";
    std::cerr << "b = " << linear_spline.b << '\n';
    std::cerr << "c = " << linear_spline.c << '\n';

    std::cerr << "\ny = x^2\n";
    std::cerr << "b = " << quadratic_spline.b << '\n';
    std::cerr << "c = " << quadratic_spline.c << '\n';



    std::cout << std::setprecision(16);

    const int n_plot = 901;


    for (int i = 0; i < n_plot; ++i) {
        const double z =
            x_min
            + (x_max - x_min)
            * static_cast<double>(i)
            / static_cast<double>(n_plot - 1);

        std::cout
            << z << ' '
            << pp::linterp(x, y, z) << ' '
            << pp::linterp_intg(x, y, z)
            << '\n';
    }

    std::cout << "\n\n";


    for (int i = 0; i < n_data; ++i) {
        std::cout
            << x[i] << ' '
            << y[i]
            << '\n';
    }

    std::cout << "\n\n";


    for (int i = 0; i < n_plot; ++i) {
        const double z =
            x_min
            + (x_max - x_min)
            * static_cast<double>(i)
            / static_cast<double>(n_plot - 1);

        std::cout
            << z << ' '
            << qcos.eval(z) << ' '
            << qcos.deriv(z) << ' '
            << qcos.intg(z)
            << '\n';
    }

    return 0;
}