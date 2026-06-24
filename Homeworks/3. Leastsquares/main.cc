#include "lsfit.h"

#include <cmath>
#include <fstream>
#include <iostream>
#include <tuple>

//Large sections modified by ChatGPT to ensure a better out.txt as well as better plots
int main() {

    pp::vec t {
        1, 2, 3, 4, 6, 9, 10, 13, 15
    };

    pp::vec y {
        117, 100, 88, 72, 53,
        29.5, 25.2, 15.2, 11.1
    };

    pp::vec dy {
        6, 5, 4, 4, 4, 3, 3, 2, 2
    };

    pp::vec log_y(y.size());
    pp::vec dlog_y(y.size());

    for (int i = 0; i < y.size(); ++i) {
        log_y[i] = std::log(y[i]);


        dlog_y[i] = dy[i] / y[i];
    }

    pp::Funcs fs {
        [](double) {
            return 1.0;
        },

        [](double z) {
            return z;
        }
    };

    const auto [c, cov] =
        pp::lsfit(fs, t, log_y, dlog_y);

    const double dc0 = std::sqrt(cov(0, 0));
    const double dc1 = std::sqrt(cov(1, 1));

    const double log_a = c[0];
    const double slope = c[1];

    const double a = std::exp(log_a);
    const double lambda = -slope;
    const double half_life = std::log(2.0) / lambda;

    const double modern_half_life = 3.66;

    double chi2 = 0.0;

    for (int i = 0; i < t.size(); ++i) {
        const double fitted_log_y =
            pp::eval(fs, c, t[i]);

        const double residual =
            (fitted_log_y - log_y[i])
            / dlog_y[i];

        chi2 += residual * residual;
    }

    const int dof =
        t.size() - static_cast<int>(fs.size());

    std::cout
        << "A. Ordinary least-squares fit by QR decomposition\n";

    std::cout
        << "------------------------------------------------\n\n";

    std::cout
        << "The linear fit is of the form: ln(y) = c0 + c1 * t\n" << " where c0 = ln(a) and c1 = -lambda\n\n";
    


    std::cout
        << "Error propagation:\n"
        << "For f(y) = ln(y), f'(y) = 1/y.\n"
        << "Therefore delta(ln(y)) is approximately dy/y.\n\n";

    std::cout
        << "c0 = ln(a) = "
        << log_a
        << '\n';

    std::cout
        << "c1 = -lambda = "
        << slope
        << '\n';

    std::cout
        << "a = "
        << a
        << '\n';

    std::cout
        << "lambda = "
        << lambda
        << " 1/day\n";

    std::cout
        << "half-life = "
        << half_life
        << " days\n";

    std::cout
        << "chi^2 = "
        << chi2
        << '\n';

    std::cout
        << "degrees of freedom = "
        << dof
        << '\n';

    const double relative_difference =
        std::abs(half_life - modern_half_life)
        / modern_half_life
        * 100.0;

    std::cout
        << "modern half-life of Ra-224 = "
        << modern_half_life
        << " days\n";

    std::cout
        << "relative difference = "
        << relative_difference
        << " %\n\n";


    std::ofstream data("fit.data");

    if (!data) {
        std::cerr << "Could not open fit.data\n";
        return 1;
    }


    for (int i = 0; i < t.size(); ++i) {
        data
            << t[i] << ' '
            << y[i] << ' '
            << dy[i] << '\n';
    }

    const int number_of_points = 200;
    const double t_min = 0.0;
    const double t_max = 16.0;


    const auto write_curve =
        [&](double curve_c0, double curve_c1) {

            data << "\n\n";

            for (int i = 0; i < number_of_points; ++i) {
                const double z =
                    t_min
                    + (t_max - t_min)
                    * i
                    / static_cast<double>(number_of_points - 1);

                const double fitted_y =
                    std::exp(curve_c0 + curve_c1 * z);

                data
                    << z << ' '
                    << fitted_y << '\n';
            }
        };

 
    write_curve(
        c[0],
        c[1]
    );


    write_curve(
        c[0] + dc0,
        c[1] + dc1
    );


    write_curve(
        c[0] + dc0,
        c[1] - dc1
    );

    write_curve(
        c[0] - dc0,
        c[1] + dc1
    );


    write_curve(
        c[0] - dc0,
        c[1] - dc1
    );

    data.close();

    std::cout
        << "B. Uncertainties of the fitting coefficients\n";

    std::cout
        << "------------------------------------------------\n\n";


    const double da =
        a * dc0;

    const double dlambda =
        dc1;


    const double dhalf_life =
        std::log(2.0)
        / (lambda * lambda)
        * dlambda;

    std::cout
        << "Covariance matrix:\n";

    for (int i = 0; i < cov.size1(); ++i) {
        for (int j = 0; j < cov.size2(); ++j) {
            std::cout
                << cov(i, j)
                << '\t';
        }

        std::cout << '\n';
    }

    std::cout << '\n';

    std::cout
        << "c0 = "
        << c[0]
        << " +/- "
        << dc0
        << '\n';

    std::cout
        << "c1 = "
        << c[1]
        << " +/- "
        << dc1
        << '\n';

    std::cout
        << "a = "
        << a
        << " +/- "
        << da
        << '\n';

    std::cout
        << "lambda = "
        << lambda
        << " +/- "
        << dlambda
        << " 1/day\n";

    std::cout
        << "half-life = "
        << half_life
        << " +/- "
        << dhalf_life
        << " days\n";

    const double lower_bound =
        half_life - dhalf_life;

    const double upper_bound =
        half_life + dhalf_life;

    std::cout
        << "one-sigma interval = ["
        << lower_bound
        << ", "
        << upper_bound
        << "] days\n";

    std::cout
        << "modern value = "
        << modern_half_life
        << " days\n";

    const bool modern_value_is_within_uncertainty =
        modern_half_life >= lower_bound
        && modern_half_life <= upper_bound;

    if (modern_value_is_within_uncertainty) {
        std::cout
            << "The modern value agrees with the fitted value "
            << "within the estimated uncertainty.\n";
    } else {
        std::cout
            << "The modern value does not agree with the fitted value "
            << "within the estimated uncertainty.\n";
    }
    {
        std::cout << "C.Evaliation of the quality of the uncertainties of the fit coefficients" << std::endl;

        std::cout << "------------------------------------------------" << std::endl;

        std::cout << "The best fit as well as the coefficients are = " << a << '\n' <<  "as well as" <<  c[1] + dc1 << std::endl;
        
    }

    return 0;
}