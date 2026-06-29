#include "integration.h"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numbers>
#include <string>
namespace {

void test_integral(
    const std::string& name,
    const pp::integrand& f,
    double a,
    double b,
    double exact,
    double acc,
    double eps
)
{
    std::size_t calls = 0;

  
    const pp::integrand counted_function =
        [&f, &calls](double x) {
            ++calls;
            return f(x);
        };

    const double result =
        pp::integrate(
            counted_function,
            a,
            b,
            acc,
            eps
        );

    const double absolute_error =
        std::abs(result - exact);


    const double accuracy_goal =
        acc + eps * std::abs(exact);

    const bool passed =
        absolute_error <= accuracy_goal;

    std::cout << name << '\n';
    std::cout << "  numerical result : "
              << result << '\n';

    std::cout << "  exact result     : "
              << exact << '\n';

    std::cout << "  absolute error   : "
              << absolute_error << '\n';

    std::cout << "  accuracy goal    : "
              << accuracy_goal << '\n';

    std::cout << "  function calls   : "
              << calls << '\n';

    std::cout << "  test             : "
              << (passed ? "PASSED" : "FAILED")
              << "\n\n";
}
void compare_integrators(
    const std::string& name,
    const pp::integrand& f,
    double a,
    double b,
    double exact,
    double acc,
    double eps
)
{
    std::size_t ordinary_calls = 0;
    std::size_t cc_calls = 0;

    const pp::integrand ordinary_function =
        [&f, &ordinary_calls](double x) {
            ++ordinary_calls;
            return f(x);
        };

    const pp::integrand cc_function =
        [&f, &cc_calls](double x) {
            ++cc_calls;
            return f(x);
        };

    const pp::integration_result ordinary =
        pp::integrate_with_error(
            ordinary_function,
            a,
            b,
            acc,
            eps
        );

    const pp::integration_result cc =
        pp::clenshaw_curtis_with_error(
            cc_function,
            a,
            b,
            acc,
            eps
        );

    std::cout << name << '\n';

    std::cout
        << std::setw(20) << "method"
        << std::setw(20) << "result"
        << std::setw(20) << "actual error"
        << std::setw(20) << "estimated error"
        << std::setw(16) << "evaluations"
        << '\n';

    std::cout
        << std::setw(20) << "ordinary"
        << std::setw(20) << ordinary.value
        << std::setw(20)
        << std::abs(ordinary.value - exact)
        << std::setw(20) << ordinary.error
        << std::setw(16) << ordinary_calls
        << '\n';

    std::cout
        << std::setw(20) << "Clenshaw-Curtis"
        << std::setw(20) << cc.value
        << std::setw(20)
        << std::abs(cc.value - exact)
        << std::setw(20) << cc.error
        << std::setw(16) << cc_calls
        << "\n\n";
}
void test_infinite_integral(
    const std::string& name,
    const pp::integrand& f,
    double a,
    double b,
    double exact,
    double acc,
    double eps
)
{
    std::size_t calls = 0;

    const pp::integrand counted_function =
        [&f, &calls](double x) {
            ++calls;
            return f(x);
        };

    const pp::integration_result result =
        pp::clenshaw_curtis_with_error(
            counted_function,
            a,
            b,
            acc,
            eps
        );

    std::cout
        << name << '\n'
        << "  numerical result : "
        << result.value << '\n'
        << "  exact result     : "
        << exact << '\n'
        << "  actual error     : "
        << std::abs(result.value - exact) << '\n'
        << "  estimated error  : "
        << result.error << '\n'
        << "  function calls   : "
        << calls << "\n\n";
}
} // namespace


int main()
{
    constexpr double pi = std::numbers::pi;
    const double infinity = std::numeric_limits<double>::infinity();
    const double acc = 1e-6;
    const double eps = 1e-6;

    std::cout << std::scientific
              << std::setprecision(12);

    std::cout
        << "A. Recursive adaptive integration\n"
        << "=================================\n\n";

    test_integral(
        "1) integral_0^1 sqrt(x) dx",
        [](double x) {
            return std::sqrt(x);
        },
        0.0,
        1.0,
        2.0 / 3.0,
        acc,
        eps
    );


    test_integral(
        "2) integral_0^1 1/sqrt(x) dx",
        [](double x) {
            return 1.0 / std::sqrt(x);
        },
        0.0,
        1.0,
        2.0,
        acc,
        eps
    );


    test_integral(
        "3) integral_-1^1 sqrt(1-x^2) dx",
        [](double x) {
            return std::sqrt(1.0 - x * x);
        },
        -1.0,
        1.0,
        pi / 2.0,
        acc,
        eps
    );

    test_integral(
        "4) integral_0^1 log(x)/sqrt(x) dx",
        [](double x) {
            return std::log(x) / std::sqrt(x);
        },
        0.0,
        1.0,
        -4.0,
        acc,
        eps
    );

    std::cout
        << "A2. Error function\n"
        << "=================\n\n";

    std::ofstream erf_data("erf.data");

    if (!erf_data) {
        std::cerr
            << "Could not open erf.data\n";

        return 1;
    }

    erf_data << std::scientific
             << std::setprecision(15);

    const double erf_acc = 1e-7;
    const double erf_eps = 1e-7;

    for (int i = -30; i <= 30; ++i) {
        const double z =
            static_cast<double>(i) / 10.0;

        const double numerical =
            pp::erf(z, erf_acc, erf_eps);

        const double reference =
            std::erf(z);

        const double difference =
            std::abs(numerical - reference);

        erf_data
            << z << ' '
            << numerical << ' '
            << reference << ' '
            << difference << '\n';
    }

    erf_data.close();

    const double numerical_erf1 =
        pp::erf(1.0, erf_acc, erf_eps);

    const double exact_erf1 =
        0.84270079294971486934;

    std::cout
        << "erf(1), numerical : "
        << numerical_erf1 << '\n';

    std::cout
        << "erf(1), reference : "
        << exact_erf1 << '\n';

    std::cout
        << "absolute error    : "
        << std::abs(numerical_erf1 - exact_erf1)
        << "\n\n";


    std::cout
        << "A3. Convergence of erf(1)\n"
        << "========================\n\n";

    std::ofstream convergence_data(
        "erf-convergence.data"
    );

    if (!convergence_data) {
        std::cerr
            << "Could not open erf-convergence.data\n";

        return 1;
    }

    convergence_data
        << std::scientific
        << std::setprecision(15);

    std::cout
        << std::setw(16) << "acc"
        << std::setw(22) << "erf(1)"
        << std::setw(22) << "absolute error"
        << '\n';


    for (int k = 1; k <= 10; ++k) {
        const double current_acc =
            std::pow(10.0, -k);

        const double current_eps = 0.0;

        const double numerical =
            pp::erf(
                1.0,
                current_acc,
                current_eps
            );

        const double absolute_error =
            std::abs(numerical - exact_erf1);


        convergence_data
            << current_acc << ' '
            << absolute_error << ' '
            << numerical << '\n';

        std::cout
            << std::setw(16) << current_acc
            << std::setw(22) << numerical
            << std::setw(22) << absolute_error
            << '\n';
    }

    convergence_data.close();
    std::cout
        << "\n\nB. Variable transformation quadratures\n"
        << "======================================\n\n";
    compare_integrators(
        "B1) integral_0^1 dx/sqrt(x)",
        [](double x) {
            return 1.0 / std::sqrt(x);
        },
        0.0,
        1.0,
        2.0,
        acc,
        eps
    );

    compare_integrators(
        "B2) integral_0^1 log(x)/sqrt(x) dx",
        [](double x) {
            return std::log(x) / std::sqrt(x);
        },
        0.0,
        1.0,
        -4.0,
        acc,
        eps
    );
    std::cout
    << "B3. Infinite integration limits\n"
    << "-------------------------------\n\n";
    test_infinite_integral(
    "integral_0^infinity exp(-x) dx",
    [](double x) {
        return std::exp(-x);
    },
    0.0,
    infinity,
    1.0,
    acc,
    eps
);
test_infinite_integral(
    "integral_-infinity^infinity exp(-x^2) dx",
    [](double x) {
        return std::exp(-x * x);
    },
    -infinity,
    infinity,
    std::sqrt(pi),
    acc,
    eps
);
test_infinite_integral(
    "integral_0^infinity dx/(1+x^2)",
    [](double x) {
        return 1.0 / (1.0 + x * x);
    },
    0.0,
    infinity,
    pi / 2.0,
    acc,
    eps
);
std::cout
    << "C. Error estimate\n"
    << "=================\n\n";

std::ofstream error_data(
    "error-estimate.data"
);

if (!error_data) {
    std::cerr
        << "Could not open error-estimate.data\n";

    return 1;
}

error_data
    << std::scientific
    << std::setprecision(15);

std::cout
    << std::setw(16) << "acc"
    << std::setw(22) << "actual error"
    << std::setw(22) << "estimated error"
    << std::setw(16) << "evaluations"
    << '\n';

for (int k = 2; k <= 9; ++k) {
    const double current_acc =
        std::pow(10.0, -k);

    std::size_t calls = 0;

    const pp::integrand difficult =
        [&calls](double x) {
            ++calls;

            return std::log(x)
                   /
                   std::sqrt(x);
        };

    const pp::integration_result result =
        pp::clenshaw_curtis_with_error(
            difficult,
            0.0,
            1.0,
            current_acc,
            0.0
        );

    const double actual_error =
        std::abs(result.value + 4.0);

    error_data
        << current_acc << ' '
        << actual_error << ' '
        << result.error << ' '
        << calls << '\n';

    std::cout
        << std::setw(16) << current_acc
        << std::setw(22) << actual_error
        << std::setw(22) << result.error
        << std::setw(16) << calls
        << '\n';
}

error_data.close();
    return 0;
}