#include "minimization.h"
#include <cmath>
#include <exception>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>


namespace {
void print_vector(const pp::vec& x)
{
    std::cout << "(";

    for (int i = 0; i < x.size(); ++i) {
        std::cout << x[i];

        if (i + 1 < x.size()) {
            std::cout << ", ";
        }
    }

    std::cout << ")";
}


void run_minimization(
    const std::string& name,
    const pp::objective_function& phi,
    const pp::vec& start,
    double acc,
    int max_steps
)
{
    const pp::minimization_result result =
        pp::newton_minimize(
            phi,
            start,
            acc,
            max_steps
        );

    const pp::vec final_gradient =
        pp::gradient(phi, result.x);

    std::cout << name << '\n';
    std::cout << std::string(name.size(), '-') << '\n';

    std::cout << "Starting point   = ";
    print_vector(start);
    std::cout << '\n';

    std::cout << "Minimum point    = ";
    print_vector(result.x);
    std::cout << '\n';

    std::cout << "Function value   = "
              << result.value
              << '\n';

    std::cout << "Gradient norm    = "
              << final_gradient.norm()
              << '\n';

    std::cout << "Newton steps     = "
              << result.steps
              << '\n';

    std::cout << "Converged        = "
              << std::boolalpha
              << result.converged
              << "\n\n";
}

} // namespace


int main()
{
    try {
        std::cout << std::setprecision(12);

        const double acc = 1e-6;
        const int max_steps = 1000;

        const pp::objective_function rosenbrock =
            [](const pp::vec& x) {
                const double x0 = x[0];
                const double x1 = x[1];

                return std::pow(1.0 - x0, 2)
                     + 100.0
                     * std::pow(x1 - x0 * x0, 2);
            };

        std::cout << "A. Newton minimization with numerical derivatives\n";
        std::cout << "------------------------------------------------" << std::endl;

        run_minimization(
            "Rosenbrock function",
            rosenbrock,
            pp::vec{-1.2, 1.0},
            acc,
            max_steps
        );


        const pp::objective_function himmelblau =
            [](const pp::vec& x) {
                const double x0 = x[0];
                const double x1 = x[1];

                return std::pow(
                           x0 * x0 + x1 - 11.0,
                           2
                       )
                     + std::pow(
                           x0 + x1 * x1 - 7.0,
                           2
                       );
            };

        const std::vector<pp::vec> himmelblau_starts = {
            pp::vec{ 3.5,  1.5},
            pp::vec{-3.0,  3.0},
            pp::vec{-4.0, -3.0},
            pp::vec{ 4.0, -2.0}
        };

        for (
            int i = 0;
            i < static_cast<int>(himmelblau_starts.size());
            ++i
        ) {
            run_minimization(
                "Himmelblau function, start "
                    + std::to_string(i + 1),
                himmelblau,
                himmelblau_starts[i],
                acc,
                max_steps
            );
        }

        return 0;
    }
    catch (const std::exception& error) {
        std::cerr
            << "Error: "
            << error.what()
            << '\n';

        return 1;
    }
}