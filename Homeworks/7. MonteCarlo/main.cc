#include "mc.h"

#include <cmath>
#include <cstddef>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numbers>
#include <random>
#include <stdexcept>
#include <string_view>
#include <vector>

namespace {

constexpr double pi =
    std::numbers::pi_v<double>;

const std::vector<std::size_t> sample_counts{
    100,
    300,
    1'000,
    3'000,
    10'000,
    30'000,
    100'000,
    300'000,
    1'000'000
};



double unit_circle(const pp::vec& x)
{
    const double radius_squared =
        x[0] * x[0]
        + x[1] * x[1];

    return radius_squared <= 1.0
        ? 1.0
        : 0.0;
}


double ellipsoid(const pp::vec& x)
{
    constexpr double axis_a = 1.0;
    constexpr double axis_b = 2.0;
    constexpr double axis_c = 3.0;

    const double radius_squared =
        x[0] * x[0] / (axis_a * axis_a)
        + x[1] * x[1] / (axis_b * axis_b)
        + x[2] * x[2] / (axis_c * axis_c);

    return radius_squared <= 1.0
        ? 1.0
        : 0.0;
}


double gaussian_2d(const pp::vec& x)
{
    return std::exp(
        -x[0] * x[0]
        -x[1] * x[1]
    );
}



double difficult_integrand(const pp::vec& x)
{
    const double denominator =
        1.0
        - std::cos(pi * x[0])
        * std::cos(pi * x[1])
        * std::cos(pi * x[2]);

    return 1.0 / denominator;
}


double gaussian_exact()
{
    const double one_dimensional =
        0.5 * std::sqrt(pi) * std::erf(1.0);

    return one_dimensional * one_dimensional;
}


void print_result(
    std::string_view method,
    const pp::mc_result& result,
    double exact
)
{
    const double actual_error =
        std::abs(result.integral - exact);

    std::cout
        << std::left
        << std::setw(22)
        << method

        << std::right
        << " result = "
        << std::setw(20)
        << result.integral

        << "  estimated error = "
        << std::setw(14)
        << result.error

        << "  actual error = "
        << actual_error
        << '\n';
}



void write_circle_errors()
{
    std::ofstream data("out.circle.txt");

    if (!data) {
        throw std::runtime_error(
            "Could not open out.circle.txt"
        );
    }

    const pp::vec lower{-1.0, -1.0};
    const pp::vec upper{ 1.0,  1.0};

    data << std::scientific
         << std::setprecision(15);

    data
        << "# N estimated_error actual_error "
        << "reference_1_over_sqrt_N integral\n";

    for (const std::size_t N : sample_counts) {
        pp::lcg generator(17);

        auto random_double =
            [&generator]()
        {
            return generator.next();
        };

        const pp::mc_result result =
            pp::plainmc(
                unit_circle,
                lower,
                upper,
                N,
                random_double
            );

        const double actual_error =
            std::abs(result.integral - pi);


        const double reference_error =
            std::sqrt(4.0 * pi - pi * pi)
            / std::sqrt(static_cast<double>(N));

        data
            << N << ' '
            << result.error << ' '
            << actual_error << ' '
            << reference_error << ' '
            << result.integral << '\n';
    }
}



void write_quasi_errors()
{
    std::ofstream data("out.quasi.txt");

    if (!data) {
        throw std::runtime_error(
            "Could not open out.quasi.txt"
        );
    }

    const pp::vec lower{0.0, 0.0};
    const pp::vec upper{1.0, 1.0};

    const double exact =
        gaussian_exact();

    const pp::Halton sequence1(2, 0);
    const pp::Halton sequence2(2, 2);

    data << std::scientific
         << std::setprecision(15);

    data
        << "# N "
        << "lcg_actual lcg_estimated "
        << "std_actual std_estimated "
        << "quasi_actual quasi_estimated "
        << "reference_sqrt reference_inverse\n";

    for (const std::size_t N : sample_counts) {


        pp::lcg lcg_generator(19);

        auto lcg_random =
            [&lcg_generator]()
        {
            return lcg_generator.next();
        };

        const pp::mc_result lcg_result =
            pp::plainmc(
                gaussian_2d,
                lower,
                upper,
                N,
                lcg_random
            );



        std::mt19937_64 engine(19);

        std::uniform_real_distribution<double> uniform(
            0.0,
            1.0
        );

        auto standard_random =
            [&engine, &uniform]()
        {
            return uniform(engine);
        };

        const pp::mc_result standard_result =
            pp::plainmc(
                gaussian_2d,
                lower,
                upper,
                N,
                standard_random
            );



        const std::size_t points_per_sequence =
            std::max<std::size_t>(1, N / 2);

        const pp::mc_result quasi_result =
            pp::quasimc(
                gaussian_2d,
                lower,
                upper,
                points_per_sequence,
                sequence1,
                sequence2
            );

        data
            << N << ' '

            << std::abs(
                lcg_result.integral - exact
            ) << ' '

            << lcg_result.error << ' '

            << std::abs(
                standard_result.integral - exact
            ) << ' '

            << standard_result.error << ' '

            << std::abs(
                quasi_result.integral - exact
            ) << ' '

            << quasi_result.error << ' '

            << 1.0
                / std::sqrt(
                    static_cast<double>(N)
                ) << ' '

            << 1.0
                / static_cast<double>(N)
            << '\n';
    }
}

} // namespace


int main()
{
    try {
        std::cout
            << std::scientific
            << std::setprecision(15);



        std::cout
            << "A. Plain Monte Carlo integration\n"
            << "============================================\n\n";

        constexpr std::size_t N_part_a =
            1'000'000;



        {
            const pp::vec lower{-1.0, -1.0};
            const pp::vec upper{ 1.0,  1.0};

            pp::lcg generator(1);

            auto random_double =
                [&generator]()
            {
                return generator.next();
            };

            const pp::mc_result result =
                pp::plainmc(
                    unit_circle,
                    lower,
                    upper,
                    N_part_a,
                    random_double
                );

            std::cout
                << "Area of the unit circle\n"
                << "N = " << N_part_a << '\n';

            print_result(
                "LCG plain MC",
                result,
                pi
            );

            std::cout << '\n';
        }



        {
            const pp::vec lower{-1.0, -2.0, -3.0};
            const pp::vec upper{ 1.0,  2.0,  3.0};

            constexpr double exact_volume =
                (4.0 / 3.0) * pi * 1.0 * 2.0 * 3.0;

            pp::lcg generator(2);

            auto random_double =
                [&generator]()
            {
                return generator.next();
            };

            const pp::mc_result result =
                pp::plainmc(
                    ellipsoid,
                    lower,
                    upper,
                    N_part_a,
                    random_double
                );

            std::cout
                << "Volume of the ellipsoid\n"
                << "N = " << N_part_a << '\n';

            print_result(
                "LCG plain MC",
                result,
                exact_volume
            );

            std::cout << '\n';
        }



        std::cout
            << "B. Quasi-random sequences\n"
            << "============================================\n\n";

        {
            const pp::vec lower{0.0, 0.0, 0.0};
            const pp::vec upper{1.0, 1.0, 1.0};

            constexpr double exact =
                1.3932039296856768591842462603255;


            constexpr std::size_t N =
                1'000'000;


            pp::lcg lcg_generator(3);

            auto lcg_random =
                [&lcg_generator]()
            {
                return lcg_generator.next();
            };

            const pp::mc_result lcg_result =
                pp::plainmc(
                    difficult_integrand,
                    lower,
                    upper,
                    N,
                    lcg_random
                );



            std::mt19937_64 engine(3);

            std::uniform_real_distribution<double> uniform(
                0.0,
                1.0
            );

            auto standard_random =
                [&engine, &uniform]()
            {
                return uniform(engine);
            };

            const pp::mc_result standard_result =
                pp::plainmc(
                    difficult_integrand,
                    lower,
                    upper,
                    N,
                    standard_random
                );



            const pp::Halton sequence1(3, 0);
            const pp::Halton sequence2(3, 3);

            const pp::mc_result quasi_result =
                pp::quasimc(
                    difficult_integrand,
                    lower,
                    upper,
                    N / 2,
                    sequence1,
                    sequence2
                );

            std::cout
                << "Difficult three-dimensional integral\n"
                << "Total evaluation budget = " << N << '\n'
                << "Exact result = " << exact << "\n\n";

            print_result(
                "LCG",
                lcg_result,
                exact
            );

            print_result(
                "Standard C++",
                standard_result,
                exact
            );

            print_result(
                "Halton",
                quasi_result,
                exact
            );

            std::cout
                << "\nThe difficult integral is not used "
                << "for the convergence test.\n\n";
        }



        std::cout
            << "C. Recursive stratified sampling\n"
            << "============================================\n\n";

        {
            const pp::vec lower{-1.0, -1.0};
            const pp::vec upper{ 1.0,  1.0};

            constexpr std::size_t N =
                200'000;

            constexpr std::size_t nmin =
                64;


            pp::lcg plain_generator(4);

            auto plain_random =
                [&plain_generator]()
            {
                return plain_generator.next();
            };

            const pp::mc_result plain_result =
                pp::plainmc(
                    unit_circle,
                    lower,
                    upper,
                    N,
                    plain_random
                );


    
            pp::lcg stratified_generator(4);

            auto stratified_random =
                [&stratified_generator]()
            {
                return stratified_generator.next();
            };

            const pp::mc_result stratified_result =
                pp::stratifiedmc(
                    unit_circle,
                    lower,
                    upper,
                    N,
                    stratified_random,
                    nmin
                );

            std::cout
                << "Unit-circle comparison\n"
                << "N = " << N
                << ", nmin = " << nmin
                << "\n\n";

            print_result(
                "Plain MC",
                plain_result,
                pi
            );

            print_result(
                "Stratified MC",
                stratified_result,
                pi
            );

            std::cout << '\n';
        }


  
        write_circle_errors();
        write_quasi_errors();

        std::cout
            << "Generated data files:\n"
            << "  out.circle.txt\n"
            << "  out.quasi.txt\n";

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