#include "ode.h"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numbers>
#include <stdexcept>
#include <string>

namespace {

constexpr double pi = std::numbers::pi_v<double>;


/*
 * Del A:
 * Test af ODE-integratoren på
 *
 *     u'' = -u
 *
 * med u(0) = 0 og u'(0) = 1.
 *
 * Den eksakte løsning er
 *
 *     u(x)  = sin(x)
 *     u'(x) = cos(x).
 */
void write_harmonic_oscillator()
{
    pp::ode_function f =
        [](double x, const pp::vec& y) {
            static_cast<void>(x);

            pp::vec dydx(2);

            dydx[0] = y[1];
            dydx[1] = -y[0];

            return dydx;
        };

    pp::vec yinit = {0.0, 1.0};

    auto [xlist, ylist] =
        pp::driver(
            f,
            0.0,
            yinit,
            10.0,
            0.05,
            1e-6,
            1e-6,
            0.05
        );

    std::ofstream out("oscillator.data");

    if (!out) {
        throw std::runtime_error(
            "Could not open oscillator.data"
        );
    }

    out << std::setprecision(15);

    double max_error = 0.0;

    for (std::size_t i = 0; i < xlist.size(); ++i) {
        const double x = xlist[i];

        const double numerical_u = ylist[i][0];
        const double numerical_du = ylist[i][1];

        const double exact_u = std::sin(x);
        const double exact_du = std::cos(x);

        const double error_u =
            std::abs(numerical_u - exact_u);

        const double error_du =
            std::abs(numerical_du - exact_du);

        max_error = std::max(
            max_error,
            std::max(error_u, error_du)
        );

        out
            << x << ' '
            << numerical_u << ' '
            << exact_u << ' '
            << numerical_du << ' '
            << exact_du << '\n';
    }

    std::cout
        << "A. Harmonic oscillator\n"
        << "Number of accepted points: "
        << xlist.size() << '\n'
        << "Maximum absolute error: "
        << max_error << "\n\n";
}


/*
 * Del B:
 *
 *     u''(phi) + u(phi) = 1 + epsilon*u(phi)^2
 *
 * Med
 *
 *     y[0] = u
 *     y[1] = u'
 *
 * får vi systemet
 *
 *     y[0]' = y[1]
 *     y[1]' = 1 - y[0] + epsilon*y[0]^2.
 */
void write_planetary_orbit(
    const std::string& filename,
    double epsilon,
    const pp::vec& yinit,
    double phi_end
)
{
    pp::ode_function orbit =
        [epsilon](double phi, const pp::vec& y) {
            static_cast<void>(phi);

            pp::vec dydphi(2);

            dydphi[0] = y[1];

            dydphi[1] =
                1.0
                - y[0]
                + epsilon * y[0] * y[0];

            return dydphi;
        };

    auto [philist, ylist] =
        pp::driver(
            orbit,
            0.0,
            yinit,
            phi_end,
            0.02,
            1e-6,
            1e-6,
            0.02
        );

    std::ofstream out(filename);

    if (!out) {
        throw std::runtime_error(
            "Could not open " + filename
        );
    }

    out << std::setprecision(15);

    /*
     * Kolonne 1: phi
     * Kolonne 2: u(phi) = 1/r(phi)
     *
     * Banen kan tegnes i gnuplot med
     *
     * x = cos(phi)/u
     * y = sin(phi)/u.
     */
    for (std::size_t i = 0; i < philist.size(); ++i) {
        const double phi = philist[i];
        const double u = ylist[i][0];

        out
            << phi << ' '
            << u << '\n';
    }

    std::cout
        << "Created " << filename << '\n'
        << "epsilon = " << epsilon << '\n'
        << "Number of accepted points: "
        << philist.size() << "\n\n";
}


/*
 * Del C:
 *
 * Tilstandsvektoren er
 *
 * z = {
 *     x1, y1, vx1, vy1,
 *     x2, y2, vx2, vy2,
 *     x3, y3, vx3, vy3
 * }.
 */
pp::vec three_body_rhs(
    double t,
    const pp::vec& z
)
{
    static_cast<void>(t);

    constexpr int number_of_bodies = 3;
    constexpr int components_per_body = 4;

    pp::vec dzdt(12);

    for (int i = 0; i < number_of_bodies; ++i) {
        const int i0 = components_per_body * i;

        const double xi = z[i0];
        const double yi = z[i0 + 1];

        /*
         * Positionsafledte er hastighederne.
         */
        dzdt[i0] = z[i0 + 2];
        dzdt[i0 + 1] = z[i0 + 3];

        double ax = 0.0;
        double ay = 0.0;

        for (int j = 0; j < number_of_bodies; ++j) {
            if (j == i) {
                continue;
            }

            const int j0 = components_per_body * j;

            const double xj = z[j0];
            const double yj = z[j0 + 1];

            const double dx = xj - xi;
            const double dy = yj - yi;

            const double distance_squared =
                dx * dx + dy * dy;

            if (distance_squared == 0.0) {
                throw std::runtime_error(
                    "Collision in three-body system"
                );
            }

            /*
             * 1/r^3 = 1/(r^2 * sqrt(r^2)).
             */
            const double inverse_distance_cubed =
                1.0
                / (
                    distance_squared
                    * std::sqrt(distance_squared)
                );

            ax += dx * inverse_distance_cubed;
            ay += dy * inverse_distance_cubed;
        }

        dzdt[i0 + 2] = ax;
        dzdt[i0 + 3] = ay;
    }

    return dzdt;
}


void write_three_body_orbit()
{
    /*
     * Initialbetingelser for figure-eight-løsningen.
     *
     * Rækkefølge for hvert legeme:
     *
     * x, y, vx, vy
     */
    pp::vec zinit = {
         0.97000436,
        -0.24308753,
         0.466203685,
         0.432365730,

        -0.97000436,
         0.24308753,
         0.466203685,
         0.432365730,

         0.0,
         0.0,
        -0.932407370,
        -0.864731460
    };

    constexpr double period = 6.32591398;

    auto [tlist, zlist] =
        pp::driver(
            three_body_rhs,
            0.0,
            zinit,
            period,
            0.01,
            1e-7,
            1e-7,
            0.01
        );

    std::ofstream out("threebody.data");

    if (!out) {
        throw std::runtime_error(
            "Could not open threebody.data"
        );
    }

    out << std::setprecision(15);

    /*
     * Kolonner:
     *
     * 1: t
     * 2: x1
     * 3: y1
     * 4: x2
     * 5: y2
     * 6: x3
     * 7: y3
     */
    for (std::size_t i = 0; i < tlist.size(); ++i) {
        const pp::vec& z = zlist[i];

        out
            << tlist[i] << ' '
            << z[0] << ' '
            << z[1] << ' '
            << z[4] << ' '
            << z[5] << ' '
            << z[8] << ' '
            << z[9] << '\n';
    }

    const pp::vec final_difference =
        zlist.back() - zinit;

    std::cout
        << "C. Three-body figure-eight orbit\n"
        << "Number of accepted points: "
        << tlist.size() << '\n'
        << "Difference between initial and final state: "
        << final_difference.norm() << "\n\n";
}

} // namespace


int main()
{
    std::cout << std::setprecision(12);

    /*
     * Del A
     */
    write_harmonic_oscillator();

    /*
     * Del B1:
     * Newtonsk cirkel.
     *
     * epsilon = 0
     * u(0)     = 1
     * u'(0)    = 0
     */
    write_planetary_orbit(
        "orbit-circular.data",
        0.0,
        pp::vec{1.0, 0.0},
        12.0 * pi
    );

    /*
     * Del B2:
     * Newtonsk ellipse.
     *
     * epsilon = 0
     * u(0)     = 1
     * u'(0)    = -0.5
     */
    write_planetary_orbit(
        "orbit-elliptic.data",
        0.0,
        pp::vec{1.0, -0.5},
        12.0 * pi
    );

    /*
     * Del B3:
     * Relativistisk præcession.
     *
     * epsilon = 0.01
     * u(0)     = 1
     * u'(0)    = -0.5
     */
    write_planetary_orbit(
        "orbit-relativistic.data",
        0.01,
        pp::vec{1.0, -0.5},
        20.0 * pi
    );

    /*
     * Del C
     */
    write_three_body_orbit();

    return 0;
}