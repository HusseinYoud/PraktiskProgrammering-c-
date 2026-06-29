#include "hydrogen.h"
#include <limits>
#include "ode.h"
#include "roots.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace pp {

namespace {

void validate_hydrogen_parameters(
    double rmin,
    double rmax,
    double ode_acc,
    double ode_eps,
    double initial_step
)
{

    if (rmin <= 0.0) {
        throw std::invalid_argument(
            "hydrogen: rmin must be positive"
        );
    }

    if (rmax <= rmin) {
        throw std::invalid_argument(
            "hydrogen: rmax must be larger than rmin"
        );
    }

    if (ode_acc < 0.0 || ode_eps < 0.0) {
        throw std::invalid_argument(
            "hydrogen: ode_acc and ode_eps "
            "must be non-negative"
        );
    }

    if (ode_acc == 0.0 && ode_eps == 0.0) {
        throw std::invalid_argument(
            "hydrogen: ode_acc and ode_eps "
            "cannot both be zero"
        );
    }

    if (initial_step <= 0.0) {
        throw std::invalid_argument(
            "hydrogen: initial_step must be positive"
        );
    }
}



ode_function radial_schrodinger_equation(
    double energy
)
{
    return [energy](
        double r,
        const vec& y
    ) {
        if (y.size() != 2) {
            throw std::invalid_argument(
                "radial_schrodinger_equation: "
                "state vector must have size 2"
            );
        }

        if (r <= 0.0) {
            throw std::invalid_argument(
                "radial_schrodinger_equation: "
                "r must be positive"
            );
        }

        const double f = y[0];
        const double df = y[1];

        const double ddf =
            -2.0
            * (energy + 1.0 / r)
            * f;

        return vec{
            df,
            ddf
        };
    };
}



vec hydrogen_initial_condition(
    double rmin
)
{
    const double f =
        rmin - rmin * rmin;

    const double df =
        1.0 - 2.0 * rmin;

    return vec{
        f,
        df
    };
}

}


std::tuple<std::vector<double>, std::vector<vec>>
hydrogen_wavefunction(
    double energy,
    double rmin,
    double rmax,
    double ode_acc,
    double ode_eps,
    double initial_step
)
{
    validate_hydrogen_parameters(
        rmin,
        rmax,
        ode_acc,
        ode_eps,
        initial_step
    );

    ode_function equation =
        radial_schrodinger_equation(energy);

    vec yinit =
        hydrogen_initial_condition(rmin);

    const double h =
        std::min(
            initial_step,
            rmax - rmin
        );


    return driver(
        equation,
        rmin,
        yinit,
        rmax,
        h,
        ode_acc,
        ode_eps
    );
}


double hydrogen_mismatch(
    double energy,
    double rmin,
    double rmax,
    double ode_acc,
    double ode_eps,
    double initial_step
)
{
    auto [r_values, y_values] =
        hydrogen_wavefunction(
            energy,
            rmin,
            rmax,
            ode_acc,
            ode_eps,
            initial_step
        );

    if (r_values.empty() || y_values.empty()) {
        throw std::runtime_error(
            "hydrogen_mismatch: ODE driver "
            "returned an empty solution"
        );
    }

    const vec& final_state =
        y_values.back();

    if (final_state.size() != 2) {
        throw std::runtime_error(
            "hydrogen_mismatch: invalid final state"
        );
    }


    return final_state[0];
}


double hydrogen_ground_state_energy(
    double rmin,
    double rmax,
    double ode_acc,
    double ode_eps,
    double root_acc,
    double energy_guess,
    double energy_step,
    double alpha_min,
    int max_iterations
)
{
    validate_hydrogen_parameters(
        rmin,
        rmax,
        ode_acc,
        ode_eps,
        0.01
    );

    if (root_acc <= 0.0) {
        throw std::invalid_argument(
            "hydrogen_ground_state_energy: "
            "root_acc must be positive"
        );
    }

    if (energy_step <= 0.0) {
        throw std::invalid_argument(
            "hydrogen_ground_state_energy: "
            "energy_step must be positive"
        );
    }

    if (energy_guess >= 0.0) {
        throw std::invalid_argument(
            "hydrogen_ground_state_energy: "
            "the initial energy guess should be negative"
        );
    }


    root_function shooting_function =
        [=](const vec& energy_vector) {
            if (energy_vector.size() != 1) {
                throw std::invalid_argument(
                    "hydrogen shooting function: "
                    "energy vector must have size 1"
                );
            }

            const double energy =
                energy_vector[0];

            const double mismatch =
                hydrogen_mismatch(
                    energy,
                    rmin,
                    rmax,
                    ode_acc,
                    ode_eps
                );

            return vec{
                mismatch
            };
        };

    const vec initial_energy{
        energy_guess
    };


    const vec dE{
        energy_step
    };

    vec root =
        newton(
            shooting_function,
            initial_energy,
            root_acc,
            alpha_min,
            max_iterations,
            dE
        );

    if (root.size() != 1) {
        throw std::runtime_error(
            "hydrogen_ground_state_energy: "
            "root finder returned an invalid result"
        );
    }

    return root[0];
}


double hydrogen_exact_wavefunction(
    double r
)
{
    if (r < 0.0) {
        throw std::invalid_argument(
            "hydrogen_exact_wavefunction: "
            "r must be non-negative"
        );
    }

    return r * std::exp(-r);
}


double hydrogen_exact_ground_state_energy()
{
    return -0.5;
}

}