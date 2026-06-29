#pragma once
#include "vector.h"
#include <tuple>
#include <vector>
//Hydrogen.h is mostly just a header function taken from eigenvalues homework files,  with changes here and there to support this homework.

namespace pp {


std::tuple<std::vector<double>, std::vector<vec>>
hydrogen_wavefunction(
    double energy,
    double rmin,
    double rmax,
    double ode_acc,
    double ode_eps,
    double initial_step = 0.01
);



double hydrogen_mismatch(
    double energy,
    double rmin,
    double rmax,
    double ode_acc,
    double ode_eps,
    double initial_step = 0.01
);



double hydrogen_ground_state_energy(
    double rmin,
    double rmax,
    double ode_acc,
    double ode_eps,
    double root_acc = 1e-8,
    double energy_guess = -0.7,
    double energy_step = 1e-4,
    double alpha_min = 1e-3,
    int max_iterations = 100
);

}