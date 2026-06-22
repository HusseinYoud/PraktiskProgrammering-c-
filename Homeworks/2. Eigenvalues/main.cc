#include "matrix.h"
#include "jacobi.h"
#include <iostream>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include <string>
#include <stdexcept>
#include <vector>


double exact_1s(double r) {
    return 2.0 * r * std::exp(-r);
}

double exact_2s(double r) {
    return r * (2.0 - r)
         * std::exp(-r / 2.0)
         / (2.0 * std::sqrt(2.0));
}

double exact_3s(double r) {
    return 2.0 * r
         * (27.0 - 18.0 * r + 2.0 * r * r)
         * std::exp(-r / 3.0)
         / (81.0 * std::sqrt(3.0));
}

void sort_eigen(pp::vec& w, pp::matrix& V) {
    const int n = w.size();

    for (int i = 0; i < n - 1; ++i) {
        int smallest = i; 

        for(int j = i + 1; j < n; ++j) {
            if (w[j] < w[smallest]) {
                smallest = j;
            }
        }
            if (smallest != i) {
                std::swap(w[i], w[smallest]);
                std::swap(V[i], V[smallest]);
            }
    } 
}

struct HydrogenSolution {
    pp::vec r;
    pp::vec energies;
    pp::matrix V;
    double dr;
};

HydrogenSolution solve_hydrogen(double rmax, double dr) {
    using namespace pp;

    if (rmax <= 0.0 || dr <= 0.0 || rmax <= dr) {
        throw std::invalid_argument(
            "rmax and dr must be positive, and rmax must be larger than dr"
        );
    }

    const int npoints =
        static_cast<int>(rmax / dr) - 1;

    if (npoints < 1) {
        throw std::invalid_argument(
            "Too few radial grid points"
        );
    }

    vec r(npoints);

    for (int i = 0; i < npoints; ++i) {
        r[i] = (i + 1) * dr;
    }

    matrix H(npoints, npoints);

    const double off_diagonal =
        -0.5 / (dr * dr);

    for (int i = 0; i < npoints; ++i) {
        H(i, i) =
            -2.0 * off_diagonal - 1.0 / r[i];
    }

    for (int i = 0; i < npoints - 1; ++i) {
        H(i, i + 1) = off_diagonal;
        H(i + 1, i) = off_diagonal;
    }

    EVD hydrogen(H);

    vec energies = hydrogen.w;
    matrix V = hydrogen.V;

    sort_eigen(energies, V);

    return {r, energies, V, dr};
}


double exact_energy(int principle_n){
    return -1.0 / (2.0 * principle_n * principle_n);
}

void output_convergence_dr(double rmax) {
    const std::vector<double> dr_values = {
        1.0,
        0.8,
        0.5,
        0.4,
        0.25,
        0.2,
        0.125,
        0.1
    };

    std::cout << std::setprecision(15);

    for (double current_dr : dr_values) {
        const HydrogenSolution solution =
            solve_hydrogen(rmax, current_dr);

        std::cout
            << current_dr << ' '
            << solution.energies[0] << ' '
            << exact_energy(1) << '\n';
    }
}

void output_convergence_rmax(double dr) {
    const std::vector<double> rmax_values = {
        2.0,
        3.0,
        4.0,
        5.0,
        6.0,
        8.0,
        10.0,
        12.0,
        15.0
    };

    std::cout << std::setprecision(15);

    for (double current_rmax : rmax_values) {
        const HydrogenSolution solution =
            solve_hydrogen(current_rmax, dr);

        std::cout
            << current_rmax << ' '
            << solution.energies[0] << ' '
            << exact_energy(1) << '\n';
    }
}



void output_wavefunctions(double rmax, double dr) {
    const HydrogenSolution solution =
        solve_hydrogen(rmax, dr);

    if (solution.energies.size() < 3) {
        throw std::runtime_error(
            "At least three eigenstates are required"
        );
    }

    const double normalization =
        1.0 / std::sqrt(dr);

    double signs[3] = {
        1.0,
        1.0,
        1.0
    };

    for (int k = 0; k < 3; ++k) {
        double overlap = 0.0;

        for (int i = 0; i < solution.r.size(); ++i) {
            const double radius = solution.r[i];

            double analytical = 0.0;

            if (k == 0) {
                analytical = exact_1s(radius);
            }
            else if (k == 1) {
                analytical = exact_2s(radius);
            }
            else {
                analytical = exact_3s(radius);
            }

            const double numerical =
                normalization * solution.V(i, k);

            overlap +=
                numerical * analytical * dr;
        }

        if (overlap < 0.0) {
            signs[k] = -1.0;
        }
    }

    std::cout << std::setprecision(15);

    for (int i = 0; i < solution.r.size(); ++i) {
        const double radius = solution.r[i];

        const double numerical_1s =
            signs[0] * normalization *
            solution.V(i, 0);

        const double numerical_2s =
            signs[1] * normalization *
            solution.V(i, 1);

        const double numerical_3s =
            signs[2] * normalization *
            solution.V(i, 2);

        std::cout
            << radius << ' '
            << numerical_1s << ' '
            << exact_1s(radius) << ' '
            << numerical_2s << ' '
            << exact_2s(radius) << ' '
            << numerical_3s << ' '
            << exact_3s(radius) << '\n';
    }
}

int main(int argc, char* argv[]){
    using namespace pp;

    if (argc == 1) {

        const int n = 5;
        const double scale = 1.0;
        const unsigned seed = 42;

        const double acc = 1e-8;
        const double eps = 1e-8;

        matrix A = rando_sym_M(n, scale, seed);
        matrix A_OG = A; 

        EVD evd(A);
        
        matrix V = evd.V;
        matrix D = diag(evd.w);

        matrix VT = V.T();
        matrix I = Id(n);

        matrix VTAV = VT * A_OG * V;
        matrix VDVT = V * D * VT;
        matrix VTV  = VT * V;
        matrix VVT  = V * VT;

        // Floating-point matrices must be compared approximately.
        bool test_VTAV = approx(VTAV, D, acc, eps);
        bool test_VDVT = approx(VDVT, A_OG, acc, eps);
        bool test_VTV  = approx(VTV,  I, acc, eps);
        bool test_VVT  = approx(VVT,  I, acc, eps);

        A_OG.print("Original symmetric matrix A:");
        D.print("Diagonal eigenvalue matrix D:");
        V.print("Eigenvector matrix V:");

        std::cout << std::boolalpha;

        std::cout << "Verification of the EVD:";
        std::cout << "V^T A V = D  : " << test_VTAV << '\n';
        std::cout << "V D V^T = A  : " << test_VDVT << '\n';
        std::cout << "V^T V = I    : " << test_VTV  << '\n';
        std::cout << "V V^T = I    : " << test_VVT  << '\n';

        const bool all_tests_passed =
            test_VTAV &&
            test_VDVT &&
            test_VTV &&
            test_VVT;

        std::cout << "All tests passed: "
                << all_tests_passed
                << '\n';

        // Return 0 if everything works, otherwise return 1.

        return  all_tests_passed ? 0 : 1;
    }
    
    std::string mode = "energies";
    double rmax = 10.0;
    double dr = 0.3; 

    for (int i = 1; i < argc; ++i) {
        const std::string argument = argv[i];

        if (argument == "--convergence-dr") {
            mode = "convergence-dr";
        }
        else if (argument == "--convergence-rmax") {
            mode = "convergence-rmax";
        }
        else if (argument == "--wavefunctions") {
            mode = "wavefunctions";
        }
        else if (
            argument == "-rmax" &&
            i + 1 < argc
        ) {
            rmax = std::stod(argv[++i]);
        }
        else if (
            argument == "-dr" &&
            i + 1 < argc
        ) {
            dr = std::stod(argv[++i]);
        }
        else {
            std::cerr
                << "Unknown argument: "
                << argument << '\n';

            return 1;
        }
    }
    if (mode == "convergence-dr") {
    output_convergence_dr(rmax);
    return 0;
    }

    if (mode == "convergence-rmax") {
        output_convergence_rmax(dr);
        return 0;
    }

    if (mode == "wavefunctions") {
        output_wavefunctions(rmax, dr);
        return 0;
    }
    std::cout << "Part B Hydrogen atom, s-wave radial schrodinger equation on a grid" << std::endl;
    const HydrogenSolution solution =
        solve_hydrogen(rmax, dr);

    const int num_states =
        std::min(3, solution.energies.size());

    std::cout << std::setprecision(12);

    std::cout << "Hydrogen atrom, s-wave states" << std::endl;
    std::cout << "rmax = " << rmax << std::endl;
    std::cout << "dr = " << dr << std::endl;
    std::cout << "N  =" << solution.r.size() << std::endl;


    std::cout << "n" << std::setw(22) << "Numerical Energy" << std::setw(22) << "Exact Energy" << std::setw(22) << "Error" << std::endl;

    for (int k = 0; k < num_states; ++k) {
        const int principle_n = k + 1;

        const double exact_energy = -1.0 / (2.0 * principle_n * principle_n);
        const double numerical_energy = std::abs(solution.energies[k] - exact_energy);

        std::cout << principle_n << std::setw(22) << solution.energies[k] << std::setw(22) << exact_energy << std::setw(22) << numerical_energy << std::endl;

    }
    return 0;
}