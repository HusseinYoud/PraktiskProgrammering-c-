#include"generalized.h"
#include"hydrogen.h"
#include"matrix.h"
#include"vector.h"
#include"minimizer.h"
#include<cmath>
#include<iomanip>
#include<numbers>
#include<iostream>
#include<functional>
#include<limits>
#include<algorithm>



pp::vec alpha_from_log_alpha(const pp::vec& log_alpha){
    pp::vec alpha(log_alpha.size());

    for (int i = 0; i < log_alpha.size(); ++i){
        alpha[i] = std::exp(log_alpha[i]);
    }
    return alpha;
}

double hydro_obj(const pp::vec& log_alpha){
    const pp::vec alpha = alpha_from_log_alpha(log_alpha);
    try {
        return pp::hydro_state(alpha).energy;
    }
    catch (const std::exception&){
        return 1e100;
    }
}

double overlap_eigenvalue_ratio(const pp::vec& alpha){
    pp::matrix H;
    pp::matrix N;
    pp::hydro_matrix(alpha, H,N);
    pp::EVD overlap_evd(N);
    double smallest = std::numeric_limits<double>::infinity();
    double largest = 0.0;
    for (int i = 0; i < overlap_evd.w.size(); ++i){
        smallest = std::min(smallest, overlap_evd.w[i]);
        largest =std::max(largest,overlap_evd.w[i]);
    }
    if (largest<=0.0){
        return 0.0;
    }
    return smallest / largest;
}
int main() {
    std::cout << std::boolalpha;
    std::cout << std::setprecision(15);
    {
        std::cout << "A. Implement a function to solve the generalized eigenvalue problem AV = BVE." << std::endl;
    const int n = 5;

    pp::matrix A = pp::rando_sym_M(n,1.0,1);
    pp::matrix R = pp::rando_M(n,n, 1.0,2 );
    pp::matrix B = R.T() * R + static_cast<double>(n) * pp::Id(n);

    pp::GEVD gevd(A,B);

    pp::matrix E = pp::diag(gevd.w);
    pp::matrix left = A * gevd.V;
    pp::matrix right = B * gevd.V * E;

    pp::matrix B_orthogonality = gevd.V.T() * B * gevd.V;
    
    std::cout << std::boolalpha;
    std::cout << "V^T B V = I: " << pp::approx(B_orthogonality, pp::Id(n),1e-8, 1e-8) << std::endl;
    std::cout << "A V = B V E: " << pp::approx(left,right, 1e-8,1e-8) << std::endl;

    gevd.w.print("Generalized eigen values : ");
    gevd.V.print("Generalized eigenvectors ");
    }
    {
        std::cout << "B. Finding the ground state of Hydrogen with the Schrödinger equation" << std::endl;

        for(int n = 1; n <= 5; ++n){
            pp::vec start_log_alpha(n);
            const double center = 0.5*static_cast<double>(n-1);

            for (int i = 0; i < n; ++i){
                const double initial_alpha = 0.3 * std::pow(4.0, static_cast<double>(i) - center);
                start_log_alpha[i] = std::log(initial_alpha);
            }
            const pp::simple_result result = pp::down_simple(hydro_obj, start_log_alpha, 0.5,1e-7,2000);
            const pp::vec alpha = alpha_from_log_alpha(result.x);

            try {
                const pp::hydrogenstate state = pp::hydro_state(alpha);
                const double exact_energy = -0.5;
                const double energy_error = std::abs(state.energy - exact_energy);
                const double overlap_ratio = overlap_eigenvalue_ratio(alpha);

                std::cout << std::left << std::setw(5) << n << std::setw(22) << std::setprecision(14) << state.energy << std::setw(22) << std::scientific << energy_error << std::setw(22) << overlap_ratio << std::defaultfloat << std::setw(14) << result.iterations << std::boolalpha << result.converged << std::endl;

                std::cout << "      The optimized alpha = ";

                for (int i = 0; i < alpha.size(); ++i){
                    std::cout << std::setprecision(10) << alpha[i] << "       ";
                }
                std::cout << std::endl;
            }
            catch(const std::exception& error) {
                std::cout << std::left << std:: setw(5) << n << "Failed " << error.what() << std::endl;
            }
        }
    }
    return 0;
}
