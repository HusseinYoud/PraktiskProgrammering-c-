#include "hydrogen.h"
#include "generalized.h"
#include <cmath>
#include <numbers>
#include <stdexcept>

namespace pp{

    void hydro_matrix(const vec& alpha, matrix& H, matrix& N){
        const int n = alpha.size();
        if (n < 1){
            throw std::invalid_argument("Hydrogen basis doesn't have a function");
        }
        H = matrix(n,n);
        N = matrix(n,n);
        const double sqrt_pi = std::sqrt(std::numbers::pi_v<double>);

        for (int i = 0; i < n; ++i){
            for (int j = i; j < n; ++j){
                const double ai = alpha[i];
                const double aj = alpha[j];
                const double sum = ai + aj;

                //We are trying to solve the integral Nij = integral (phi_i (r)) phi_j(r) dr
                const double Nij = sqrt_pi / (4.0 * std::pow(sum,1.5));

                const double Tij = 3.0 * sqrt_pi * ai * aj / (4.0 * std::pow(sum,2.5));

                const double Vij = -1.0 / (2.0 *sum);

                const double Hij = Tij + Vij;
                H(i,j) = Hij;
                H(j,i) = Hij;
                N(i,j) = Nij;
                N(j,i) = Nij;
            }
        }
    }
    hydrogenstate hydro_state(const vec& alpha){
        const int n = alpha.size();
        matrix H;
        matrix N;
        hydro_matrix(alpha, H,N);

        GEVD gevd(H,N);

        int ground_state_index = 0;
        
        for (int i = 1; i < n; ++i){
            if (gevd.w[i] < gevd.w[ground_state_index]) {
                ground_state_index = i;
            }
        } 
        const double energy = gevd.w[ground_state_index];

        vec coefficients(n);

        for (int i = 0; i < n; ++i){
            coefficients[i] = gevd.V(i,ground_state_index);
        }

        const vec Ncoefficients = N*coefficients;
        const double normal = std::sqrt(dot(coefficients, Ncoefficients));

        coefficients /= normal;

        //since the sign in eigenvectors are arbitrary i choose the sign such that u(1) is positive.
        if(hydro_wave(1.0,alpha,coefficients) < 0.0) {coefficients *= -1.0;} return {energy, alpha, coefficients};
    }
    double hydro_wave(double r, const vec& alpha, const vec& coefficients){
        double value = 0.0;
        for (int i = 0; i < alpha.size(); ++i){
            value += coefficients[i] * r * std::exp(-alpha[i] *r * r);
        }
        return value;
    }
}//namespace pp