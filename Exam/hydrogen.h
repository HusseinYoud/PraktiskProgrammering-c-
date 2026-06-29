#pragma once
#include "matrix.h"
#include "vector.h"

namespace pp{
    //Hydrogenstate copy pasted from homework 2 but rewritten to fit exam better
    struct hydrogenstate {
        double energy;
        vec alpha;
        vec coefficients;
    };

    void hydro_matrix(const vec& alpha, matrix& H, matrix& N);
    
    hydrogenstate hydro_state(const vec& alpha);

    double hydro_wave(double r, const vec& alpha, const vec& coefficients);

}