#pragma once
#include "vector.h"

namespace pp{
    int binsearch(const vec& x, double z);

    double linterp(const vec& x, const vec& y, double z);
    double linterp_intg(const vec& x, const vec& y, double z);

    struct qspline{vec x; vec y; vec b; vec c; 
        qspline(const vec& x, const vec& y);

        int binsearch(double z) const;
        double eval(double z) const;
        double deriv(double z) const;
        double intg(double z) const;
        
        
    };



}