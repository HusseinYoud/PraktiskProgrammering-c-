#pragma once
#include "vector.h"
#include<functional>
//Debugging helped with by ChatGPT 
namespace pp{
    struct simple_result{
        vec x;
        double fx;
        int iterations;
        bool converged;
    };
    simple_result down_simple(const std::function<double(const vec&)>& f,const vec& start, double step = 0.5, double tol = 1e-7, int max_iterations = 1000);
}//namespace pp