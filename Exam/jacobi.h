#pragma once
#include "matrix.h"
#include<tuple>

namespace pp{
    struct EVD{
        vec w;
        matrix V;

        EVD(matrix A);

        static void timesJ(pp::matrix& A, int p, int q, double theta);
        static void Jtimes(pp::matrix& A, int p, int q, double theta);



    };





}