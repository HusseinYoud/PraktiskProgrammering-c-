#pragma once
#include "matrix.h"
#include<tuple>

//Entire jacobi.h file is taken from the homework 2, with some modifications to be closer to what i know about the subject.
namespace pp{
    struct EVD{
        vec w;
        matrix V;

        EVD(matrix A);

        static void timesJ(pp::matrix& A, int p, int q, double theta);
        static void Jtimes(pp::matrix& A, int p, int q, double theta);



    };





}