#pragma once

#include"matrix.h"
#include"vector.h"
#include"jacobi.h"

namespace pp{
    struct GEVD{
        vec w;
        matrix V;

        GEVD(const matrix& A, const matrix& B);


    };



} //namespace pp