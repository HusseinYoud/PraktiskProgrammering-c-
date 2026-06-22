#include"generalized.h"
#include"matrix.h"
#include<iostream>


int main() {
    const int n = 5;

    pp::matrix A = pp::rando_sym_M(n,1.0,1);
    pp::matrix R = pp::rando_M(n,n, 1.0,2 );
    pp::matrix B = R.T() * R + static_cast<double>(n) * pp::Id(n);

    pp::GEVD gevd(A,B);

    pp::matrix E = pp::diag(gevd.w);
    pp::matrix left = A * gevd.V;
    pp::matrix right = B * gevd.V * E;

    std::cout << std::boolalpha;
    std:: cout << "V^T B V = I: " << pp::approx(B_orthogonality, pp::Id(n),1e-8, 1e-8) << std::endl;
    std::cout << "A V = B V E: " << pp::approx(left,right, 1e-8,1e-8) << std::endl;

    gevd.w.print("Generalized eigen values : ");
    gevd.V.print("Generalized eigenvectors ");

    return 0;
}
