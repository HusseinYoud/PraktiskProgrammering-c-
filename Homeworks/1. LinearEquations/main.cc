#include "sfuns.h"
#include <iostream>
#include <string>


int main(int argc, char** argv){
    if (argc == 3 && std::string(argv[1]) == "-size"){
        int N = std::stoi(argv[2]);

        pp::matrix A = pp::rando_M(N, N, 1.0, static_cast<unsigned>(N));

        pp::qr decomp(A);

        return 0;
    }
    {
    std::cout << "Generating random matrix A";
    int n = 7, m = 6;
    pp::matrix A = pp::rando_M(n, m, 1.0, 11);
    A.print("Random Matrix A:");

        std::cout << "Performing QR decomposition on A" << std::endl;
    pp::qr decomp(A);
    decomp.Q.print("Q:");
    decomp.R.print("R:");

        std::cout << "Checking if Q is orthogonal (Q^T * Q should be identity)" << std::endl;
    pp::matrix QtQ = decomp.Q.T() * decomp.Q;
    QtQ.print("Q^T * Q:");

        std::cout << "Checking if A = Q * R" << std::endl;
    pp::matrix QR = decomp.Q * decomp.R;
    QR.print("Q * R:");
    
    std::cout << "Checking if QR = A" << pp::approx(QR,A,1e-8,1e-8) << std::endl;

        std::cout << "Checking if R is upper triangular" << std::endl;
    pp::uptrig(decomp.R) ? std::cout << "R is upper triangular." << std::endl : std::cout << "R is not upper triangular." << std::endl;    
}

    {
        int rows = 6;

        pp::matrix A = pp::rando_M(rows, rows, 1.0, 15);
        pp::vec b = pp::rando_vec(rows, 1.0, 16);

        pp::qr decomp(A);
        pp::vec x = decomp.solve(b);
        pp::vec Ax = A * x;

        b.print("b:");
        x.print("Solution x:");
        Ax.print("A * x:");

        std::cout << "Checking if A * x is approximately equal to b" << pp::approx(Ax,b) << std::endl;
        std::cout << "The determinant of A from R: " << decomp.det() << std::endl;
    }

    {
        int rows = 5;
        pp::matrix A = pp::rando_M(rows, rows, 1.0, 20);

        A.print("Original Matrix A:");
        pp::qr decomp(A);
        pp::matrix A_inv = decomp.inverse();
        pp::matrix AB = A*A_inv;

        pp::matrix I(rows,rows);
        I.setid();
        I.print("Identity Matrix I:");

        A_inv.print("Inverse of A:");
        AB.print("A * A^-1:");
        std::cout << "Checking if A * A^-1 is approximately equal to the identity matrix" << pp::approx(AB, I) << std::endl;
        std::cout << "Checking if A*A^{-1} is identity matrix" << pp::approx(AB, I, 1e-8, 1e-8) << std::endl;


    }
    {
        pp::matrix A = pp::rando_M(4, 4, 1.0, 25);
        A.print("Original Matrix A:");
        pp::qr decomp(A);
    }
    return 0;
}