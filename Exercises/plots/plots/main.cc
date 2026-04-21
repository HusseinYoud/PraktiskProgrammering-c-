#include<iostream>
#include<cmath>
#include<iomanip>
#include "sfuns.h"
#include<string>   
#include<fstream> //Used to overwrite and read files. 

int main(int argc, char* argv[]){
    std::cout << std::scientific;
    std::string fun = "erf"; // default function
    std::string outfile = "data.txt"; // default output file

    double xmin=0, xmax=10, dx=0.125;
    for (int i=0; i < argc; i++){

        std::string arg = argv[i];
        std::cerr<<"i= " << i << " arg= " << arg << "\n";

        if(arg=="-xmin" && i+1<argc)xmin=std::stod(argv[i+1]);//hvor stod er string to double
        if(arg=="-xmax" && i+1<argc)xmax=std::stod(argv[i+1]);
        if(arg=="-dx" && i+1<argc)dx=std::stod(argv[i+1]);
        if(arg=="-fun" && i+1<argc)fun=argv[i+1];
        if(arg=="-out" && i+1<argc)outfile=argv[i+1];

    }
    std::cerr<<"fun= "<<fun<<"\n";
    std::cerr<<"xmin= "<<xmin<<"\n";
    std::cerr<<"xmax= "<<xmax<<"\n";
    std::cerr<<"dx= "<<dx<<"\n";
    std::cerr<<"outfile= "<<outfile<<"\n";
    


    std::ofstream out(outfile);
    out << std::scientific;

    for (double x=xmin; x<=xmax; x+=dx){
        double y = 0.0;
        bool ok = true;
        
        if(fun == "erf"){
            y = sfuns::erf(x);
        } else if(fun == "sgamma"){
            y = sfuns::sgamma(x);
        } else if(fun == "lngamma"){
            y = sfuns::lngamma(x);
        } else {
            std::cerr << "Øhhhh " << fun << "\n";
            return 1;
        }

        out << x << " " << y << "\n";
    }

    out << "\n\n\n";

    // Reference values for testing. 
    if (fun == "erf") {
        out << 0.0 << " " << 0.0                << "\n";
        out << 0.5 << " " << 0.5204998778130465 << "\n";
        out << 1.0 << " " << 0.8427007929497149 << "\n";
        out << 1.5 << " " << 0.9661051464753108 << "\n";
        out << 2.0 << " " << 0.9953222650189527 << "\n";
    }

    if (fun == "sgamma") {
        // Gamma(n) = (n-1)!
        out << 1.0 << " " << 1.0   << "\n";
        out << 2.0 << " " << 1.0   << "\n";
        out << 3.0 << " " << 2.0   << "\n";
        out << 4.0 << " " << 6.0   << "\n";
        out << 5.0 << " " << 24.0  << "\n";
        out << 6.0 << " " << 120.0 << "\n";
    }

    if (fun == "lngamma") {
        // ln(Gamma(n)) = ln((n-1)!)
        out << 1.0 << " " << std::log(1.0)   << "\n";
        out << 2.0 << " " << std::log(1.0)   << "\n";
        out << 3.0 << " " << std::log(2.0)   << "\n";
        out << 4.0 << " " << std::log(6.0)   << "\n";
        out << 5.0 << " " << std::log(24.0)  << "\n";
        out << 6.0 << " " << std::log(120.0) << "\n";
    }

    return 0;
}
