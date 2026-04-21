#include<iostream>
#include<cmath>
#include<vector>
#include<fstream>

// int main (int argc, char* argv[]) {
// 	std::vector<double> numbers;
// 	std::string infile = "input.txt";
// 	std::string outfile = "output.txt";
// 	for(int i=0;i<argc;++i){
// 		std::string arg=argv[i];
// 		if(arg=="--input" && i+1<argc) infile =argv[i+1];
// 		if(arg=="--output" && i+1<argc) outfile =argv[i+1];
// 	}

// for (double x : numbers ){
//     std::cout << x << " " << std::sin(x) << " " << std::cos(x) << "\n";
// }

// if(infile.empty()) infile="input.txt";
// if(outfile.empty()) outfile="output.txt";


// std::ifstream myinput(infile);
// std::ofstream myoutput(outfile);


// std::cout << "# x sins(x) cos(x) \n";

// double x;
// 	while(std::cin >> x){
// 		std::cout << x << " " << std::sin(x) << " " << std::cos(x) << std::endl;

// 	}
// // double y;
// // 	std::ofstream output_file("output.txt");
// // 	while(std::cin >> y){
// // 		output_file << y << " " << std::sin(y) << " " << std::cos(y) << std::endl;
// // 	}
// // 	output_file.close();


// 	return EXIT_SUCCESS;
// }

int main (int argc, char *argv[]) {
	std::string infile="", outfile="";
	for(int i=0;i<argc;i++){
		std::string arg=argv[i];
		if(arg=="--input" && i+1 < argc) infile=argv[i+1];
		if(arg=="--output" && i+1 < argc) outfile=argv[i+1];
	}
	//Her definere jeg hvad input og output filen skal hedde, og hvis der ikke er givet nogen argumenter så sætter jeg dem til default værdierne "input.txt" og "output.txt"
	if (infile.empty()) infile="input.txt";
	if (outfile.empty()) outfile="output.txt";

	std::ifstream myinput(infile);
	std::ofstream myoutput(outfile);

double x;
if( myinput.is_open() && myoutput.is_open() ){
	while( myinput >> x ){
		myoutput << x <<" "<<std::sin(x)<<" "<<std::cos(x)<<std::endl;
		}
	}
else{
	std::cerr << "Error opening files: " << infile << outfile << std::endl;
	return EXIT_FAILURE;
    }

	//synes jeg lider lidt af at den laver input filen og ikke output filen, men det her virker?
myinput.close();
myoutput.close();
exit(EXIT_SUCCESS);
}