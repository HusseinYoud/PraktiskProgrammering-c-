CXX = g++
CXXFLAGS += -Wall -std=c++23 -O2 -Wextra 
LBLIBS += -lstdc++ -lm

all: output.txt
N=10

output.txt: main input.txt
	./main --input input.txt --output output.txt

input.txt:
	seq $(N) > input.txt

main.o : main.cc
	$(CXX) $(CXXFLAGS) -o main.o -c main.cc

main : main.o 
	$(CXX) -o main main.o $(LBLIBS)




.PHONY: clean
clean:
	$(RM) main *.o *.txt input.txt output.txt
