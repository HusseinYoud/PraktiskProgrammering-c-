CXX = g++
CXXFLAGS += -Wall -std=c++23 -O2 -Wextra 
LBLIBS += -lstdc++ -lm

all: data figs

Out.txt: main
	./main > Out.txt



main.o : main.cc sfuns.h
	$(CXX) $(CXXFLAGS) -o main.o -c main.cc

main : main.o sfuns.o
	$(CXX) -o main main.o sfuns.o $(LDLIBS)

sfuns.o: sfuns.cc sfuns.h
	$(CXX) $(CXXFLAGS) -o sfuns.o -c sfuns.cc


data: data_erf.txt data_gamma.txt

data_erf.txt: main Makefile
	./main -xmin -2 -xmax 2 -erf > "$@"

data_gamma.txt: main Makefile
	./main -xmin -5 -xmax 5 -gamma -dx 0.0625 > "$@"



.PHONY: clean
clean:
	$(RM) main *.o data.txt *.svg

