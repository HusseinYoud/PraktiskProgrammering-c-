CXX = g++

CXXFLAGS = -Wall -Wextra -std=c++23 -O2
LDLIBS = -lm

OBJECTS = main.o minimization.o lsfit.o matrix.o vector.o


all: main out.txt


main: $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(OBJECTS) -o $@ $(LDLIBS)


main.o: main.cc minimization.h matrix.h vector.h
	$(CXX) $(CXXFLAGS) -c main.cc -o main.o


minimization.o: minimization.cc minimization.h lsfit.h matrix.h vector.h
	$(CXX) $(CXXFLAGS) -c minimization.cc -o minimization.o


lsfit.o: lsfit.cc lsfit.h matrix.h vector.h
	$(CXX) $(CXXFLAGS) -c lsfit.cc -o lsfit.o


matrix.o: matrix.cc matrix.h vector.h
	$(CXX) $(CXXFLAGS) -c matrix.cc -o matrix.o


vector.o: vector.cc vector.h
	$(CXX) $(CXXFLAGS) -c vector.cc -o vector.o


out.txt: main
	./main > out.txt


clean:
	$(RM) main $(OBJECTS) out.txt


.PHONY: all clean