CXX = g++
CXXFLAGS += -Wall -std=c++23 -O2 -Wextra 
LDLIBS += -lm -lstdc++ -lpthread


SOURCES = main.cc generalized.cc jacobi.cc matrix.cc vector.cc
OBJECTS = $(SOURCES:.cc=.o)

all: Out.txt

Out.txt: main.exe
	./main.exe > Out.txt


main.exe: $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(OBJECTS) -o main.exe $(LDLIBS)

main.o: main.cc generalized.h jacobi.h matrix.h vector.h
generalized.o: generalized.cc generalized.h jacobi.h matrix.h vector.h
jacobi.o:jacobi.cc jacobi.h matrix.h vector.h
matrix.o: matrix.cc matrix.h vector.h
vector.o: vector.cc vector.h

clean:
	rm -f main.o generalized.o jacobi.o matrix.o vector.o
	rm -f main.exe Out.txt

.PHONY: all run clean