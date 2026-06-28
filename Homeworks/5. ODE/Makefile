CXX = g++
CXXFLAGS = -Wall -Wextra -std=c++23 -O2
LDLIBS = -lm

.PHONY: all clean

all: Out.txt oscillator.svg


Out.txt: main
	./$< > $@


oscillator.svg: Out.txt Makefile
	echo '\
		set terminal svg background rgb "white"; \
		set output "$@"; \
		set xlabel "x"; \
		set ylabel "Solution"; \
		set grid; \
		plot \
			"$<" using 1:2 with lines linewidth 2 \
				title "Numerical u(x)", \
			"$<" using 1:3 with lines dashtype 2 \
				title "Exact sin(x)", \
			"$<" using 1:4 with lines linewidth 2 \
				title "Numerical derivative", \
			"$<" using 1:5 with lines dashtype 2 \
				title "Exact cos(x)"; \
	' | gnuplot


main: main.o ode.o vector.o
	$(CXX) $(LDFLAGS) $^ -o $@ $(LDLIBS)


main.o: main.cc ode.h vector.h

ode.o: ode.cc ode.h vector.h

vector.o: vector.cc vector.h


clean:
	$(RM) main *.o Out.txt oscillator.svg