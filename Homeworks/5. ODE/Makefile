CXX = g++
CXXFLAGS = -Wall -Wextra -std=c++23 -O2
LDLIBS = -lm

DATA = \
	oscillator.data \
	orbit-circular.data \
	orbit-elliptic.data \
	orbit-relativistic.data \
	threebody.data

PLOTS = \
	oscillator.svg \
	orbit-circular.svg \
	orbit-elliptic.svg \
	orbit-relativistic.svg \
	threebody.svg

.PHONY: all clean

all: Out.txt $(PLOTS)


Out.txt $(DATA) &: main
	./main > Out.txt


main: main.o ode.o vector.o
	$(CXX) $(LDFLAGS) $^ -o $@ $(LDLIBS)


main.o: main.cc ode.h vector.h

ode.o: ode.cc ode.h vector.h

vector.o: vector.cc vector.h


oscillator.svg: oscillator.data Makefile
	echo '\
		set terminal svg background rgb "white"; \
		set output "$@"; \
		set xlabel "x"; \
		set ylabel "u(x)"; \
		set grid; \
		plot \
			"$<" using 1:2 with lines \
				title "Numerical u(x)", \
			"$<" using 1:3 with lines dashtype 2 \
				title "Exact sin(x)"; \
	' | gnuplot


orbit-%.svg: orbit-%.data Makefile
	echo '\
		set terminal svg background rgb "white"; \
		set output "$@"; \
		set size ratio -1; \
		set xlabel "x"; \
		set ylabel "y"; \
		set grid; \
		set title "$* orbit"; \
		plot \
			"$<" using \
				(cos($$1)/$$2):(sin($$1)/$$2) \
				with lines notitle; \
	' | gnuplot


threebody.svg: threebody.data Makefile
	echo '\
		set terminal svg background rgb "white"; \
		set output "$@"; \
		set size ratio -1; \
		set xlabel "x"; \
		set ylabel "y"; \
		set grid; \
		set title "Three-body figure-eight orbit"; \
		plot \
			"$<" using 2:3 with lines title "Body 1", \
			"$<" using 4:5 with lines title "Body 2", \
			"$<" using 6:7 with lines title "Body 3"; \
	' | gnuplot


clean:
	$(RM) main *.o Out.txt $(DATA) $(PLOTS)