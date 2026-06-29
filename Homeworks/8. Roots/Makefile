CXX = g++
CXXFLAGS = -Wall -Wextra -std=c++23 -O2
LDLIBS = -lm

PROGRAM = main

OBJECTS = \
	main.o \
	hydrogen.o \
	roots.o \
	ode.o \
	matrix.o \
	vector.o

DATAFILES = \
	hydrogen-wavefunction.data \
	convergence-rmax.data \
	convergence-rmin.data \
	convergence-acc.data \
	convergence-eps.data

FIGURES = \
	hydrogen-wavefunction.svg \
	convergence-rmax.svg \
	convergence-rmin.svg \
	convergence-acc.svg \
	convergence-eps.svg


all: data $(FIGURES)


data: $(PROGRAM)
	./$(PROGRAM) > Out.txt


$(PROGRAM): $(OBJECTS)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDLIBS)


main.o: main.cc hydrogen.h roots.h vector.h

hydrogen.o: hydrogen.cc hydrogen.h roots.h ode.h vector.h

roots.o: roots.cc roots.h matrix.h vector.h

ode.o: ode.cc ode.h vector.h

matrix.o: matrix.cc matrix.h vector.h

vector.o: vector.cc vector.h


hydrogen-wavefunction.svg: data Makefile
	echo '\
		set terminal svg background rgb "white"; \
		set output "$@"; \
		set grid; \
		set key top right; \
		set xlabel "r / a_0"; \
		set ylabel "Reduced radial wave-function f(r)"; \
		set title "Hydrogen ground-state wave-function"; \
		plot \
			"hydrogen-wavefunction.data" using 1:2 \
				with lines linewidth 2 \
				title "Numerical shooting solution", \
			"hydrogen-wavefunction.data" using 1:3 \
				with lines dashtype 2 linewidth 2 \
				title "Exact: r exp(-r)"; \
	' | gnuplot


convergence-rmax.svg: data Makefile
	echo '\
		set terminal svg background rgb "white"; \
		set output "$@"; \
		set grid; \
		set logscale y; \
		set xlabel "r_max / a_0"; \
		set ylabel "|E - E_0|"; \
		set title "Convergence with respect to r_max"; \
		plot \
			"convergence-rmax.data" using 1:3 \
				with linespoints linewidth 2 \
				title "Energy error"; \
	' | gnuplot


convergence-rmin.svg: data Makefile
	echo '\
		set terminal svg background rgb "white"; \
		set output "$@"; \
		set grid; \
		set logscale xy; \
		set format x "10^{%L}"; \
		set format y "10^{%L}"; \
		set xlabel "r_min / a_0"; \
		set ylabel "|E - E_0|"; \
		set title "Convergence with respect to r_min"; \
		plot \
			"convergence-rmin.data" using 1:3 \
				with linespoints linewidth 2 \
				title "Energy error"; \
	' | gnuplot


convergence-acc.svg: data Makefile
	echo '\
		set terminal svg background rgb "white"; \
		set output "$@"; \
		set grid; \
		set logscale xy; \
		set format x "10^{%L}"; \
		set format y "10^{%L}"; \
		set xlabel "ODE absolute accuracy acc"; \
		set ylabel "|E - E_0|"; \
		set title "Convergence with respect to ODE acc"; \
		plot \
			"convergence-acc.data" using 1:3 \
				with linespoints linewidth 2 \
				title "Energy error"; \
	' | gnuplot


convergence-eps.svg: data Makefile
	echo '\
		set terminal svg background rgb "white"; \
		set output "$@"; \
		set grid; \
		set logscale xy; \
		set format x "10^{%L}"; \
		set format y "10^{%L}"; \
		set xlabel "ODE relative accuracy eps"; \
		set ylabel "|E - E_0|"; \
		set title "Convergence with respect to ODE eps"; \
		plot \
			"convergence-eps.data" using 1:3 \
				with linespoints linewidth 2 \
				title "Energy error"; \
	' | gnuplot


clean:
	$(RM) \
		$(PROGRAM) \
		*.o \
		Out.txt \
		$(DATAFILES) \
		$(FIGURES)


.PHONY: all data clean