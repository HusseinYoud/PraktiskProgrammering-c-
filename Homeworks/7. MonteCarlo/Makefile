CXX = g++
CXXFLAGS = -Wall -Wextra -O2 -std=c++23
LDLIBS = -lm

OBJECTS = main.o mc.o vector.o

all: Out.txt circle-errors.svg quasi-errors.svg


main: $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(OBJECTS) -o $@ $(LDLIBS)


main.o: main.cc mc.h vector.h
	$(CXX) $(CXXFLAGS) -c main.cc -o main.o


mc.o: mc.cc mc.h vector.h
	$(CXX) $(CXXFLAGS) -c mc.cc -o mc.o


vector.o: vector.cc vector.h
	$(CXX) $(CXXFLAGS) -c vector.cc -o vector.o


Out.txt out.circle.txt out.quasi.txt &: main
	./main > Out.txt


circle-errors.svg: out.circle.txt Makefile
	echo '\
		set terminal svg size 900,600 background rgb "white"; \
		set output "$@"; \
		set logscale xy; \
		set grid; \
		set tics in; \
		set xlabel "Number of sampling points N"; \
		set ylabel "Absolute error"; \
		set title "Plain Monte Carlo integration of the unit circle"; \
		set key top right; \
		plot \
			"$<" using 1:2 with linespoints \
				title "Estimated error", \
			"$<" using 1:3 with linespoints \
				title "Actual error |Q_N - pi|", \
			"$<" using 1:4 with lines \
				linewidth 2 dashtype 2 \
				title "Expected 1/sqrt(N) scaling"; \
	' | gnuplot


quasi-errors.svg: out.quasi.txt Makefile
	echo '\
		set terminal svg size 1000,650 background rgb "white"; \
		set output "$@"; \
		set logscale xy; \
		set grid; \
		set tics in; \
		set xlabel "Total number of function evaluations N"; \
		set ylabel "Actual absolute error"; \
		set title "Pseudo-random and quasi-random convergence"; \
		set key top right; \
		plot \
			"$<" using 1:2 with linespoints \
				title "LCG actual error", \
			"$<" using 1:4 with linespoints \
				title "Standard C++ actual error", \
			"$<" using 1:6 with linespoints \
				title "Halton actual error", \
			"$<" using 1:8 with lines \
				linewidth 2 dashtype 2 \
				title "1/sqrt(N)", \
			"$<" using 1:9 with lines \
				linewidth 2 dashtype 3 \
				title "1/N"; \
	' | gnuplot


clean:
	$(RM) main
	$(RM) *.o
	$(RM) Out.txt
	$(RM) out.circle.txt
	$(RM) out.quasi.txt
	$(RM) circle-errors.svg
	$(RM) quasi-errors.svg


.PHONY: all clean