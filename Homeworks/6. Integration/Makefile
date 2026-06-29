CXX = g++
CXXFLAGS = -Wall -Wextra -std=c++23 -O2
LDLIBS = -lm

all: Out.txt erf.svg erf-convergence.svg error-estimate.svg

main: main.o integration.o
	$(CXX) $(LDFLAGS) $^ $(LDLIBS) -o $@

main.o: main.cc integration.h

integration.o: integration.cc integration.h

Out.txt erf.data erf-convergence.data error-estimate.data &: main
	./main > Out.txt

erf.svg: erf.data Makefile
	echo '\
		set terminal svg background rgb "white"; \
		set output "$@"; \
		set xlabel "z"; \
		set ylabel "erf(z)"; \
		set grid; \
		set key top left; \
		plot \
			"$<" using 1:2 with points pointtype 7 \
				title "integral", \
			"$<" using 1:3 with lines \
				title "std::erf"; \
	' | gnuplot

erf-convergence.svg: erf-convergence.data Makefile
	echo '\
		set terminal svg background rgb "white"; \
		set output "$@"; \
		set logscale xy; \
		set xlabel "acc"; \
		set ylabel "|erf(1) - exact erf(1)|"; \
		set grid; \
		set key top left; \
		plot \
			"$<" using 1:2 with linespoints \
				title "Observed error"; \
	' | gnuplot

error-estimate.svg: error-estimate.data Makefile
	echo '\
		set terminal svg background rgb "white"; \
		set output "$@"; \
		set logscale xy; \
		set xlabel "Requested absolute accuracy"; \
		set ylabel "Error"; \
		set grid; \
		set key top left; \
		plot \
			"$<" using 1:2 with linespoints \
				title "Actual error", \
			"$<" using 1:3 with linespoints \
				title "Estimated error"; \
	' | gnuplot

clean:
	rm -f main main.o integration.o \
		Out.txt \
		erf.data \
		erf-convergence.data \
		error-estimate.data \
		erf.svg \
		erf-convergence.svg \
		error-estimate.svg

.PHONY: all clean