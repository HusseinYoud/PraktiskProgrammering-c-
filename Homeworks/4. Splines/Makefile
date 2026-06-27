CXX = g++
CXXFLAGS += -Wall -std=c++23 -O2 -Wextra 
LDLIBS += -lm -lstdc++ -lpthread

all: Out.txt linear.svg integral.svg quadratic.svg quadratic_integral.svg quadratic_deriv.svg

main: main.o Spliner.o vector.o
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDLIBS)


main.o: main.cc Spliner.h vector.h

Spliner.o: Spliner.cc Spliner.h vector.h

vector.o: vector.cc vector.h


out.data Out.txt &: main
	./main > out.data 2> Out.txt


linear.svg: out.data Makefile
	echo '\
		set terminal svg background rgb "white"; \
		set output "$@"; \
		set tics in; \
		set grid; \
		set title "Linear spline interpolation"; \
		set xlabel "x"; \
		set ylabel "y"; \
		plot \
			"$<" index 0 using 1:2 with lines linewidth 2 \
				title "Linear spline", \
			"$<" index 1 using 1:2 with points pointtype 7 pointsize 0.8 \
				title "Data points", \
			cos(x) with lines dashtype 2 \
				title "cos(x)"; \
	' | gnuplot


integral.svg: out.data Makefile
	echo '\
		set terminal svg background rgb "white"; \
		set output "$@"; \
		set tics in; \
		set title "Integral of linear spline"; \
		set xlabel "z"; \
		set ylabel "Integral from 0 to z"; \
		set grid; \
		plot \
			"$<" index 0 using 1:3 with lines linewidth 2 \
				title "Spline integral", \
			sin(x) with lines dashtype 2 \
				title "Exact integral sin(x)"; \
	' | gnuplot

quadratic.svg: out.data Makefile
	echo '\
		set terminal svg background rgb "white"; \
		set output "$@"; \
		set tics in; \
		set title "Quadratic spline"; \
		set xlabel "x"; \
		set ylabel "y"; \
		set grid; \
		plot \
			"$<" index 2 using 1:2 with lines linewidth 2 \
				title "Quadratic spline", \
			cos(x) with lines dashtype 2 \
				title "cos(x)"; \
	' | gnuplot

quadratic_deriv.svg: out.data Makefile
	echo '\
		set terminal svg background rgb "white"; \
		set output "$@"; \
		set tics in; \
		set title "Derivative of quadratic spline"; \
		set xlabel "x"; \
		set ylabel "dy/dx"; \
		set grid; \
		plot \
			"$<" index 2 using 1:3 with lines linewidth 2 \
				title "Quadratic spline derivative", \
			-sin(x) with lines dashtype 2 \
				title "-sin(x)"; \
	' | gnuplot

quadratic_integral.svg: out.data Makefile
	echo '\
		set terminal svg background rgb "white"; \
		set output "$@"; \
		set tics in; \
		set title "Integral of quadratic spline"; \
		set xlabel "z"; \
		set ylabel "Integral from 0 to z"; \
		set grid; \
		plot \
			"$<" index 2 using 1:4 with lines linewidth 2 \
				title "Quadratic spline integral", \
			sin(x) with lines dashtype 2 \
				title "Exact integral sin(x)"; \
	' | gnuplot

clean:
	$(RM) main *.o out.data Out.txt linear.svg integral.svg quadratic.svg quadratic_integral.svg quadratic_deriv.svg


.PHONY: all clean