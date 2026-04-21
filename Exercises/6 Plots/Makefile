CXX = g++
CXXFLAGS += -Wall -std=c++23 -O2 -Wextra 
LDLIBS += -lm

all: fig1.svg fig2.svg fig3.svg

Out.txt: main
	./main > Out.txt

main.o : main.cc sfuns.h
	$(CXX) $(CXXFLAGS) -o main.o -c main.cc

main : main.o sfuns.o
	$(CXX) -o main main.o sfuns.o $(LDLIBS)

sfuns.o: sfuns.cc sfuns.h
	$(CXX) $(CXXFLAGS) -o sfuns.o -c sfuns.cc

data_erf.txt: main
	./main -fun erf -xmin -2 -xmax 2 -dx 0.125 -out $@

data_gamma.txt: main
	./main -fun sgamma -xmin -5 -xmax 5 -dx 0.0625 -out $@

data_lngamma.txt: main
	./main -fun lngamma -xmin 0.05 -xmax 5 -dx 0.0625 -out $@


fig1.svg: data_erf.txt Makefile
	echo '\
		set terminal svg background rgb "white";\
		set output "$@";\
		set key top left;\
		set xlabel "x";\
		set ylabel "y";\
		set xrange [-2:2];\
		set yrange [-1.2:1.2];\
		set title "Plot of the error function";\
		set grid;\
		plot [] [-2:2]\
		"data_erf.txt" index 0 using 1:2 with points pt 7 ps 0.5 title "Approximate points",\
		"data_erf.txt" index 1 using 1:2 with lines pt 4 ps 1.5 title "erf(x)";\
		' | gnuplot

fig2.svg: data_gamma.txt Makefile
	echo '\
		set terminal svg background rgb "white";\
		set output "$@";\
		set key top left;\
		set xlabel "x";\
		set ylabel "Gamma(x)";\
		set xrange [-6:6];\
		set yrange [-5:5];\
		set title "Plot of the gamma function";\
		set grid;\
		plot \
		"data_gamma.txt" index 0 using 1:2 with points lw 2 title "Approximate points",\
		"data_gamma.txt" index 1 using 1:2 with lines pt 4 ps 1.5 title "Gamma(x)"\
		' | gnuplot

fig3.svg: data_lngamma.txt Makefile
	echo '\
		set terminal svg background rgb "white";\
		set output "$@";\
		set key top left;\
		set xlabel "x";\
		set ylabel "ln(Gamma(x))";\
		set xrange [-6:6];\
		set yrange [-5:5];\
		set title "Plot of the logarithm of the gamma function";\
		set grid;\
		plot \
		"data_lngamma.txt" index 0 using 1:2 with points lw 2 title "Approximate points",\
		"data_lngamma.txt" index 1 using 1:2 with lines pt 4 ps 1.5 title "Ln(gamma(x))"\
		' | gnuplot




.PHONY: clean
clean:
	$(RM) main *.o data_erf.txt data_gamma.txt data_lngamma.txt  *.svg

