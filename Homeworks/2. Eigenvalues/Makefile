CXX = g++
CXXFLAGS += -Wall -std=c++23 -O2 -Wextra 
LDLIBS += -lm -lstdc++ -lpthread

OBJECTS = main.o matrix.o vector.o jacobi.o

CONVERGENCE_DATA = out.convergence.dr.txt out.convergence.rmax.txt
CONVERGENCE_PLOTS = convergence_dr.svg convergence_rmax.svg

.PHONY: all convergence-parallel timing-sequential timing-parallel clean

all: Out.txt \
	wavefunctions.svg \
	timing.txt \
	convergence-parallel

convergence-parallel: main
	+$(MAKE) --no-print-directory -j2 $(CONVERGENCE_PLOTS)

timing-sequential: main
	$(RM) $(CONVERGENCE_DATA)
	@echo "Sequential convergence calculations:"
	@/usr/bin/time -f "Real time: %e seconds" \
		$(MAKE) --no-print-directory -j1 $(CONVERGENCE_DATA)


timing-parallel: main
	$(RM) $(CONVERGENCE_DATA)
	@echo "Parallel convergence calculations:"
	@/usr/bin/time -f "Real time: %e seconds" \
		$(MAKE) --no-print-directory -j2 $(CONVERGENCE_DATA)

main: $(OBJECTS)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDLIBS)


Out.txt: main
	{ \
	echo "Exercise 1 "; \
	echo "-------------"; \
	echo; \
	./main; \
	echo; \
	echo "-------------"; \
	echo "Exercise 2 "; \
	echo "-------------"; \
	echo; \
	./main -rmax 10 -dr 0.3; \
	echo; \
	echo "-------------"; \
	} > $@

timing.txt: main
	@$(RM) $(CONVERGENCE_DATA)
	@echo "Timing of convergence calculations" > $@
	@echo "Sequential:" >> $@
	@/usr/bin/time -f "%e seconds" -o $@ -a \
		$(MAKE) --no-print-directory -j1 $(CONVERGENCE_DATA)
	@$(RM) $(CONVERGENCE_DATA)
	@echo "Parallel:" >> $@
	@/usr/bin/time -f "%e seconds" -o $@ -a \
		$(MAKE) --no-print-directory -j2 $(CONVERGENCE_DATA)


		
out.convergence.dr.txt : main
	./main --convergence-dr -rmax 10 > $@


convergence_dr.svg: out.convergence.dr.txt Makefile
	echo '\
	set terminal svg background rgb "white";\
	set output "$@";\
	set tics in; \
	set grid; \
	set xlabel "Grid spacing (dr)";\
	set ylabel "Ground-state energy"; \
	set title "Ground state convergence with respect to dr"; \
	plot \
		"$<" using 1:2 with linespoints title "Numerical energy",\
		"$<" using 1:3 with lines dashtype 2 title "Exact Energy E_1 = -0.5"; \
		' | gnuplot

out.convergence.rmax.txt : main
	./main --convergence-rmax -dr 0.2 > $@

convergence_rmax.svg: out.convergence.rmax.txt Makefile
	echo '\
	set terminal svg background rgb "white";\
	set output "$@";\
	set tics in; \
	set grid; \
	set xlabel "r_{max}"; \
	set ylabel "Ground-state energy"; \
	set title "Ground state convergence with respect to r_{max}"; \
	plot \
		"$<" using 1:2 with linespoints title "Numerical energy", \
		"$<" using 1:3 with lines dashtype 2 title "Exact Energy E_1 = -0.5"; \
		' | gnuplot

out.wavefunctions.txt: main
	./main --wavefunctions -rmax 30 -dr 0.3 > $@

wavefunctions.svg: out.wavefunctions.txt Makefile
	echo '\
	set terminal svg background rgb "white";\
	set output "$@";\
	set multiplot layout 3,1 \
		title "Wave functions, lowest s-wave eigenfunctions.";\
	set tics in; \
	set grid; \
	set xrange [0:30]; \
	set xlabel "r"; \
	set ylabel "f_1(r)";\
	plot \
		"$<" using 1:2 with lines title "Numerical 1s", \
		"$<" using 1:3 with lines dashtype 2 title "Analytical 1s"; \
	set ylabel "f_2(r)";\
	set xlabel "r"; \
	plot \
		"$<" using 1:4 with lines title "Numerical 2s", \
		"$<" using 1:5 with lines dashtype 2 title "Analytical 2s"; \
	set ylabel "f_3(r)";\
	set xlabel "r"; \
	plot \
		"$<" using 1:6 with lines title "Numerical 3s", \
		"$<" using 1:7 with lines dashtype 2 title "Analytical 3s"; \
		' | gnuplot


clean:
	rm -f main
	rm -f *.o
	rm -f Out.txt
	rm -f out.convergence.dr.txt
	rm -f out.convergence.rmax.txt
	rm -f out.wavefunctions.txt
	rm -f convergence_dr.svg
	rm -f convergence_rmax.svg
	rm -f wavefunctions.svg
	rm -f timing.txt