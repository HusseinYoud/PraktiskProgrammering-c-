CXX = g++
CXXFLAGS += -Wall -std=c++23 -O2 -Wextra 
LDLIBS += -lm -lstdc++ -lpthread

all: main Out.txt fig1.svg

fig1.svg: times.data Makefile
	echo '\
		set terminal svg background rgb "white";\
		set output "$@";\
		set tics in; \
		set xlabel "Matrix size";\
		set ylabel "Time (s)";\
		a=1e-3; b=1e-6; \
		f(x) = a*x*x*x + b; \
		fit f(x) "times.data" using 1:2 via a,b; \
		plot \
			"$<" index 0 using 1:2 with linespoints title "Time taken", \
			f(x)  title "Fitted curve (O(n^3))"; \
		' | gnuplot

times.data: main Makefile
	>$@
	>err
	for N in $$(seq 100 15 1000); do \
		/usr/bin/time \
			--format "$$N\t%U" \
			--output $@ \
			--append \
			./$< -size $$N >/dev/null 2>>err; \
	done

Out.txt: main 
	./main > $@


main : main.o 
	$(CXX) -o $@ $^ $(LDLIBS) 

main.o : main.cc sfuns.h
	$(CXX) $(CXXFLAGS) -c $< -o $@



.PHONY: clean
clean: 
	rm out times.data fig1.svg err main main.o sfuns.o Out.txt fit.log
test: 
	@echo $(CXX) 
	@echo $(CXXFLAGS) 
	@echo $(RM)