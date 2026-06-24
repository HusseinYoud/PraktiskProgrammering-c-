CXX = g++
CXXFLAGS += -Wall -std=c++23 -O2 -Wextra 
LDLIBS += -lm -lstdc++ -lpthread

OBJECTS = main.o matrix.o vector.o lsfit.o
TARGET = main


all: Out.txt fit1.svg

$(TARGET) : $(OBJECTS)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDLIBS)
	
Out.txt: $(TARGET)
	./$(TARGET) > $@

main.o: main.cc lsfit.h matrix.h vector.h
matrix.o: matrix.cc matrix.h vector.h
vector.o: vector.cc vector.h
lsfit.o: lsfit.cc lsfit.h vector.h matrix.h

fit.data: Out.txt
	@test -f $@ ||(echo "Error: did not create $@" && exit 1)

fit1.svg: fit.data Makefile
	echo '\
	set terminal svg background rgb "white"; \
	set output "$@"; \
	set tics in; \
	set grid; \
	set xlabel "Time t (Days)";\
	set ylabel "Activity";\
	set title "Radioactive decay of ThX"; \
	set key outside right; \
	plot \
	"$<" index 0 using 1:2:3 with yerrorbars pointtype 7 title "Experimental data",\
	"$<" index 1 using 1:2 with lines linewidth 2 title "Best fit",\
	"$<" index 2 using 1:2 with lines dashtype 2 title "c0+dc0, c1+dc1",\
	"$<" index 3 using 1:2 with lines dashtype 3 title "c0+dc0, c1-dc1", \
	"$<" index 4 using 1:2 with lines dashtype 4 title "c0-dc0, c1+dc1", \
	"$<" index 5 using 1:2 with lines dashtype 5 title "c0-dc0, c1-dc1";\
	'| gnuplot


.PHONY: all clean
clean:
	rm -f main
	rm -f *.o
	rm -f Out.txt
	rm -f fit.data
	rm -f fit1.svg
