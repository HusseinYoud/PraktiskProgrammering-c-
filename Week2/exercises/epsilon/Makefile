CXX = g++
CXXFLAGS += -Wall -std=c++23 -O2 -Wextra 
LBLIBS += -lstdc++ -lm

Out.txt: main
	./main > Out.txt


sfuns.o : sfuns.h sfuns.cc
	$(CXX) $(CXXFLAGS) -o sfuns.o -c sfuns.cc
main.o : sfuns.h main.cc
	$(CXX) $(CXXFLAGS) -o main.o -c main.cc

main : main.o sfuns.o
	$(CXX) -o main main.o sfuns.o $(LDLIBS)

.PHONY: clean
clean:
	$(RM) main main.o Out.txt

