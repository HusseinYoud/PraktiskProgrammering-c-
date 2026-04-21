CXX = g++
CXXFLAGS += -Wall -std=c++23 -O2 -Wextra 
LBLIBS += -lstdc++ -lm

Out.txt: main
	./main > Out.txt


vec.o: vec.h vec.cc
	$(CXX) $(CXXFLAGS) -o vec.o -c vec.cc

main.o : vec.h main.cc
	$(CXX) $(CXXFLAGS) -o main.o -c main.cc

main : main.o vec.o
	$(CXX) -o main main.o vec.o $(LDLIBS)

.PHONY: clean
clean:
	$(RM) main *.o Out.txt

