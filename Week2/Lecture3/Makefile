CXXFLAGS = -Wall -Wextra -std=c++17 -O
LINK.o = $(CXX) $(LDFLAGS)
LDLIBS += -LSTDC++ -lm

Out.txt: main
	./main > Out.txt

sfuns.o : sfuns.h sfuns.cc
	$(CXX) $(CXXFLAGS) -c sfuns.cc -o sfuns.o

main.o : sfuns.h main.cc
	$(CXX) $(CXXFLAGS) -c main.cc -o main.o

