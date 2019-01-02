CXX=clang++
#CXX=g++
CXXFLAGS=-Wall -std=c++11 -Iboost_system
LDFLAGS=-lboost_random

test: test_teleport test_add

test_teleport: qc.o test_teleport.cc
	$(CXX) $(CXXFLAGS) qc.cc test_teleport.cc $(LDFLAGS) -o test_teleport

test_add: qc.o test_add.cc
	$(CXX) $(CXXFLAGS) qc.cc test_add.cc $(LDFLAGS) -o test_add

clean:
	rm -f qc.o

qc.o: qc.cc
