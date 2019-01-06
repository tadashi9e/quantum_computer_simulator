#CXX=clang++
CXX=g++
CXXFLAGS=-Wall -O3 -std=c++11 -Iboost_system
LDFLAGS=-lboost_random

test: test_teleport test_add test_cluster

test_teleport: qc.o test_teleport.cc
	$(CXX) $(CXXFLAGS) qc.o test_teleport.cc $(LDFLAGS) -o test_teleport

test_add: qc.o test_add.cc
	$(CXX) $(CXXFLAGS) qc.o test_add.cc $(LDFLAGS) -o test_add

test_cluster: qc.o test_cluster.cc
	$(CXX) $(CXXFLAGS) qc.o test_cluster.cc $(LDFLAGS) -o test_cluster

clean:
	rm -f qc.o

qc.o: qc.cc
