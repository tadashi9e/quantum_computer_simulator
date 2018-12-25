#CXX=clang++
CXX=g++
CXXFLAGS=-Wall -Iboost_system
LDFLAGS=-lboost_random

test_teleport: qc.o test_teleport.cc
	$(CXX) $(CXXFLAGS) qc.cc test_teleport.cc $(LDFLAGS) -o test_teleport

qc.o: qc.cc
