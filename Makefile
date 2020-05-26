ifeq (,$(findstring intercept,$(CXX)))
CXX=g++
endif
#ifneq ($(CC),intercept-c)
ifeq (,$(findstring intercept,$(CC)))
CC=gcc
endif
CXXFLAGS= -Wall -std=c++17 -O3
CFLAGS= -Wall -O3
HEADERS=./quartic.hpp ./pvector.hpp
LDFLAGS=-lm -llapack -lblas 

all: quartic

quartic: quartic.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -o quartic quartic.cpp  

clean:
	rm -f quartic  *.o
