ifeq (,$(findstring intercept,$(CXX)))
CXX=g++
endif
#ifneq ($(CC),intercept-c)
ifeq (,$(findstring intercept,$(CC)))
CC=gcc
endif
BOOST_LIB=-L /usr/local/lib -lmpc -lmpfr
CXXFLAGS= -Wall -std=c++17 -g
CFLAGS= -Wall -O3 
HEADERS=./quartic.hpp ./pvector.hpp 
LDFLAGS=-lm -llapack -lblas $(BOOST_LIB) 
all: quartic quartic_mp

quartic: quartic.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -o quartic quartic.cpp  

quartic_mp: quartic_mp.cpp $(HEADERS) 
	$(CXX) $(CXXFLAGS) -I /usr/local/include/ $(LDFLAGS) -o quartic_mp quartic_mp.cpp  

clean:
	rm -f quartic quartic_mp *.o
