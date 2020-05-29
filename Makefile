ifeq (,$(findstring intercept,$(CXX)))
CXX=g++
endif
#ifneq ($(CC),intercept-c)
ifeq (,$(findstring intercept,$(CC)))
CC=gcc
endif
BOOST_LIB=-L /usr/local/lib -lmpc -lmpfr
CXXFLAGS= -Wall -std=c++17 -O3 -I /usr/local/include/ $(LDFLAGS) 
CFLAGS= -Wall -O3 
HEADERS=./quartic.hpp ./pvector.hpp 
LDFLAGS=-lm -llapack -lblas $(BOOST_LIB) 
all: quartic quartic_mp accuracytest statanalysis

quartic: quartic.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -o quartic quartic.cpp  

quartic_mp: quartic_mp.cpp $(HEADERS) 
	$(CXX) $(CXXFLAGS) -I /usr/local/include/ $(LDFLAGS) -o quartic_mp quartic_mp.cpp  

accuracytest: accuracytest.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -I /usr/local/include/ $(LDFLAGS) -o accuracytest accuracytest.cpp  

statanalysis: statanalysis.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -I /usr/local/include/ $(LDFLAGS) -o statanalysis statanalysis.cpp  

clean:
	rm -f quartic quartic_mp *.o
