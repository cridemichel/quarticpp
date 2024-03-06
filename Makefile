BREWINST=.brew_packages_installed_
ifneq ($(MAKECMDGOALS),quartic)
ifeq ("$(wildcard $(BREWINST))","")
ifeq ($(shell command -v brew 2> /dev/null),)
  $(error Please install homebrew first!)
endif
# check if homebrew is installed 
HBPACK=$(shell brew ls --versions gmp gcc boost) 
ifeq ($(shell echo $(HBPACK)),)
  $(warning Please install gmp, boost and gcc homebrew packages through the command:)
  $(warning > brew install gmp boost gcc)
  $(error Aborting...)
endif
ifeq ($(shell echo $(HBPACK)|grep gmp),)
  $(warning Please install gmp homebrew packages through the command:)
  $(warning > brew install gmp)
  $(error Aborting...)
  GMPPAK=0
else	
  GMPPAK=1
endif
ifeq ($(shell echo $(HBPACK)|grep gcc),)
  $(warning Please install gcc homebrew packages through the command:)
  $(warning > brew install gcc)
  $(error Aborting...)
  GCCPAK=0
else	
  $(shell CC=$[$CC+1])
  GCCPAK=1
endif
ifeq ($(shell echo $(HBPACK)|grep boost),)
  $(warning Please install boost homebrew packages through the command:)
  $(warning > brew install boost)
  $(error Aborting...)
  BOOSTPAK=0
else
  BOOSTPAK=1
endif
ifeq ($(GMPPAK),1)
  ifeq ($(GCCPAK),1)
    ifeq ($(BOOSTPAK),1) 
      $(shell touch $(BREWINST))
    endif
  endif
endif
endif
endif
ifneq ($(shell command -v brew 2> /dev/null),)
ifeq ($(HBDIR),)
  HBDIR=$(shell brew --prefix)
endif
endif
ifeq (,$(findstring intercept,$(CXX)))
  #check if GNU g++ is installed
  ifneq ("$(wildcard $(HBDIR))","")
    CXXHB=$(shell ls $(HBDIR)/bin/g++-[0-9]* | sort -t - -k 2 -n | tail -1)
    ifneq ("$(wildcard $(CXXHB))","")
       	CXX=$(shell basename $(CXXHB))
        $(info Found GNU $(CXX) from homebrew)
    else	
	CXX=g++
    endif
  else
    $(info Homebrew not installed, I am going to use g++)
    CXX=g++
  endif
endif
ifneq ($(HBDIR),)
HBLIBS=-L $(HBDIR)/lib -lmpc -lmpfr -lgmp -lgmpxx
HBHDRS=-I $(HBDIR)/include
else
HBLIBS=
HBHDRS=
endif
LIBS=$(HBLIBS) 
CXXFLAGS= -Wall -std=c++17 -O3 
CXXFLAGSMP=$(CXXFLAGS) $(HBHDRS)
HEADERS=./quartic.hpp ./pvector.hpp 
LDFLAGS=-lm $(LIBS) 
all: quartic quartic_mp quartic_cmplx accuracytest statanalysis

quartic: quartic.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -o quartic quartic.cpp  

quartic_mp: quartic_mp.cpp $(HEADERS) 
	$(CXX) $(CXXFLAGSMP) $(LDFLAGS) -o quartic_mp quartic_mp.cpp  

quartic_cmplx: quartic_cmplx.cpp $(HEADERS) 
	$(CXX) $(CXXFLAGSMP) $(LDFLAGS) -o quartic_cmplx quartic_cmplx.cpp  

accuracytest: accuracytest.cpp $(HEADERS)
	$(CXX) $(CXXFLAGSMP) $(LDFLAGS) -o accuracytest accuracytest.cpp  

statanalysis: statanalysis.cpp $(HEADERS)
	$(CXX) $(CXXFLAGSMP) $(LDFLAGS) -o statanalysis statanalysis.cpp  

clean:
	rm -f quartic quartic_mp quartic_cmplx accuracytest statanalysis *.o $(BREWINST)
