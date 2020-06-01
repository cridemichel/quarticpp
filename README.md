# QUARTICPP
a very accurate and every efficient quartic solver based on the paper 
ACM Transactions on Mathematical Software May 2020 Article No.: 20 https://doi.org/10.1145/3386241.
In addition to headers you will find some .cpp files with examples on how to use this class.
Before compiling change the directories related to boost in the Makefile (BOOSTLIB and BOOSTHDR).

By issuing the command:

```bash
make all
```

you obtain the following executables:

**quartic**: example of usage of quartic.hpp class for solving quartic equations

**quartic_mp**:  example of usage of quartic.hpp to solve quartic in multiple precisione (using boost)

**quartic_cmplx**: solution of a complex quartic

**statanalysis**: it performs the statistical analyses shown in Figs. 2-4 of the ACM paper. 
The syntax is the following (where '>' is the shell prompt string):

```bash
> statanalysis <trials> <output> <sample> <solver>
```

where:

*trials*: it is the number of roots (samples A-E) or coefficients (sample F) to generate

*output*: every <output> trials save the the probability distribution function P(eps_rel) 
	in the file named P_of_eps_rel-XXX.dat and the cumulative distribution function F(eps_rel) 
	in the file named F_of_eps_rel-XXX.dat, where XXX can be dbl (for double precision) or mp (for multiprecision)   

*sample*: it is an integer between 0 and 5 which specifies the sample to generate 
	according to Table 4 with 0={sample A}, 1={sample B}, 2={sample C}, 3={sample D}, 4={sample E} and
	5={sample F}

*solver*: it is 0 or 1 which specifies the precision to use  
	for the analysis, where 0=dbl (double) and 1=mp (multiprecision),
	(-1 performance both tests)

**accuracytest**: It performs the accuracy tests shown in Table 1 of the ACM paper.
The syntax is (where '>' is the shell prompt string):

```bash
> accuracytest <case>
```

where 
*case* is an integer between 1 and 24, which corresponds to the 24
cases shown in Table 1 of the ACM paper.

