# quarticpp

**[ For a generic polynomial solver which includes also this class please see https://github.com/cridemichel/polypp ]**

The C++ class *quartic.hpp* provides a very accurate and very efficient quartic solver based on the paper 
ACM Transactions on Mathematical Software May 2020 Article No.: 20 https://doi.org/10.1145/3386241.
With this class you can solve both real and complex quartics in double or multiple precision (see below).
In addition to headers you will find some .cpp files with examples on how to use this class.
Multiprecision is implemented through boost multiprecision libraries (https://www.boost.org/doc/libs/1_73_0/libs/multiprecision/doc/html/index.html) and you need to have both boost and gmp (https://gmplib.org/) installed.
Boost and gmp are conveniently provided by *boost* and *gmp* homebrew packages (https://brew.sh/). 
Note that homebrew not only supports Mac OSX but also Linux and Windows (see https://docs.brew.sh/Homebrew-on-Linux).
To build the sources you need a c++ compiler which complies with *c++17*. It is strongly recommended to install g++ from homebrew. The package is called gcc and it will provide the g++ executable.

The class itself can be used without boost and gmp and it does not require gcc from homebrew. 
This means that the *quartic* example can be straightforwardly compiled by the command:
```shell
make quartic
```
without installing anything else.

Instead, by issuing the command:

```bash
make all
```

you obtain the following executables:

**quartic**: example of usage of quartic.hpp class for solving quartic equations

**quartic_mp**:  example of usage of quartic.hpp to solve quartic in multiple precision (using boost)

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
	(-1 perform both tests)

**accuracytest**: It performs the accuracy tests shown in Table 1 of the ACM paper.
The syntax is (where '>' is the shell prompt string):

```bash
> accuracytest <case>
```

where 
*case* is an integer between 1 and 24, which corresponds to the 24
cases shown in Table 1 of the ACM paper.

**Copyright © 2025 Cristiano De Michele**
Permission to make digital or hard copies of all or part of this work for personal or classroom use is granted without fee provided that copies are not made or distributed for profit or commercial advantage and that copies bear this notice and the full citation on the first page. Copyrights for components of this work owned by others than Cristiano De Michele must be honored. Abstracting with credit is permitted. To copy otherwise, or republish, to post on servers or to redistribute to lists, requires prior specific permission and/or a fee. Request permissions from cristiano.demichele@uniroma1.it.
	
**Citing**

If you use quarticpp in your research, please consider giving proper attribution by citing the following publication:
	
- A.G. Orellana and C. De Michele, Algorithm 1010: Boosting Efficiency in Solving Quartic Equations with No Compromise in Accuracy,
*ACM Trans. Math. Softw.* **46**, 20:1-20:28 (2020) https://doi.org/10.1145/3386241


  
