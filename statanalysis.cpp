#define WP 60
#define WPQ 50 
// N.B. you can use either CPP, GMP or MPC backend by
// defining CPP_MP, GMP_MP or MPC_MP
#define MPC_MP
#ifdef CPP_MP
#include <boost/multiprecision/cpp_bin_float.hpp> 
#include <boost/multiprecision/cpp_complex.hpp>
using namespace boost;
using namespace boost::multiprecision;
using namespace boost::multiprecision::backends;
using mpreal = number<cpp_bin_float<WP>>;
using mpcmplx = cpp_complex<WP>;
using mprealQ = number<cpp_bin_float<WPQ>>;
using mpcmplxQ = cpp_complex<WPQ>;
#elif defined(GMP_MP)
#include <boost/multiprecision/gmp.hpp>
#include <boost/multiprecision/complex_adaptor.hpp>
using namespace boost;
using namespace boost::multiprecision;
using namespace boost::multiprecision::backends;
using mpreal=number<gmp_float<WP>>;
using mpcmplx=number<complex_adaptor<gmp_float<WP>>>;
using mprealQ=number<gmp_float<WPQ>>;
using mpcmplxQ=number<complex_adaptor<gmp_float<WPQ>>>;
#elif defined(MPC_MP)
#include <boost/multiprecision/mpc.hpp>
#include <boost/multiprecision/mpfr.hpp>
using namespace boost;
using namespace boost::multiprecision;
using namespace boost::multiprecision::backends;
//we set 100 digits working precision!
using mpreal=number<mpfr_float_backend<WP>>;
using mpcmplx=number<mpc_complex_backend<WP>>;
using mprealQ=number<mpfr_float_backend<WPQ>>;
using mpcmplxQ=number<mpc_complex_backend<WPQ>>;
#endif
#include<iostream>
#include<fstream>
#include"./quartic.hpp"
double ranf(void)
{
  return drand48();
}
#define PEPTS 500
int cmplxreal=0, restart, dojust=-1;
mpreal cumPEmp[PEPTS],PEmp[PEPTS]; 
mpreal cumPE[PEPTS], PE[PEPTS];
pvector<mpcmplx,4> csol;
pvector<mpcmplx,4> exsol;

int perm[24][4]={
      {0, 1, 2, 3}, {0, 1, 3, 2}, {0, 2, 1, 3}, {0, 2, 3, 1}, {0, 3, 1, 2}, {0, 3, 2, 1}, 
      {1, 0, 2, 3}, {1, 0, 3, 2}, {1, 2, 0, 3}, {1, 2, 3, 0}, {1, 3, 0, 2}, {1, 3, 2, 0}, 
      {2, 0, 1, 3}, {2, 0, 3, 1}, {2, 1, 0, 3}, {2, 1, 3, 0}, {2, 3, 0, 1}, {2, 3, 1, 0},
      {3, 0, 1, 2}, {3, 0, 2, 1}, {3, 1, 0, 2}, {3, 1, 2, 0}, {3, 2, 0, 1}, {3, 2, 1, 0}};

void save_PE(long long int numtrials, int numpts, mpreal dlogdE, mpreal logdEmin)
{
  fstream f;
  int k, kk;
  char fname[255];
  if (dojust == -1 || dojust == 0)
    {
      sprintf(fname,"P_of_eps_rel-dbl.dat");
      f.open(fname, ios::trunc|ios::out);
      for (k=0; k < numpts; k++)
        {
          if (PE[k]==0) 
            continue;
          f <<  k*dlogdE+logdEmin << " " <<  PE[k]/(mpreal(numtrials))/4. << "\n";
        }
      f.close();

    }
  if (dojust==-1 || dojust==1)
    {
      f.open("P_of_eps_rel-mp.dat", ios::trunc|ios::out);
      for (k=0; k < numpts; k++)
	{
	  if (PEmp[k]==0)
	    continue;
	  f <<  k*dlogdE+logdEmin << " " <<  PEmp[k]/(mpreal(numtrials))/4. << "\n";
	}
       f.close();
    }
  if (dojust == -1 || dojust==0)
    {  
      for (k=0; k < numpts; k++)
        {
          cumPE[k] = 0.0;
          for (kk=k; kk < numpts; kk++) 
            {
              cumPE[k] += PE[kk]/(mpreal(numtrials))/4.0;
            }
	}
    }
  if (dojust==-1 || dojust == 1)
    {
      for (k=0; k < numpts; k++)
	{
	  cumPEmp[k] = 0.0;
	  for (kk=k; kk < numpts; kk++) 
	    {
	      cumPEmp[k] += PEmp[kk]/(mpreal(numtrials))/4.0;
	    }
	}
    }
  if (dojust == -1 || dojust == 0)
    { 
      sprintf(fname,"F_of_eps_rel-dbl.dat");
      f.open(fname, ios::trunc|ios::out);
      for (k=0; k < numpts; k++)
        {
          if (cumPE[k]==0)
            continue;
          f << k*dlogdE+logdEmin << " " << cumPE[k] << "\n";
        }
      f.close();
    }
  if (dojust==-1 || dojust == 1)
    {
      f.open("F_of_eps_rel-mp.dat", ios::trunc|ios::out);
      if (cmplxreal!=5)
	{
	  for (k=0; k < numpts; k++)
	    {
	      if (cumPEmp[k]==0)
		continue;
              f <<  k*dlogdE+logdEmin << " " << cumPEmp[k] << "\n";
	    }
	}
      else
	f << "-16 1\n&\n";
      f.close();
    }
}

void sort_sol_opt(pvector<mpcmplx,4>& sol, pvector<mpcmplx,4>& exsol)
{
  int k1, k2, k1min=0;
  mpreal v, vmin;
  pvector<mpcmplx,4> solt;
  for (k1=0; k1 < 24; k1++)
    {
      v = 0;
      for (k2=0; k2 < 4; k2++)
	{
	  v += (exsol[k2]==0.0)?abs(sol[perm[k1][k2]]-exsol[k2]):abs((sol[perm[k1][k2]]-exsol[k2])/exsol[k2]);
	}
      if (k1==0 || v < vmin)
	{
	  k1min=k1;
	  vmin = v;
	}
    } 
  for (k2=0; k2 < 4; k2++)
    solt[k2] = sol[k2];

  for (k2=0; k2 < 4; k2++)
    sol[k2] = solt[perm[k1min][k2]];
}
void sort_sol_opt(pvector<complex<double>,4>& sol, pvector<complex<double>,4>& exsol)
{
  int k1, k2, k1min;
  double v, vmin;
  pvector<complex<double>,4> solt;
  for (k1=0; k1 < 24; k1++)
    {
      v = 0;
      for (k2=0; k2 < 4; k2++)
	{
	  v += (exsol[k2]==0.0)?abs(sol[perm[k1][k2]]-exsol[k2]):abs((sol[perm[k1][k2]]-exsol[k2])/exsol[k2]);
	}
      if (k1==0 || v < vmin)
	{
	  k1min=k1;
	  vmin = v;
	}
    } 
  for (k2=0; k2 < 4; k2++)
    solt[k2] = sol[k2];

  for (k2=0; k2 < 4; k2++)
    sol[k2] = solt[perm[k1min][k2]];
}

int main(int argc, char **argv)
{
  mpcmplx x1c, x2c, x3c, x4c; 
  mpreal logdE, dlogdE, logdEmax, logdEmin, sig, sig2;
  mpreal dE, x1;
  mpreal y1;
  pvector<double,5> c;
  pvector<mprealQ,5> cmpQ;
  pvector<mpreal,5> cmp;
  pvector<complex<double>,4> r;
  pvector<mpcmplx,4> csolmp;
  pvector<mpcmplxQ,4> rmpQ;
  quartic<double> Q;
  quartic<mprealQ,mpcmplxQ> Qmp;
  long long int dEzeromp=0, dEzero=0;
  long long int numtrials, its, numout, itsI;
  int numpts, ilogdE;
  int k, k2, ic=0, nsample;

  srand48(4242);
  
  for (k=0; k < PEPTS; k++)
    PE[k] = 0;
  for (k=0; k < PEPTS; k++)
    PEmp[k] = 0;

  sig = 1.0;
  sig2= 1.0;
  logdEmax=10.0;
  numpts = PEPTS; 
  logdEmin=-WPQ-1;
  //logdEmin=-18;
  dlogdE = (logdEmax -logdEmin)/numpts;
  //cout << "dlogdE=" << dlogdE << "\n";
  if (argc>=2)
    numtrials=atoll(argv[1]);
  else 
    numtrials=1000000000;

  restart = 0;
  itsI = 0;
  if (argc>=3)
    numout=atoll(argv[2]);
  else
    numout=100;
  nsample=cmplxreal;
  if (argc>=4)
    cmplxreal = nsample = atoi(argv[3]);
  if (cmplxreal < 0 || cmplxreal > 5)
    {
      cout << "cmplxreal must be between 0 and 5!\n";
      exit(-1);
    }
  if (cmplxreal==3)
    {
      sig = 1.0;
      sig2= 1E6;
      cmplxreal=1;
    }
  else if (cmplxreal==4)
    {
      sig = 1E6;
      sig2 = 1E6;
      cmplxreal = 2;
    } 
  if (argc  >= 5)
    dojust=atoi(argv[4]);
  if (dojust > 8)
    {
      cout << "which test should I have to perform?!?\n Last arg is too big!\n";
      exit(-1);
    }
  if (numtrials < 0)
    {
      cout << "number of trials must be a positive integer!\n";
      exit(-1);
    }
  else
    {
      cout << "numtrials=" << numtrials << " numout=" << numout << " cmplxreal=" << cmplxreal << " dojust=" << dojust << "\n"; 
    }

  x1c=x2c=x3c=x4c=0;
  for (its=itsI; its < numtrials; its++)
    {
      if (its > 0 && (its % (numtrials/numout) == 0))
	{
          if (nsample == 0)
	    printf("[SAMPLE A sig=%G %G]>>> its=%lld/%lld\n", double(sig), double(sig2), its, numtrials);
          else if (nsample==1)
	    printf("[SAMPLE B sig=%G %G]>>> its=%lld/%lld\n", double(sig), double(sig2), its, numtrials);
	  else if (nsample==2)
	    printf("[SAMPLE C sig=%G %G]>>> its=%lld/%lld\n", double(sig), double(sig2), its, numtrials);
          else if (nsample==3)
	    printf("[SAMPLE D sig=%G %G]>>> its=%lld/%lld\n", double(sig), double(sig2), its, numtrials);
	  else if (nsample==4)
            printf("[SAMPLE E sig=%G %G]>>> its=%lld/%lld\n", double(sig), double(sig2), its, numtrials);
          else if (nsample==5)
            printf("[SAMPLE F sig=%G %G]>>> its=%lld/%lld\n", double(sig), double(sig2), its, numtrials);
	  save_PE(its, numpts, dlogdE, logdEmin);
	  sync();  
	}
      /* generate 4 random roots */
      if (cmplxreal==2) /* 4 complex */
	{
	  x1 = sig2*(ranf()-0.5);
	  y1 = sig2*(ranf()-0.5);
	  x1c = mpcmplx(x1,y1);
	  x2c = mpcmplx(x1, -y1);
	  x1 = sig*(ranf()-0.5);
	  y1 = sig*(ranf()-0.5);
	  x3c = mpcmplx(x1, y1);
	  x4c = mpcmplx(x1, -y1);
	}
      else if (cmplxreal==1) /* two complex two real */
	{
	  x1 = sig2*(ranf()-0.5);
	  y1 = sig2*(ranf()-0.5);
	  x1c = mpcmplx(x1, y1);
	  x2c = mpcmplx(x1, -y1);
	  x1 = sig*(ranf()-0.5);
	  y1 = sig*(ranf()-0.5);
	  x3c = mpcmplx(x1);
	  x4c = mpcmplx(y1);
	}
      else if (cmplxreal==0)/* four real */
	{
#if 0
          x1c=x2c=x3c=1000.0*mpcmplx(sig*(ranf()-0.5));
          x4c = 1E-15*mpcmplx(sig*(ranf()-0.5));
#else
	  x1c = mpcmplx(sig*(ranf()-0.5));
	  x2c = mpcmplx(sig*(ranf()-0.5));
	  x3c = mpcmplx(sig*(ranf()-0.5));
	  x4c = mpcmplx(sig*(ranf()-0.5));
#endif
	}
      
      if (cmplxreal == 5)
	{
	  c[4]=1.0;
	  c[3]=ranf()-0.5;
	  c[2]=ranf()-0.5;
	  c[1]=ranf()-0.5;
	  c[0]=ranf()-0.5;
          for (int i=0; i < 5; i++)
            cmp[i] = mpreal(c[i]);
	}
      else
	{
	  cmp[4] = 1.0;
	  cmp[3] = mpcmplx(-(x1c+x2c+x3c+x4c)).real();
	  cmp[2] = mpcmplx(x1c*x2c + (x1c+x2c)*(x3c+x4c) + x3c*x4c).real(); 
	  cmp[1] = mpcmplx(-x1c*x2c*(x3c+x4c) - x3c*x4c*(x1c+x2c)).real();
	  cmp[0] = mpcmplx(x1c*x2c*x3c*x4c).real();

          for (int i=0; i < 5; i++)
            c[i] = double(cmp[i]); 
        
          exsol[0] = x1c;
	  exsol[1] = x2c;
	  exsol[2] = x3c;
	  exsol[3] = x4c;
	}

      if (dojust==-1 || dojust == 0)
	{
          Q.set_coeff(c);
          Q.find_roots(r);
          for (int i=0; i < 4; i++)
            csol[i] = mpcmplx(r[i]);
	}
      ic++;  
      if (dojust==-1 || dojust == 1 || cmplxreal==5)
	{
          for (int i=0; i < 5; i++)
            cmpQ[i] = mprealQ(cmp[i]);
          Qmp.set_coeff(cmpQ);
          Qmp.find_roots(rmpQ);
          for (int i=0; i < 4; i++)
            csolmp[i] = mpcmplx(rmpQ[i]);
	}
      if (cmplxreal==5)
	{
	  for (k2=0; k2 < 4; k2++) 
	    exsol[k2] = csolmp[k2];
	}
      if (dojust == -1 || dojust == 1)
	sort_sol_opt(csolmp, exsol);

      if (dojust==-1 || dojust==1) 
	{
	  for (k=0; k < 4; k++)
	    {
	      dE = (exsol[k]!=0)?abs((csolmp[k]- exsol[k])/exsol[k]):abs(csolmp[k]- exsol[k]); 
	      if (dE > 0.0)
		{
		  logdE=log10(dE)-logdEmin;
		  ilogdE=(int)(logdE/dlogdE);
		  if (ilogdE >= 0 && ilogdE < numpts)
		    {
		      (PEmp[ilogdE])++;
		    }
		}
              else 
                dEzeromp++;
	    }
	}
      if (dojust == -1 || dojust==0)
        {
          sort_sol_opt(csol, exsol);
          for (k=0; k < 4; k++)
            {	
              dE = (exsol[k]!=0)?abs((csol[k] - exsol[k])/exsol[k]):
                abs(csol[k] - exsol[k]); 
              if (dE > 0.0)
                {
                  logdE=log10(dE)-logdEmin;
                  ilogdE=(int)(logdE/dlogdE);
                  if (ilogdE >= 0 && ilogdE < numpts)
                    {
                      (PE[ilogdE])++;
                    }
                }
              else
                dEzero++;
            }
        }
    }
  if (dojust==-1 || dojust==0)
    cout << "Fraction of exact roots (double)=" << dEzero/(mpreal(numtrials)*4.0) << "\n";
  if (dojust==-1 || dojust==1)
    cout << "Fraction of exact roots (multiprecision)=" << dEzeromp/(mpreal(numtrials)*4.0) << "\n";

  save_PE(numtrials, numpts, dlogdE, logdEmin);
  cout << "Finished\n";
  sync();
  exit(0);
}
