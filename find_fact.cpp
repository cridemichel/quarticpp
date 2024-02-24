#include <ctime>
#include <complex>
#define WP 128
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
#elif defined(GMP_MP)
#include <boost/multiprecision/gmp.hpp>
#include <boost/multiprecision/complex_adaptor.hpp>
using namespace boost;
using namespace boost::multiprecision;
using namespace boost::multiprecision::backends;
using mpreal=number<gmp_float<WP>>;
using mpcmplx=number<complex_adaptor<gmp_float<WP>>>;
#elif defined(MPC_MP)
#include <boost/multiprecision/mpc.hpp>
#include <boost/multiprecision/mpfr.hpp>
using namespace boost;
using namespace boost::multiprecision;
using namespace boost::multiprecision::backends;
//we set 100 digits working precision!
using mpreal= number<mpfr_float_backend<WP>>;
using mpcmplx= number<mpc_complex_backend<WP>>;
#endif
#include"./quartic.hpp"
int perm[24][4]={
      {0, 1, 2, 3}, {0, 1, 3, 2}, {0, 2, 1, 3}, {0, 2, 3, 1}, {0, 3, 1, 2}, {0, 3, 2, 1}, 
      {1, 0, 2, 3}, {1, 0, 3, 2}, {1, 2, 0, 3}, {1, 2, 3, 0}, {1, 3, 0, 2}, {1, 3, 2, 0}, 
      {2, 0, 1, 3}, {2, 0, 3, 1}, {2, 1, 0, 3}, {2, 1, 3, 0}, {2, 3, 0, 1}, {2, 3, 1, 0},
      {3, 0, 1, 2}, {3, 0, 2, 1}, {3, 1, 0, 2}, {3, 1, 2, 0}, {3, 2, 0, 1}, {3, 2, 1, 0}};

#if 1
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
	  v += (abs(exsol[k2])==0.0)?abs(sol[perm[k1][k2]]-exsol[k2]):abs((sol[perm[k1][k2]]-exsol[k2])/exsol[k2]);
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
#endif
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
mpreal print_accuracy_at(pvector<mpcmplx,4> csol, pvector<mpcmplx,4> exsol)
{
  /* we follow FLocke here */
  int k1;
  mpreal relerr, relerrmax;
  for (k1=0; k1 < 4; k1++)
    {
      relerr=abs((csol[k1] - exsol[k1])/exsol[k1]); 
      if (k1==0 || relerr > relerrmax)
        {
          relerrmax=abs((csol[k1] - exsol[k1])/exsol[k1]); 
        }
    }
  cout << "relative accuracy="<< relerrmax << "\n";
  return relerrmax;
}
#if 1
double print_accuracy_at(pvector<complex<double>,4> csol, pvector<complex<double>,4> exsol)
{
  /* we follow FLocke here */
  int k1;
  double relerr, relerrmax;
  for (k1=0; k1 < 4; k1++)
    {
      relerr=abs((csol[k1] - exsol[k1])/exsol[k1]); 
      if (k1==0 || relerr > relerrmax)
        {
          relerrmax=abs((csol[k1] - exsol[k1])/exsol[k1]); 
        }
    }
  cout << "relative accuracy="<< relerrmax << "\n";
  return relerrmax;
}
#endif
void print_roots(const char *str, mpcmplx x1c, mpcmplx x2c,
                 mpcmplx x3c, mpcmplx x4c)
{
  cout << str;
  cout << setprecision(50) << "root #1=" << x1c.real() << "+I*(" << x1c.imag() <<  ")\n";
  cout << setprecision(50) << "root #2=" << x2c.real() << "+I*(" << x2c.imag() <<  ")\n";
  cout << setprecision(50) << "root #3=" << x3c.real() << "+I*(" << x3c.imag() <<  ")\n";
  cout << setprecision(50) << "root #4=" << x4c.real() << "+I*(" << x4c.imag() <<  ")\n";
}

int main(int argc, char **argv)
{
  quartic<double> Q;
  quartic<mpreal,mpcmplx> Qmp;
  pvector<double,5> cpv;
  pvector<mpreal,5> cpvmp;
  pvector<complex<double>,4> r;
  pvector<mpcmplx,4> rmp;
  mpcmplx x1c, x2c, x3c, x4c; 
  mpreal c[5], epsfact=1.0, tmpfact, A, B, err, C, D;
  pvector<complex<double>,4> csol, csolREF;
  pvector<mpcmplx,4> csolmp, csolREFmp;
  static const mpcmplx I = mpcmplx(0,1);
  //srand48(time(NULL));
  srand48(0);
  int i, fine = 0, k1;
  int maxtrials = (argc > 1)?atoi(argv[1]):100000;
  for (i=0; i < maxtrials && !fine; i++)
    {
      A = mpreal((drand48()-0.5)*1);
      B = mpreal((drand48()-0.5)*1);
      C = mpreal((drand48()-0.5)*1);
      D = mpreal((drand48()-0.5)*1);
    
      //printf("CASE 26\n");

      cpvmp[4]=mpreal(1.0);
      cpvmp[3]=C;
      cpvmp[2]=A;
      cpvmp[1]=D;
      cpvmp[0]=B;

      Qmp.set_coeff(cpvmp);

      Qmp.set_check_always_d20(true);
      Qmp.find_roots(rmp);
#if 0
      csolREFmp = rmp;
      tmpfact = Qmp.maxfact;
      Qmp.set_check_always_d20(false);
      //cout << "second quar\n";
      Qmp.find_roots(rmp);
      csolmp = rmp; 

      sort_sol_opt(csolmp, csolREFmp);

      err = 0.0;
      for (k1=0; k1 < 4; k1++)
        {
          err += abs(csolREFmp[k1])!=mpreal(0.0)?(abs(csolmp[k1]-csolREFmp[k1])/abs(csolREFmp[k1])):abs(csolmp[k1]-csolREFmp[k1]);
        }
#endif
#if 0
      cout << setprecision(30) << "err=" << err << "\n";
      Qmp.show("poly= ");
      for (k1=0; k1 < 4; k1++)
        {
          cout << "sol#" << k1 << setprecision(WP) << " " <<  csolmp[k1] << "\n";
        }
#endif
#if 0
      if (err > 0.1)
        {
          if (tmpfact > epsfact)
            epsfact = tmpfact;
#if 0
          if (tmpfact > 0.1)
            {
              printf("BAD POLYNOMIAL\n");
              //printf("x^4 + (%.16G)*x^3 + (%.16G)*x^2 + (%.16G)*x + %.16G\n", c[3], c[2], c[1], c[0]);
              cout << setprecision(WP) << "x^4+(" << cpvmp[3] << ")*x^3 + (" << cpvmp[2] << ")*x^2 + (" << cpvmp[1] << ")*x + (" <<  cpvmp[0] << ") == 0\n";
              print_roots("found roots", csolmp[0], csolmp[1], csolmp[2], csolmp[3]);
              printf("\n");
              print_roots("exact roots", csolREFmp[0], csolREFmp[1], csolREFmp[2], csolREFmp[3]);
              //printf("err=%16G\n", err);
              cout << "tmpfact=" << tmpfact << " err=" << err << "\n";

              exit(-1);
            }
#endif
          fine = 0;
        }
#endif
    }
  //cout << setprecision(20) << "maxfact=" <<  epsfact << "\n";;
  exit(-1);
}
