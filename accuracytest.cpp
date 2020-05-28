#define WP 200
#include <boost/multiprecision/mpc.hpp>
#include <boost/multiprecision/mpfr.hpp>
using namespace boost;
using namespace boost::multiprecision;
using namespace boost::multiprecision::backends;
//we set 100 digits working precision!
using mpreal=number<mpfr_float_backend<WP>>;
using mpcmplx=number<mpc_complex_backend<WP>>;
#include"./quartic.hpp"

double print_accuracy_at(char *str, double complex *csol, complex double exsol[4])
{
  /* we follow FLocke here */
  int k1;
  double relerr, relerrmax;
  for (k1=0; k1 < 4; k1++)
    {
      relerr=cabs((csol[k1] - exsol[k1])/exsol[k1]); 
      if (k1==0 || relerr > relerrmax)
        {
          relerrmax=cabs((csol[k1] - exsol[k1])/exsol[k1]); 
        }
    }
  printf("[%s] relative accuracy=%.16G\n", str, relerrmax);
  return relerrmax;
}
void print_roots(char *str, mpcmplx x1c, mpcmplx x2c,
                 mpcmplx x3c, mpcmplx x4c)
{
  cout << str;
  cout << "exact root #1=" << x1c.real() << "+I*(" << x1c.imag() <<  ")\n";
  cout << "exact root #2=" << x2c.real() << "+I*(" << x2c.imag() <<  ")\n";
  cout << "exact root #3=" << x3c.real() << "+I*(" << x3c.imag() <<  ")\n";
  cout << "exact root #4=" << x4c.real() << "+I*(" << x4c.imag() <<  ")\n";
}

int main(void)
{
  quartic<double> Q;
  pvector<double,5> cpv;
  pvector<complex<double>,4> r;
  mpcmplx x1c, x2c, x3c, x4c; 
  cmplex<double> csol[4];
  complex<double> csolREF[4];
  mpreal c[5], S;
  int num, k1, okHQR, caso;
  if (argc == 2)
    {
      caso = atoi(argv[1]);
    }
  else
    {
      caso = 1;
    }
  if (caso <= 0 || caso > 24)
    {
      printf("caso must be between 1 and 20\n");
      exit(-1);
    }
  x1c=x2c=x3c=x4c=0.0;
  if (caso > 0)
    {
      switch (caso)
        {
          /* Strobach */
        case 14:
          x1c = x2c = x3c = x4c = 1000;
          print_roots("CASE 14", x1c, x2c, x3c, x4c);
          break;
        case 15:
          x1c = x2c = x3c = 1000;
          x4c = 1E-15;
          print_roots("CASE 15", x1c, x2c, x3c, x4c);
          break;
        case 16:
          x3c = 1E16 + I*1E7;
          x4c = 1E16 - I*1E7;
          x1c = 1 + 0.1*I;
          x2c = 1 - 0.1*I;
          print_roots("CASE 16", x1c, x2c, x3c, x4c);
          break;
        case 17:
          x1c=10000;
          x2c=10001;
          x3c=10010;
          x4c=10100;
          print_roots("CASE 17", x1c, x2c, x3c, x4c);
          break;
        case 18:
          x1c=4E5+I*3E2;
          x2c=4E5-I*3E2;
          x3c=3E4+I*7E3;
          x4c=3E4-I*7E3;
          print_roots("CASE 18", x1c, x2c, x3c, x4c);
          break;
        case 1:
          x1c = 1E9;
          x2c = 1E6;
          x3c = 1E3;
          x4c = 1;
          print_roots("CASE 1", x1c, x2c, x3c, x4c);
          break;
        case 2:
          x1c = 2.003;
          x2c = 2.002;
          x3c = 2.001;
          x4c = 2;
          print_roots("CASE 2", x1c, x2c, x3c, x4c);
          break;
        case 3:
          x1c = 1E53;
          x2c = 1E50;
          x3c = 1E49;
          x4c = 1E47;
          print_roots("CASE 3", x1c, x2c, x3c, x4c);
          break;
        case 4:
          x1c = 1E14;
          x2c = 2;
          x3c = 1;
          x4c = -1;
          print_roots("CASE 4", x1c, x2c, x3c, x4c);
          break;
        case 5:
          x1c = -2E7;
          x2c = 1E7;
          x3c = 1;
          x4c = -1;
          print_roots("CASE 5", x1c, x2c, x3c, x4c);
          break;
        case 6:
          x1c = 1E7;
          x2c = -1E6;
          x3c = 1+I;
          x4c = 1-I;
          print_roots("CASE 6", x1c, x2c, x3c, x4c);
          break;
        case 7:
          x1c = -7;
          x2c = -4;
          x3c = -1E6+I*1E5;
          x4c = -1E6-I*1E5;
          print_roots("CASE 7", x1c, x2c, x3c, x4c);
          break;
        case 8:
          x1c = 1E8;
          x2c = 11;
          x3c = 1E3+I;
          x4c = 1E3-I;
          print_roots("CASE 8", x1c, x2c, x3c, x4c);
          break;
        case 9:
          x1c = 1E7+I*1E6;
          x2c = 1E7-I*1E6;
          x3c = 1+2*I;
          x4c = 1-2*I;
          print_roots("CASE 9", x1c, x2c, x3c, x4c);
          break;
        case 10:
          x1c = 1E4+3*I;
          x2c = 1E4-3*I;
          x3c = -7+1E3*I;
          x4c = -7-1E3*I;
          print_roots("CASE 10", x1c, x2c, x3c, x4c);
          break;
        case 11:
          x1c = 1.001+4.998*I;
          x2c = 1.001-4.998*I;
          x3c = 1.000+5.001*I;
          x4c = 1.000-5.001*I;
          print_roots("CASE 11", x1c, x2c, x3c, x4c);
          break;
        case 12:
          x1c = 1E3+3*I;
          x2c = 1E3-3*I;
          x3c = 1E3+I;
          x4c = 1E3-I;
          print_roots("CASE 12", x1c, x2c, x3c, x4c);
          break;
        case 13:
          x1c = 2+1E4*I;
          x2c = 2-1E4*I;
          x3c = 1+1E3*I;
          x4c = 1-1E3*I;
          print_roots("CASE 13", x1c, x2c, x3c, x4c);
          break;
        case 19:
          x1c = 1E44;
          x2c = 1E30;
          x3c = 1E30;
          x4c = 1.0;
          print_roots("CASE 19", x1c, x2c, x3c, x4c);
          break;
        case 20:
          x1c = 1E14;
          x2c = 1E7;
          x3c = 1E7;
          x4c = 1.0;
          print_roots("CASE 20", x1c, x2c, x3c, x4c);
          break;
        case 21:
          x1c = 1E15;
          x2c = 1E7;
          x3c = 1E7;
          x4c = 1.0;
          print_roots("CASE 21", x1c, x2c, x3c, x4c);
          break;
        case 22:
          x1c = 1E154;
          x2c = 1E152;
          x3c = 10.0;
          x4c = 1.0;
          print_roots("CASE 22", x1c, x2c, x3c, x4c);
          break;
        case 23:
          /* condition */
          c[4] = 1.0;
          c[3] = 1.0;
          c[2] = 1.0;
          c[1] = 3.0/8.0;
          c[0] = 0.001;
          printf("CASE 23\n");
          break;
        case 24:
          S = 1E30;
          c[4] = 1.0;
          c[3] = -(1.0+1.0/S);
          c[2] = 1.0/S - S*S;
          c[1] = S*S + S;
          c[0] = -S;
          printf("CASE 24\n");
          break;
        }
    }
  if (caso <=22)
    {
      csolREF[0]=x1c;
      csolREF[1]=x2c;
      csolREF[2]=x3c;
      csolREF[3]=x4c;
      c[4] = 1.0;
      c[3] = mpcmplx(-(x1c+x2c+x3c+x4c)).real();
      c[2] = mpcmplx(x1c*x2c + (x1c+x2c)*(x3c+x4c) + x3c*x4c).real(); 
      c[1] = mpcmplx(-x1c*x2c*(x3c+x4c) - x3c*x4c*(x1c+x2c)).real();
      c[0] = mccmplx(x1c*x2c*x3c*x4c).real();
    }
  else
    csolREF[0]=csolREF[1]=csolREF[2]=csolREF[3]=0.0;
  cout << "(" << c[4] << ")*x^4+(" << c[3] << ")*x^3+(" << c[2] << ")*x^2+(" << c[1] << ")*x+(" << c[0]<< ")==0\n";

  cpv << c[0], c[1], c[2], c[3], c[3];

  Q.set_coeff(cpv);
  Q.find_roots(r);

  sort_sol_opt(csol, csolREF);
  for (k1=0; k1 < 4; k1++)
    {
      if (caso <= 22)
        printf("[ODM] root #%d=  %.15G+I*(%.15G) [%.15G + I*(%.15G)]\n", 
               k1, creal(csol[k1]), cimag(csol[k1]), creal(csolREF[k1]), cimag(csolREF[k1]));
      else
        printf("[ODM] root #%d=  %.15G+I*(%.15G)\n", 
               k1, creal(csol[k1]), cimag(csol[k1]));
    }
  if (caso <=22)
    print_accuracy_at("ODM", csol, csolREF);

  exit(-1);
} 
