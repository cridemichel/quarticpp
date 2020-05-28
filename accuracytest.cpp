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

int perm[24][4]={
      {0, 1, 2, 3}, {0, 1, 3, 2}, {0, 2, 1, 3}, {0, 2, 3, 1}, {0, 3, 1, 2}, {0, 3, 2, 1}, 
      {1, 0, 2, 3}, {1, 0, 3, 2}, {1, 2, 0, 3}, {1, 2, 3, 0}, {1, 3, 0, 2}, {1, 3, 2, 0}, 
      {2, 0, 1, 3}, {2, 0, 3, 1}, {2, 1, 0, 3}, {2, 1, 3, 0}, {2, 3, 0, 1}, {2, 3, 1, 0},
      {3, 0, 1, 2}, {3, 0, 2, 1}, {3, 1, 0, 2}, {3, 1, 2, 0}, {3, 2, 0, 1}, {3, 2, 1, 0}};

void sort_sol_opt(pvector<mpcmplx,4>& sol, pvector<mpcmplx,4>& exsol)
{
  int k1, k2, k1min;
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
void print_roots(const char *str, mpcmplx x1c, mpcmplx x2c,
                 mpcmplx x3c, mpcmplx x4c)
{
  cout << str;
  cout << "exact root #1=" << x1c.real() << "+I*(" << x1c.imag() <<  ")\n";
  cout << "exact root #2=" << x2c.real() << "+I*(" << x2c.imag() <<  ")\n";
  cout << "exact root #3=" << x3c.real() << "+I*(" << x3c.imag() <<  ")\n";
  cout << "exact root #4=" << x4c.real() << "+I*(" << x4c.imag() <<  ")\n";
}

int main(int argc, char** argv)
{
  quartic<double> Q;
  quartic<mpreal,mpcmplx> Qmp;
  pvector<double,5> cpv;
  pvector<mpreal,5> cpvmp;
  pvector<complex<double>,4> r;
  pvector<mpcmplx,4> rmp;
  mpcmplx x1c, x2c, x3c, x4c; 
  mpreal c[5], S;
  pvector<complex<double>,4> csol, csolREF;
  pvector<mpcmplx,4> csolmp, csolREFmp;
   static const mpcmplx  I = mpcmplx(0,1);
  int k1, caso;
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
  /* This cases are those discussed in  ACM Trans. Math. Softw. 46, 2, Article 20 (May 2020),
   * https://doi.org/10.1145/3386241
   */
  if (caso > 0)
    {
      switch (caso)
        {
          /* Strobach */
        case 14:
          x1c = x2c = x3c = x4c = mpcmplx("1000");
          print_roots("CASE 14", x1c, x2c, x3c, x4c);
          break;
        case 15:
          x1c = x2c = x3c = mpcmplx("1000");
          x4c = mpcmplx("1E-15");
          print_roots("CASE 15", x1c, x2c, x3c, x4c);
          break;
        case 16:
          x3c = mpcmplx("1E16", "1E7");
          x4c = mpcmplx("1E16", "-1E7");
          x1c = mpcmplx("1", "0.1");
          x2c = mpcmplx("1", "-0.1");
          print_roots("CASE 16", x1c, x2c, x3c, x4c);
          break;
        case 17:
          x1c=mpcmplx("10000");
          x2c=mpcmplx("10001");
          x3c=mpcmplx("10010");
          x4c=mpcmplx("10100");
          print_roots("CASE 17", x1c, x2c, x3c, x4c);
          break;
        case 18:
          x1c=mpcmplx("4E5","3E2");
          x2c=mpcmplx("4E5","-3E2");
          x3c=mpcmplx("3E4","7E3");
          x4c=mpcmplx("3E4","-7E3");
          print_roots("CASE 18", x1c, x2c, x3c, x4c);
          break;
        case 1:
          x1c = mpcmplx("1E9");
          x2c = mpcmplx("1E6");
          x3c = mpcmplx("1E3");
          x4c = mpcmplx("1");
          print_roots("CASE 1", x1c, x2c, x3c, x4c);
          break;
        case 2:
          x1c = mpcmplx("2.003");
          x2c = mpcmplx("2.002");
          x3c = mpcmplx("2.001");
          x4c = mpcmplx("2");
          print_roots("CASE 2", x1c, x2c, x3c, x4c);
          break;
        case 3:
          x1c = mpcmplx("1E53");
          x2c = mpcmplx("1E50");
          x3c = mpcmplx("1E49");
          x4c = mpcmplx("1E47");
          print_roots("CASE 3", x1c, x2c, x3c, x4c);
          break;
        case 4:
          x1c = mpcmplx("1E14");
          x2c = mpcmplx("2");
          x3c = mpcmplx("1");
          x4c = mpcmplx("-1");
          print_roots("CASE 4", x1c, x2c, x3c, x4c);
          break;
        case 5:
          x1c = mpcmplx("-2E7");
          x2c = mpcmplx("1E7");
          x3c = mpcmplx("1");
          x4c = mpcmplx("-1");
          print_roots("CASE 5", x1c, x2c, x3c, x4c);
          break;
        case 6:
          x1c = mpcmplx("1E7");
          x2c = mpcmplx("-1E6");
          x3c = mpcmplx("1","1");
          x4c = mpcmplx("1","-1");
          print_roots("CASE 6", x1c, x2c, x3c, x4c);
          break;
        case 7:
          x1c = mpcmplx("-7");
          x2c = mpcmplx("-4");
          x3c = mpcmplx("-1E6","1E5");
          x4c = mpcmplx("-1E6","-1E5");
          print_roots("CASE 7", x1c, x2c, x3c, x4c);
          break;
        case 8:
          x1c = mpcmplx("1E8");
          x2c = mpcmplx("11");
          x3c = mpcmplx("1E3","1");
          x4c = mpcmplx("1E3","-1");
          print_roots("CASE 8", x1c, x2c, x3c, x4c);
          break;
        case 9:
          x1c = mpcmplx("1E7","1E6");
          x2c = mpcmplx("1E7","-1E6");
          x3c = mpcmplx("1","2");
          x4c = mpcmplx("1","-2");
          print_roots("CASE 9", x1c, x2c, x3c, x4c);
          break;
        case 10:
          x1c = mpcmplx("1E4","3");
          x2c = mpcmplx("1E4","-3");
          x3c = mpcmplx("-7","1E3");
          x4c = mpcmplx("-7","-1E3");
          print_roots("CASE 10", x1c, x2c, x3c, x4c);
          break;
        case 11:
          x1c = mpcmplx("1.001","4.998");
          x2c = mpcmplx("1.001","-4.998");
          x3c = mpcmplx("1.000","5.001");
          x4c = mpcmplx("1.000","-5.001");
          print_roots("CASE 11", x1c, x2c, x3c, x4c);
          break;
        case 12:
          x1c = mpcmplx("1E3","3");
          x2c = mpcmplx("1E3","-3");
          x3c = mpcmplx("1E3","1");
          x4c = mpcmplx("1E3","-1");
          print_roots("CASE 12", x1c, x2c, x3c, x4c);
          break;
        case 13:
          x1c = mpcmplx("2","1E4");
          x2c = mpcmplx("2","-1E4");
          x3c = mpcmplx("1","1E3");
          x4c = mpcmplx("1","-1E3");
          print_roots("CASE 13", x1c, x2c, x3c, x4c);
          break;
        case 19:
          x1c = mpcmplx("1E44");
          x2c = mpcmplx("1E30");
          x3c = mpcmplx("1E30");
          x4c = mpcmplx("1.0");
          print_roots("CASE 19", x1c, x2c, x3c, x4c);
          break;
        case 20:
          x1c = mpcmplx("1E14");
          x2c = mpcmplx("1E7");
          x3c = mpcmplx("1E7");
          x4c = mpcmplx("1.0");
          print_roots("CASE 20", x1c, x2c, x3c, x4c);
          break;
        case 21:
          x1c = mpcmplx("1E15");
          x2c = mpcmplx("1E7");
          x3c = mpcmplx("1E7");
          x4c = mpcmplx("1.0");
          print_roots("CASE 21", x1c, x2c, x3c, x4c);
          break;
        case 22:
          x1c = mpcmplx("1E154");
          x2c = mpcmplx("1E152");
          x3c = mpcmplx("10.0");
          x4c = mpcmplx("1.0");
          print_roots("CASE 22", x1c, x2c, x3c, x4c);
          break;
        case 23:
          /* condition */
          c[4] = mpreal("1.0");
          c[3] = mpreal("1.0");
          c[2] = mpreal("1.0");
          c[1] = mpreal("3.0/8.0");
          c[0] = mpreal("0.001");
          printf("CASE 23\n");
          break;
        case 24:
          S = mpreal("1E30");
          c[4] = mpreal("1.0");
          c[3] = -(mpreal("1.0")+1.0/S);
          c[2] = 1.0/S - S*S;
          c[1] = S*S + S;
          c[0] = -S;
          printf("CASE 24\n");
          break;
        }
    }
  if (caso <=22)
    {
      csolREF[0]=complex<double>(x1c);
      csolREF[1]=complex<double>(x2c);
      csolREF[2]=complex<double>(x3c);
      csolREF[3]=complex<double>(x4c);
      csolREFmp[0]=x1c;
      csolREFmp[1]=x2c;
      csolREFmp[2]=x3c;
      csolREFmp[3]=x4c;
      c[4] = 1.0;
      c[3] = mpcmplx(-(x1c+x2c+x3c+x4c)).real();
      c[2] = mpcmplx(x1c*x2c + (x1c+x2c)*(x3c+x4c) + x3c*x4c).real(); 
      c[1] = mpcmplx(-x1c*x2c*(x3c+x4c) - x3c*x4c*(x1c+x2c)).real();
      c[0] = mpcmplx(x1c*x2c*x3c*x4c).real();
    }
  else
    {
      csolREF[0]=csolREF[1]=csolREF[2]=csolREF[3]=0.0;
      csolREF[0]=csolREF[1]=csolREF[2]=csolREF[3]=0.0;
    }
  cout << "(" << c[4] << ")*x^4+(" << c[3] << ")*x^3+(" << c[2] << ")*x^2+(" << c[1] << ")*x+(" << c[0]<< ")==0\n";


  cpv << double(c[0]), double(c[1]), double(c[2]), double(c[3]), double(c[4]);
  cpvmp << c[0],c[1],c[2],c[3],c[4];

  Q.set_coeff(cpv);
  Q.find_roots(r);
  csol = r; 
  sort_sol_opt(csol, csolREF);
  for (k1=0; k1 < 4; k1++)
    {
      if (caso <= 22)
        cout << setprecision(15) << "root #"<< k1 <<  "=  " << csol[k1].real() << "+I*(" << csol[k1].imag() << ") [" << csolREF[k1].real() <<  " + I*(" << csolREF[k1].imag() << ")]\n"; 
      else
        cout << setprecision(15) << "root #"<< k1 <<  "=  " << csol[k1].real() << "+I*(" << csol[k1].imag() << ")\n";
    }
  if (caso <=22)
    print_accuracy_at(csol, csolREF);
  cout << "====================================================\n";
  cout << "Solve same quartic but use 200 digits precision and print roots with 50 digits!\n"; 
  Qmp.set_coeff(cpvmp);
  Qmp.find_roots(rmp);
  csolmp = rmp;
  sort_sol_opt(csolmp, csolREFmp);
  for (k1=0; k1 < 4; k1++)
    {
      if (caso <= 22)
        cout << setprecision(50) << "root #"<< k1 <<  "=  " << csolmp[k1].real() << "+I*(" << csolmp[k1].imag() << ") [" << csolREFmp[k1].real() <<  " + I*(" << csolREFmp[k1].imag() << ")]\n"; 
      else
        cout << setprecision(50) << "root #"<< k1 <<  "=  " << csolmp[k1].real() << "+I*(" << csolmp[k1].imag() << ")\n";
    }
  if (caso <=22)
    print_accuracy_at(csolmp, csolREFmp);

  exit(-1);
} 
