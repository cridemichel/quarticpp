#define WP 50
#include <boost/multiprecision/mpc.hpp>
#include <boost/multiprecision/mpfr.hpp>
using namespace boost;
using namespace boost::multiprecision;
using namespace boost::multiprecision::backends;
//we set 100 digits working precision!
using mpreal=number<mpfr_float_backend<WP>>;
using mpcmplx=number<mpc_complex_backend<WP>>;
#include "./quartic.hpp"
int main(void)
{
  quartic<mpreal,mpcmplx> Q;
  pvector<mpcmplx,5> c;
  pvector<mpcmplx,4> r;
  mpcmplx x1c, x2c, x3c, x4c;
  // This is the Case 2 among accuracy tests in ACM Trans. Math. Softw. 46, 2, Article 20 (May 2020),
  // https://doi.org/10.1145/3386241
  x3c = mpcmplx("10.0", "1.0");
  x4c = mpcmplx("1.0", "10");
  x1c = mpcmplx("-1.0", "0.1");
  x2c = mpcmplx("2.0", "-0.001");
  c[4] = 1.0;
  c[3] = -(x1c+x2c+x3c+x4c);
  c[2] = x1c*x2c + (x1c+x2c)*(x3c+x4c) + x3c*x4c; 
  c[1] = -x1c*x2c*(x3c+x4c) - x3c*x4c*(x1c+x2c);
  c[0] = x1c*x2c*x3c*x4c;

  Q.set_coeff(c);
  Q.find_roots(r);
  //r.show("roots");
  int cc=0;
  for (auto& r0: r)
    {
      cout << setprecision(WP) << "root #" << cc <<  "=" << r0 << "\n";
      cout << setprecision(WP) << "p(#" << cc << ")=" << Q.evalpoly(r0) << "\n\n";
      cc++;
    }
  cout << "p=" << Q.evalpoly(mpreal(1.0)) << "\n";
  //cout << setprecision(WP) << "MPSOLVE p(4.6356013430199667501755766197408175062413493523496385264695)=" << Q.evalpoly(mpreal("4.6356013430199667501755766197408175062413493523496385264695")) << "\n";
  return 0;
}
