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
  pvector<mpreal,5> c;
  pvector<mpcmplx,4> r;
  c << mpreal("16.048044012"),mpreal("-32.072044006"),mpreal("24.036011"),mpreal("-8.006000000000000"),mpreal("1.0"); 

  Q.set_coeff(c);
  Q.find_roots(r);
  r.show("roots");
  int cc=0;
  for (auto& r0: r)
    {
      cout << setprecision(WP) << "root #" << cc <<  "=" << r0 << "\n";
      cout << setprecision(WP) << "p(#" << cc << ")=" << Q.evalpoly(r0) << "\n\n";
      cc++;
    }
  //cout << setprecision(WP) << "MPSOLVE p(4.6356013430199667501755766197408175062413493523496385264695)=" << Q.evalpoly(mpreal("4.6356013430199667501755766197408175062413493523496385264695")) << "\n";
  return 0;
}
