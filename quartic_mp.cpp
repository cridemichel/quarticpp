#define WP 50
#include <boost/multiprecision/mpc.hpp>
#include <boost/multiprecision/mpfr.hpp>
using namespace boost;
using namespace boost::multiprecision;
using namespace boost::multiprecision::backends;
//we set 100 digits working precision!
using mpdbl=number<mpfr_float_backend<WP>>;
using mpcmplx=number<mpc_complex_backend<WP>>;
#include "./quartic.hpp"
int main(void)
{
  quartic<mpdbl,mpcmplx> Q;
  pvector<mpdbl,5> c;
  pvector<mpcmplx,4> r;
  c << 24,-50,35,-10,1.1;
  Q.set_coeff(c);
  Q.find_roots(r);
  r.set_show_digits(WP);
  r.show("roots");
  cout << setprecision(WP) << "OQS p(1.644551988810059402960370803461652927757816938645)=" 
    << Q.evalpoly(mpdbl("1.644551988810059402960370803461652927757816938645")) << "\n";
  cout << setprecision(WP) << "MPSOLVE p(1.6445519888100595858459604805012253420901966931458908262601)=" << Q.evalpoly(mpdbl("1.6445519888100595858459604805012253420901966931458908262601")) << "\n";
  return 0;
}
