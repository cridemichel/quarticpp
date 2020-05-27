#define OQS_MULTIPLE_PRECISION
#include "./quartic.hpp"
int main(void)
{
  quartic<mpdbl,mpcmplx> Q;
  pvector<mpdbl,5> c;
  pvector<mpcmplx,4> r;
  c << 24,-50,35,-10,1.1;
  Q.set_coeff(c);
  Q.find_roots(r);
  r.set_show_digits(50);
  r.show("roots");
  cout << setprecision(50) << "OQS p(1.644551988810059402960370803461652927757816938645)=" 
    << Q.evalpoly(mpdbl("1.644551988810059402960370803461652927757816938645")) << "\n";
  cout << setprecision(50) << "MPSOLVE p(1.6445519888100595858459604805012253420901966931458908262601)=" << Q.evalpoly(mpdbl("1.6445519888100595858459604805012253420901966931458908262601")) << "\n";
  return 0;
}
