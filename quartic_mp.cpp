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
  return 0;
}
