#include"./quartic.hpp"
int main(void)
{
  quartic<double> Q;
  pvector<double,5> c;
  pvector<complex<double>,4> r;
  c << 24,-50,35,-10,1.1;
  Q.set_coeff(c);
  Q.find_roots(r);
  r.show("roots");
  return 0;
}
