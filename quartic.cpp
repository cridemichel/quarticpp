#include"./quartic.hpp"
int main(void)
{
  quartic<double> Q;
  pvector<double,5> c;
  pvector<complex<double>,4> r;
  c << 1,1,1,1,1;
  Q.set_coeff(c);
  Q.find_roots(r);
  return 0;
}
