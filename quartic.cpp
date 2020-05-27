#include"./quartic.hpp"
int main(void)
{
  quartic<double> Q;
  pvector<double,5> c;
  pvector<complex<double>,4> r;
  c << 16.048044012,-32.072044006,24.036011,-8.006,1.0; 
  Q.set_coeff(c);
 
  Q.find_roots(r);
  r.show("roots");
  int cc=0;
  for (auto& r0: r)
    {
      cout << setprecision(16) << "root #" << cc <<  "=" << r0 << "\n";
      cout << setprecision(16) << "p(#" << cc << ")=" << Q.evalpoly(r0) << "\n\n";
      cc++;
    }

  return 0;
}
