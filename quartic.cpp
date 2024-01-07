#include"./quartic.hpp"
using dtype=long double;
int main(void)
{
  quartic<dtype,complex<dtype>,true> Q;
  pvector<dtype,-1> c(5);
  pvector<complex<dtype>,-1> r(4);
  // This is the Case 2 among accuracy tests in ACM Trans. Math. Softw. 46, 2, Article 20 (May 2020),
  // https://doi.org/10.1145/3386241
  c << 16.048044012,-32.072044006,24.036011,-8.006,1.0; 
  //c << 1.0,-17.678187643398402,-14.480985471938862,2.2459773428819827,1.0;
  // polynomial reported by tyuuni by github
  //c << 10147.0, 0.0, -30266.0, 0.0, -40413.0; 
  //c << 17654.0, 0.0, -15252.0, 0.0, -32906.0;
  //Q.set_check_always_d20(true);
  Q.set_coeff(c);
  Q.show("p(x)="); 
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
