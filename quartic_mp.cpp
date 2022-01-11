#define WP 80
// N.B. you can use either CPP, GMP or MPC backend by
// defining CPP_MP, GMP_MP or MPC_MP
#define MPC_MP
#ifdef CPP_MP
#include <boost/multiprecision/cpp_bin_float.hpp> 
#include <boost/multiprecision/cpp_complex.hpp>
using namespace boost;
using namespace boost::multiprecision;
using namespace boost::multiprecision::backends;
using mpreal = number<cpp_bin_float<WP>>;
using mpcmplx = cpp_complex<WP>;
#elif defined(GMP_MP)
#include <boost/multiprecision/gmp.hpp>
#include <boost/multiprecision/complex_adaptor.hpp>
using namespace boost;
using namespace boost::multiprecision;
using namespace boost::multiprecision::backends;
using mpreal=number<gmp_float<WP>>;
using mpcmplx=number<complex_adaptor<gmp_float<WP>>>;
#elif defined(MPC_MP)
#include <boost/multiprecision/mpc.hpp>
#include <boost/multiprecision/mpfr.hpp>
using namespace boost;
using namespace boost::multiprecision;
using namespace boost::multiprecision::backends;
//we set 100 digits working precision!
using mpreal=number<mpfr_float_backend<WP>>;
using mpcmplx=number<mpc_complex_backend<WP>>;
#endif

#include "./quartic.hpp"
int main(void)
{
  quartic<mpreal,mpcmplx> Q;
  pvector<mpreal,5> c;
  pvector<mpcmplx,4> r;
  // This is the Case 2 among accuracy tests in ACM Trans. Math. Softw. 46, 2, Article 20 (May 2020),
  // https://doi.org/10.1145/3386241
  c << mpreal("16.048044012"),mpreal("-32.072044006"),mpreal("24.036011"),mpreal("-8.006000000000000"),mpreal("1.0"); 
  //c << mpreal("1E-06"), mpreal("-1000000000.000000003"), mpreal("3000000.000000000003"), mpreal("-3000.000000000000001"), mpreal("1.0"); 
  // 
  // four multiple roots equal to 1E20 
  //c << mpreal("1E80"), mpreal("-4E60"), mpreal("6E40"), mpreal("-4E20"), mpreal("1.0");
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
