/// Time-Stamp: <2015-01-02 21:40:58 misho>

#include <fstream>
#include <iostream>

#include "../vendor/diagnostic_ignore.hpp"
#include "../vendor/slhaea/slhaea.h"
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/bindings/traits/ublas_symmetric.hpp>
#include <boost/numeric/bindings/lapack/syev.hpp>
#include "../vendor/diagnostic_ignore_end.hpp"

namespace ublas  = boost::numeric::ublas;
namespace lapack = boost::numeric::bindings::lapack;
typedef boost::numeric::ublas::symmetric_matrix<double> sym_matrix;
typedef boost::numeric::ublas::matrix<double>           matrix;
typedef boost::numeric::ublas::vector<double>           vector;

std::ostream& warning(){ std::cerr << "[Warning] "; return std::cerr; }
std::ostream& error  (){ std::cerr << "[ERROR] ";   return std::cerr; }

double get(const SLHAea::Coll& slha, const std::string& block, int key, int loc = 1){
  try {
    return SLHAea::to<double>(slha.at(block).at(key).at(loc));
  }catch(const std::out_of_range& err){
    error() << "Key " << block << " " << key << " not found." << std::endl;
    throw err;
  }
}

class Angle {
  double th_, s_, c_, t_;
  void init(double th) { th_ = th; s_ = sin(th); c_ = cos(th); t_ = s_/c_; }
public:
  //enum TYPE { SIN, COS, TAN };
  Angle() { init(0); };
  Angle(double th) { init(th); }
  //Angle(Angle::TYPE t, double v) {
  //  init(t == SIN ? asin(v) : t == COS ? acos(v) : t == TAN ? atan(v) : throw);
  //}
  Angle(const Angle& a) : th_(a.th_), s_(a.s_), c_(a.c_), t_(a.t_) {};
  Angle& operator=(const Angle& a){
    if(&a != this) { th_ = a.th_; s_ = a.s_; c_ = a.c_; t_ = a.t_; }
    return *this;
  }
  Angle& set(double v) { init(v); return *this; }
  double s() const { return s_; }
  double c() const { return c_; }
  double t() const { return t_; }
  double s2() const { return s_*s_; }
  double c2() const { return c_*c_; }
  double t2() const { return t_*t_; }
};

class SMparams {
  constexpr static const double PI   = 3.14159265358979;
  constexpr static const double sinW = sqrt(0.23116); // we don't use GF, but use this value.
  double ainv, GF, a_s, mz, mb, mt, mtau;
  Angle wangle;
public:
  SMparams(const SLHAea::Coll& coll) {
    try{
      ainv = get(coll, "SMINPUTS", 1); // MS-bar at mZ
    //GF   = get(coll, "SMINPUTS", 2); // *** UNUSED ***
      a_s  = get(coll, "SMINPUTS", 3); // MS-bar at mZ
      mz   = get(coll, "SMINPUTS", 4); // pole
      mb   = get(coll, "SMINPUTS", 5); // MS-bar at mb
      mt   = get(coll, "SMINPUTS", 6); // pole
      mtau = get(coll, "SMINPUTS", 7); // pole
    }catch(const std::out_of_range& er){
      error() << "SMINPUTS not sufficient";
      throw er;
    }
    wangle.set(asin(sinW));
  }
private:
  double ainv_Y_mZ() const { return wangle.c2() * ainv; }
  double ainv_w_mZ() const { return wangle.s2() * ainv; }
  double ainv_3_mZ() const { return 1/a_s; }
  double g_scaled(double ainv_mZ, double scale, double b1) const {
    double ainv = ainv_mZ - (0.5/PI) * b1 * log(scale / mz);
    return sqrt(4*PI/ainv);
  }

public:
  double gY(double scale) const { return g_scaled(ainv_Y_mZ(), scale,  41.0/6.0); }
  double g2(double scale) const { return g_scaled(ainv_w_mZ(), scale, -19.0/6.0); }
  double g3(double scale) const { return g_scaled(ainv_3_mZ(), scale, -7.0); }
  const Angle& w() const { return wangle; }
  double mZ() const { return mz; }
  double mW() const { return mz * wangle.c(); }
};

sym_matrix neutralino_matrix (const SLHAea::Coll& in){
  SMparams sm(in);
  Angle beta(atan(get(in, "MINPAR", 3)));
  double cb = beta.c(), sb = beta.s(), sw = sm.w().s(), cw = sm.w().c(), mZ = sm.mZ();

  sym_matrix m(4);
  m(0,0) = get(in, "EXTPAR", 1);
  m(1,1) = get(in, "EXTPAR", 2);
  m(2,2) = 0;
  m(3,3) = 0;
  m(0,1) = 0;
  m(0,2) = -cb*sw*mZ;
  m(0,3) =  sb*sw*mZ;
  m(1,2) =  cb*cw*mZ;
  m(1,3) = -sb*cw*mZ;
  m(2,3) = -get(in, "EXTPAR", 23);
  return m;
}
#include <boost/static_assert.hpp>
#include <boost/type_traits.hpp>
#include <cassert>

namespace traits = boost::numeric::bindings::traits;

std::pair<vector, matrix> neutralino (const SLHAea::Coll& in){
  matrix m = neutralino_matrix(in);

  vector eigenvalue(4);
  matrix ev = m;
  int info = lapack::syev('V', 'L', ev, eigenvalue);
  BOOST_UBLAS_CHECK(info == 0, ublas::internal_logic());
  

  //  return std::pair<vector, matrix>(v, result);
}

/*

void calculate_treelevel_matrix(const SLHAea::Coll& in){



class SLHA : public SLHAea::Coll {
  using namespace SLHAea;
  double mZ, cw, sw, tb, sb, cb;
public:
  SLHA(std::istream& ifs) : Coll(ifs) { init(); }


  void init(){
    mZ = assert_and_get("SMINPUTS", 4);
    tb = assert_and_get("MINPAR", 3);
    double
      alpha   = assert_and_get("SMINPUTS", 1),
      GF      = assert_and_get("SMINPUTS", 2),
      sin2wSQ = 2*sqrt(2)*3.14159265358979*alpha/GF/mZ/mZ;
    sw = sqrt( (1.0 - sqrt(1.0 - sin2wSQ) ) / 2.0);
    cw = sqrt( (1.0 + sqrt(1.0 - sin2wSQ) ) / 2.0);
    cb = 1.0/sqrt(tb*tb+1);
    sb = cb*tb;
  }

  void neutralino_mass_matrix (SLHAea::Coll& in){
    BSM<double> matrix(4), z(4);
    BM<double> k(4,4);
    double
      M1 = assert_and_get("EXTPAR", 1),
      M2 = assert_and_get("EXTPAR", 2);
    BV<double> v(4), d(6);
    matrix(0,0) = assert_and_get("EXTPAR", 1);
    matrix(1,1) = assert_and_get("EXTPAR", 2);
    matrix(2,2) = 0;
    matrix(3,3) = 0;
    matrix(0,1) = 0;
    matrix(0,2) = -cb*sw*mZ;
    matrix(0,3) =  sb*sw*mZ;
    matrix(1,2) =  cb*cw*mZ;
    matrix(1,3) = -sb*cw*mZ;
    matrix(2,3) = -assert_and_get("EXTPAR", 23);
    v(0) = assert_and_get("EXTPAR", 1);
    v(1) = assert_and_get("EXTPAR", 2);
    v(2) = 0;
    v(3) = 0;
    d(0) = 0;
    d(1) = -cb*sw*mZ;
    d(2) =  sb*sw*mZ;
    d(3) =  cb*cw*mZ;
    d(4) = -sb*cw*mZ;
    d(5) = -assert_and_get("EXTPAR", 23);

      BOOST_UBLAS_CHECK(info == 0, ublas::internal_logic());    }

  };
};

*/

void test_neutralino(const SLHAea::Coll& in){
  matrix neut = neutralino_matrix(in);
  std::pair<vector, matrix> r = neutralino(in);
  matrix m(4,4);
  for (int i=0; i < 4; i++) m(i,i) = r.first[i];
  std::cout << m << std::endl;
  std::cout << neut << std::endl;
  std::cout << ublas::prod(ublas::trans(r.second), neut, r.second) << std::endl;

}


int main(){
  std::ifstream ifs("slha1.txt");
  SLHAea::Coll input(ifs);
  neutralino(input);

  test_neutralino(input);
  
  std::cout << "tan(beta) = "        << input["MINPAR"][3][1] << std::endl << std::endl;
  std::cout << "m_top(pole) line:\n" << input["SMINPUTS"][6]  << std::endl << std::endl;
  std::cout << "SMINPUTS block:\n"   << input["SMINPUTS"];
  //  neutralino_mass_matrix(input);
}


