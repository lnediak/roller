#ifndef ROLLER_UTIL_HPP_
#define ROLLER_UTIL_HPP_

#include <cmath>

#include "vector.hpp"

namespace roller {

struct SliceDirs {

  v::DVec<3> c;       /// camera
  v::DVec<3> r, u, f; /// right, up, forward
  double rm, um;      /// r and u multipliers (at f dist. 1) (aspect ratio+fov)
  double fm;          /// forward multiplier (render distance)
};

struct DMat3x3 {
  v::DVec<3> a, b, c; /// rows

  DMat3x3 &operator*=(double d) {
    a *= d;
    b *= d;
    c *= d;
    return *this;
  }

  v::DVec<3> operator*(const v::DVec<3> &v) const {
    return {v::dot(a, v), v::dot(b, v), v::dot(c, v)};
  }

  DMat3x3 operator*(DMat3x3 o) const {
    o = o.transpose();
    DMat3x3 ret = {*this * o.a, *this * o.b, *this * o.c};
    return ret.transpose();
  }

  DMat3x3 transpose() const {
    return {{a[0], b[0], c[0]}, {a[1], b[1], c[1]}, {a[2], b[2], c[2]}};
  }

  double det() const {
    return a[0] * b[1] * c[2] + a[1] * b[2] * c[0] + a[2] * b[0] * c[1] -
           a[0] * b[2] * c[1] - a[1] * b[0] * c[2] - a[2] * b[1] * c[0];
  }

  DMat3x3 adjugate() const {
    return {{b[1] * c[2] - b[2] * c[1], c[1] * a[2] - c[2] * a[1],
             a[1] * b[2] - a[2] * b[1]},
            {b[2] * c[0] - b[0] * c[2], c[2] * a[0] - c[0] * a[2],
             a[2] * b[0] - a[0] * b[2]},
            {b[0] * c[1] - b[1] * c[0], c[0] * a[1] - c[1] * a[0],
             a[0] * b[1] - a[1] * b[0]}};
  }

  DMat3x3 inverse() const { return adjugate() *= (1 / det()); }
};

v::DVec<16> genProjMat(const SliceDirs &sd) {
  DMat3x3 viewspan = {sd.r, sd.u, sd.f};
  DMat3x3 toV = viewspan.transpose().inverse();
  v::DVec<3> mC = -(toV * sd.c);

  DMat3x3 a = {
      {1 / sd.rm, 0, 0}, {0, 1 / sd.um, 0}, {0, 0, (sd.fm + 1) / (sd.fm - 1)}};
  v::DVec<3> rh{0, 0, -2 * sd.fm / (sd.fm - 1)};

  DMat3x3 aa = a * toV;
  v::DVec<3> mmC = a * mC + rh;
  // clang-format off
  return { aa.a[0],  aa.a[1],  aa.a[2], mmC[0],
           aa.b[0],  aa.b[1],  aa.b[2], mmC[1],
           aa.c[0],  aa.c[1],  aa.c[2], mmC[2],
          toV.c[0], toV.c[1], toV.c[2],  mC[2]};
  // clang-format on
}

// ---------------------- ACTUAL UTILITY FUNCTIONS BELOW -----------------------

#define GETFLOOR_ENDIANESS_INDEX 0

/// do not use values out of the range of a short
int fastFloor(float val) {
  val -= 0.5;
  val += 8388608 * 1.5;
  return ((short *)&val)[GETFLOOR_ENDIANESS_INDEX];
}

int fastFloorD(double val) {
  val -= 0.5;
  val += 6755399441055744;
  return ((int *)&val)[GETFLOOR_ENDIANESS_INDEX];
}

/// This is basically floor(val - small_number)
int harshFloor(float val) {
  val -= 0.50001;
  val += 8388608 * 1.5;
  return ((short *)&val)[GETFLOOR_ENDIANESS_INDEX];
}

/// do not use values out of the range of a short
int fastRound(float val) {
  val += 8388608 * 1.5;
  return ((short *)&val)[GETFLOOR_ENDIANESS_INDEX];
}

/// basically only for sdOptimize
double safeAbsR(double a, double b) {
  if (a == 0 && b == 0) {
    return 0;
  }
  return std::abs(a / b);
}
/// minimizing c1x+c2y+c3z under:
/// z>=0,-rm<=x/z<=rm,-um<=y/z<=um,x^2+y^2+z^2<=fm^2
double sdOptimize(double c1, double c2, double c3, double rm, double um,
                  double fm) {
  double curMin = INFINITY;
  // applying lagrange multipliers:
  // 2lx=c1,2ly=c2,2lz=c3,x^2+y^2+z^2=fm^2
  // r=1/2l implies obj=r(c1^2+c2^2+c3^2)
  // while fm^2=x^2+y^2+z^2=r^2(c1^2+c2^2+c3^2)
  double nor = c1 * c1 + c2 * c2 + c3 * c3;
  double fm2 = fm * fm;
  double fatr = -std::sqrt(fm2 / nor);
  /* std::cout << "x,y,z(sphere): " << fatr * c1 << " " << fatr * c2 << " "
            << fatr * c3 << std::endl; */
  if (c3 * fatr >= 0 && safeAbsR(c1, c3) <= rm && safeAbsR(c2, c3) <= um) {
    curMin = fatr * nor;
  }

  // great circles: x=-zrm,x=zrm,y=-zum,y=zum
  // take a given circle, say x=zrm. then
  // it can be written in just y and z.
  // obj=c1rmz+c2y+c3z, so z's coef is cz=(c1rm+c3).
  // x^2+y^2+z^2=y^2+z^2(1+rm^2)=fm^2.
  // let z0=z(1+rm^2)^0.5. cz0=cz/(1+rm^2)^0.5.
  // applying lagrange multipliers:
  // 2ly=c2,2lz0=cz0,y^2+z0^2=fm^2
  // r=1/2l implies obj=r(c2^2+cz0^2)
  // while fm^2=y^2+z0^2=r^2(c2^2+cz0^2)
  double rm12 = 1 + rm * rm;
  double cz = c3 + c1 * rm;
  nor = c2 * c2 + cz * cz / rm12;
  fatr = -std::sqrt(fm2 / nor);
  double c3mod = cz / rm12;
  /* std::cout << "x,y,z(great circle): " << fatr * c3mod / rm << " "
            << fatr * c2 << " " << fatr * c3mod << std::endl; */
  if (cz * fatr >= 0 && safeAbsR(c2, c3mod) <= um) {
    curMin = std::min(curMin, fatr * nor);
  }
  cz = c3 - c1 * rm;
  nor = c2 * c2 + cz * cz / rm12;
  fatr = -std::sqrt(fm2 / nor);
  c3mod = cz / rm12;
  /* std::cout << "x,y,z(great circle): " << fatr * c3mod / rm << " "
            << fatr * c2 << " " << fatr * c3mod << std::endl; */
  if (cz * fatr >= 0 && safeAbsR(c2, c3mod) <= um) {
    curMin = std::min(curMin, fatr * nor);
  }
  double um12 = 1 + um * um;
  cz = c3 + c2 * um;
  nor = c1 * c1 + cz * cz / um12;
  fatr = -std::sqrt(fm2 / nor);
  c3mod = cz / um12;
  /* std::cout << "x,y,z(great circle): " << fatr * c1 << " "
            << fatr * c3mod / um << " " << fatr * c3mod << std::endl; */
  if (cz * fatr >= 0 && safeAbsR(c1, c3mod) <= rm) {
    curMin = std::min(curMin, fatr * nor);
  }
  cz = c3 - c2 * um;
  nor = c1 * c1 + cz * cz / um12;
  fatr = -std::sqrt(fm2 / nor);
  c3mod = cz / um12;
  /* std::cout << "x,y,z(great circle): " << fatr * c1 << " "
            << fatr * c3mod / um << " " << fatr * c3mod << std::endl; */
  if (cz * fatr >= 0 && safeAbsR(c1, c3mod) <= rm) {
    curMin = std::min(curMin, fatr * nor);
  }

  // finally, we can do the basic version
  double z = fm;
  double x = rm * fm;
  double y = um * fm;
  double mod = std::sqrt(x * x + y * y + z * z) / fm;
  z /= mod;
  x /= mod;
  y /= mod;
  double pru = c1 * x + c2 * y + c3 * z;
  double pr0u = c1 * x - c2 * y + c3 * z;
  double p0ru = -c1 * x + c2 * y + c3 * z;
  double p0r0u = -c1 * x - c2 * y + c3 * z;
  return std::min(
      curMin,
      std::min(0., std::min(pru, std::min(pr0u, std::min(p0ru, p0r0u)))));
}

v::IVec<3> sdMins(const SliceDirs &sd) {
  return {
      fastFloorD(sd.c[0] +
                 sdOptimize(sd.r[0], sd.u[0], sd.f[0], sd.rm, sd.um, sd.fm) -
                 0.0001),
      fastFloorD(sd.c[0] +
                 sdOptimize(sd.r[1], sd.u[1], sd.f[1], sd.rm, sd.um, sd.fm) -
                 0.0001),
      fastFloorD(sd.c[0] +
                 sdOptimize(sd.r[2], sd.u[2], sd.f[2], sd.rm, sd.um, sd.fm) -
                 0.0001)};
}

v::IVec<3> sdMaxs(const SliceDirs &sd) {
  return {
      fastFloorD(sd.c[0] -
                 sdOptimize(-sd.r[0], -sd.u[0], -sd.f[0], sd.rm, sd.um, sd.fm) +
                 1.0001),
      fastFloorD(sd.c[1] -
                 sdOptimize(-sd.r[1], -sd.u[1], -sd.f[1], sd.rm, sd.um, sd.fm) +
                 1.0001),
      fastFloorD(sd.c[2] -
                 sdOptimize(-sd.r[2], -sd.u[2], -sd.f[2], sd.rm, sd.um, sd.fm) +
                 1.0001)};
}

// ---------- NOW FOR SOME STUFF TO SUPPLEMENT VVEC ----------

template <class T> struct FloorU {
  static int apply(T a) {
    int i = a;
    return i - (a < i);
  }
};
struct FastFloorU {
  static int apply(float a) { return fastFloor(a); }
};
struct FastFloorDU {
  static int apply(double a) { return fastFloorD(a); }
};
struct FastRoundU {
  static int apply(float a) { return fastRound(a); }
};
struct Int2FloatU {
  static float apply(int a) { return a; }
};
struct Int2DoubleU {
  static double apply(int a) { return a; }
};

template <class T> struct AbsU {
  static T apply(T a) { return a > 0 ? a : -a; }
};
template <class T> struct L1NormU {
  static T apply(T a, T b) { return a + AbsU<T>::apply(b); }
};

template <class T> struct EI {
  T a;
  std::size_t b;

  EI<T> get(const EI<T> &c) const {
    if (c.a < a) {
      return c;
    }
    return *this;
  }
};
template <class T, std::size_t N, class A> struct MinI {

  typedef MinI<T, N - 1, A> smaller;
  const A &a;
  EI<T> evaluate() const {
    return smaller{a}.evaluate().get({a[N - 1], N - 1});
  }
};
template <class T, class A> struct MinI<T, 1, A> {
  const A &a;
  EI<T> evaluate() const { return {a[0], 0}; }
};

template <class T, std::size_t N, class A>
using Floor = v::UnaryOp<int, N, FloorU<T>, A>;
template <std::size_t N, class A>
using FastFloor = v::UnaryOp<int, N, FastFloorU, A>;
template <std::size_t N, class A>
using FastFloorD = v::UnaryOp<int, N, FastFloorDU, A>;
template <std::size_t N, class A>
using FastRound = v::UnaryOp<int, N, FastRoundU, A>;
template <std::size_t N, class A>
using Int2Float = v::UnaryOp<float, N, Int2FloatU, A>;
template <std::size_t N, class A>
using Int2Double = v::UnaryOp<double, N, Int2DoubleU, A>;

template <class T, std::size_t N, class A>
using Abs = v::UnaryOp<T, N, AbsU<T>, A>;
template <class T, std::size_t N, class A>
using L1Norm = v::Combine<T, N, L1NormU<T>, AbsU<T>, A>;

template <class A, class = typename rem_cvr<A>::thisisavvec>
Floor<typename A::value_type, A::size, A> vfloor(const A &a) {
  return {a};
}
template <class A, class = typename rem_cvr<A>::thisisavvec,
          class = typename std::enable_if<
              std::is_same<typename A::value_type, float>::value>::type>
FastFloor<A::size, A> vfastFloor(const A &a) {
  return {a};
}
template <class A, class = typename rem_cvr<A>::thisisavvec,
          class = typename std::enable_if<
              std::is_same<typename A::value_type, double>::value>::type>
FastFloorD<A::size, A> vfastFloorD(const A &a) {
  return {a};
}
template <class A, class = typename rem_cvr<A>::thisisavvec,
          class = typename std::enable_if<
              std::is_same<typename A::value_type, float>::value>::type>
FastRound<A::size, A> vfastRound(const A &a) {
  return {a};
}
template <class A, class = typename rem_cvr<A>::thisisavvec,
          class = typename std::enable_if<
              std::is_same<typename A::value_type, int>::value>::type>
Int2Float<A::size, A> vint2float(const A &a) {
  return {a};
}
template <class A, class = typename rem_cvr<A>::thisisavvec,
          class = typename std::enable_if<
              std::is_same<typename A::value_type, int>::value>::type>
Int2Double<A::size, A> vint2double(const A &a) {
  return {a};
}

template <class A, class = typename rem_cvr<A>::thisisavvec>
EI<typename A::value_type> vminI(const A &a) {
  return MinI<typename A::value_type, A::size, A>{a}.evaluate();
}

template <class A, class = typename rem_cvr<A>::thisisavvec>
typename A::value_type vabs(const A &a) {
  return Abs<typename A::value_type, A::size, A>(a).evaluate();
}
template <class A, class = typename rem_cvr<A>::thisisavvec>
typename A::value_type l1norm(const A &a) {
  return L1Norm<typename A::value_type, A::size, A>(a).evaluate();
}

// ------------ LERPING!!!!!! -------------

template <std::size_t N> struct BitsVec {
  typedef void thisisavvec;
  typedef std::int32_t value_type;
  static const std::size_t size = N;

  std::size_t value;

  value_type operator[](std::size_t i) const { return (value & (1 << i)) >> i; }
};
template <std::size_t N> struct SplitVec {
  typedef void thisisavvec;
  typedef double value_type;
  static const std::size_t size = N;

  const v::IVec<N> &inds;
  const double *vecs[2];

  value_type operator[](std::size_t i) const { return vecs[inds[i]][i]; }
};
template <std::size_t N, std::size_t M, std::size_t I, class A, class B>
struct Lerper {
  typedef void thisisavvec;
  typedef double value_type;
  static const std::size_t size = 1 << (M - 1);

  A tolerp;
  B lerp;

  value_type operator[](std::size_t i) const {
    double x = tolerp[i];
    double y = tolerp[i + size];
    return x + lerp[I] * (y - x);
  }
};
template <std::size_t N, std::size_t M, std::size_t I, class A, class B>
struct LerperT {
  typedef LerperT<N, M + 1, I + 1, A, B> nextlerpert;
  typedef Lerper<N, M, I, typename nextlerpert::nextlerper, const B &>
      nextlerper;

  nextlerper operator()(A tolerp, const B &lerp) const {
    return {nextlerpert{}(tolerp, lerp), lerp};
  }
};
template <std::size_t N, std::size_t I, class A, class B>
struct LerperT<N, N, I, A, B> {
  typedef Lerper<N, N, I, A, B> nextlerper;

  nextlerper operator()(A tolerp, const B &lerp) const {
    return {tolerp, lerp};
  }
};

/// lerp[0] applies for least sigdig of index in vecpow, etc
template <
    class A, class B, class = typename A::thisisavvec,
    class = typename B::thisisavvec,
    class = typename std::enable_if<
        std::is_same<typename A::value_type, typename B::value_type>::value &&
        A::size == (1 << B::size)>::type>
typename A::value_type vlerpAll(const A &vecpow, const B &lerp) {
  return LerperT<B::size, 1, 0, const A &, B>{}(vecpow, lerp)[0];
}

} // namespace roller

#endif // ROLLER_UTIL_HPP_

