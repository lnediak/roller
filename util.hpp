#ifndef ROLLER_UTIL_HPP_
#define ROLLER_UTIL_HPP_

#include "vector.hpp"

namespace roller {

struct SliceDirs {

  v::DVec<3> c;       /// camera
  v::DVec<3> r, u, f; /// right, up, forward
  double fm;          /// forward multiplier (render distance)
  double rm, um;      /// r and u multipliers (at f dist. 1) (aspect ratio+fov)
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

