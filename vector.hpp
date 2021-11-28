#ifndef V_VECTOR_HPP_
#define V_VECTOR_HPP_

#include <cstdint>
#include <initializer_list>
#include <iostream>
#include <type_traits>

template <class T>
using rem_cvr =
    typename std::remove_cv<typename std::remove_reference<T>::type>::type;

namespace v {

// ------------------------BASIC OPERATIONS-------------

template <class T, std::size_t N, class U, class A, class B> struct BinaryOp {

  const A &a;
  const B &b;

  typedef BinaryOp<T, N - 1, U, A, B> smaller;
  typedef void thisisavvec;
  typedef T value_type;
  static const std::size_t size = N;

  BinaryOp(const A &a, const B &b) : a(a), b(b) {}

  const smaller gsmaller() const { return smaller(a, b); }

  T operator[](std::size_t i) const { return U::apply(a[i], b[i]); }

  void evaluateInto(T *data) const {
    gsmaller().evaluateInto(data);
    data[N - 1] = operator[](N - 1);
  }
};

template <class T, class U, class A, class B> struct BinaryOp<T, 0, U, A, B> {

  BinaryOp(const A &, const B &) {}

  void evaluateInto(T *) const {}
};

template <class T, std::size_t N, class U, class A> struct ScalarOp {

  const A &a;
  T s;

  typedef ScalarOp<T, N - 1, U, A> smaller;
  typedef void thisisavvec;
  typedef T value_type;
  static const std::size_t size = N;

  ScalarOp(const A &a, T s) : a(a), s(s) {}

  const smaller gsmaller() const { return smaller(a, s); }

  T operator[](std::size_t i) const { return U::apply(a[i], s); }

  void evaluateInto(T *data) const {
    gsmaller().evaluateInto(data);
    data[N - 1] = operator[](N - 1);
  }
};

template <class T, class U, class A> struct ScalarOp<T, 0, U, A> {

  ScalarOp(const A &, T) {}

  void evaluateInto(T *) const {}
};

template <class T, std::size_t N, class U, class A> struct ReverseScalarOp {

  const A &a;
  T s;

  typedef ReverseScalarOp<T, N - 1, U, A> smaller;
  typedef void thisisavvec;
  typedef T value_type;
  static const std::size_t size = N;

  ReverseScalarOp(const A &a, T s) : a(a), s(s) {}

  const smaller gsmaller() const { return smaller(a, s); }

  T operator[](std::size_t i) const { return U::apply(s, a[i]); }

  void evaluateInto(T *data) const {
    gsmaller().evaluateInto(data);
    data[N - 1] = operator[](N - 1);
  }
};

template <class T, class U, class A> struct ReverseScalarOp<T, 0, U, A> {

  ReverseScalarOp(const A &, T) {}

  void evaluateInto(T *) const {}
};

template <class T, std::size_t N, class U1, class U2, class A> struct Combine {

  const A &a;

  typedef Combine<T, N - 1, U1, U2, A> smaller;

  Combine(const A &a) : a(a) {}

  const smaller gsmaller() const { return smaller(a); }

  T evaluate() const { return U1::apply(gsmaller().evaluate(), a[N - 1]); }
};

template <class T, class U1, class U2, class A>
struct Combine<T, 1, U1, U2, A> {

  const A &a;

  Combine(const A &a) : a(a) {}

  T evaluate() const { return U2::apply(a[0]); }
};

template <std::size_t N, class A, class B> struct IsEqual {
  bool operator()(const A &a, const B &b) const noexcept {
    return a[N - 1] == b[N - 1] && IsEqual<N - 1, A, B>{}(a, b);
  }
};

template <class A, class B> struct IsEqual<1, A, B> {
  bool operator()(const A &a, const B &b) const noexcept {
    return a[0] == b[0];
  }
};

template <class A, class B, class = typename rem_cvr<A>::thisisavvec,
          class = typename rem_cvr<B>::thisisavvec,
          class = typename std::enable_if<
              std::is_same<typename rem_cvr<A>::value_type,
                           typename rem_cvr<B>::value_type>::value &&
              rem_cvr<A>::size == rem_cvr<B>::size>::type>
struct EqualFunctor {
  bool operator()(const A &a, const B &b) const noexcept {
    return v::IsEqual<rem_cvr<A>::size, A, B>{}(a, b);
  }
};

template <class T, std::size_t N, class U, class A> struct UnaryOp {

  const A &a;

  typedef UnaryOp<T, N - 1, U, A> smaller;
  typedef void thisisavvec;
  typedef T value_type;
  static const std::size_t size = N;

  UnaryOp(const A &a) : a(a) {}

  const smaller gsmaller() { return smaller(a); }

  T operator[](std::size_t i) const { return U::apply(a[i]); }

  void evaluateInto(T *data) const {
    gsmaller().evaluateInto(data);
    data[N - 1] = operator[](N - 1);
  }
};

template <class T, class U, class A> struct UnaryOp<T, 0, U, A> {

  UnaryOp(const A &) {}

  void evaluateInto(T *) {}
};

// ------------------------TYPEDEFS------------------

template <class T> struct IdU {
  static T apply(T a) { return a; }
};

template <class T> struct AddU {
  static T apply(T a, T b) { return a + b; }
};
template <class T, std::size_t N, class A, class B>
using Add = BinaryOp<T, N, AddU<T>, A, B>;
template <class T, std::size_t N, class A>
using SAdd = ScalarOp<T, N, AddU<T>, A>;
template <class T, std::size_t N, class A>
using Sum = Combine<T, N, AddU<T>, IdU<T>, A>;

template <class T> struct SubU {
  static T apply(T a, T b) { return a - b; }
};
template <class T, std::size_t N, class A, class B>
using Sub = BinaryOp<T, N, SubU<T>, A, B>;
template <class T, std::size_t N, class A>
using SSub = ScalarOp<T, N, SubU<T>, A>;
template <class T, std::size_t N, class A>
using ReverseSSub = ReverseScalarOp<T, N, SubU<T>, A>;

template <class T> struct MultU {
  static T apply(T a, T b) { return a * b; }
};
template <class T, std::size_t N, class A, class B>
using Mult = BinaryOp<T, N, MultU<T>, A, B>;
template <class T, std::size_t N, class A>
using SMult = ScalarOp<T, N, MultU<T>, A>;
template <class T, std::size_t N, class A>
using Product = Combine<T, N, MultU<T>, IdU<T>, A>;

template <class T> struct DivU {
  static T apply(T a, T b) { return a / b; }
};
template <class T, std::size_t N, class A, class B>
using Div = BinaryOp<T, N, DivU<T>, A, B>;
template <class T, std::size_t N, class A>
using SDiv = ScalarOp<T, N, DivU<T>, A>;
template <class T, std::size_t N, class A>
using ReverseSDiv = ReverseScalarOp<T, N, DivU<T>, A>;

template <class T> struct NegU {
  static T apply(T a) { return -a; }
};
template <class T, std::size_t N, class A>
using UnaryNeg = UnaryOp<T, N, NegU<T>, A>;
template <class T> struct AbsU {
  static T apply(T a) { return a > 0 ? a : -a; }
};
template <class T, std::size_t N, class A>
using Abs = UnaryOp<T, N, AbsU<T>, A>;

template <class T> struct MinU {
  static T apply(T a, T b) { return a > b ? b : a; }
};
template <class T, std::size_t N, class A, class B>
using ElementwiseMin = BinaryOp<T, N, MinU<T>, A, B>;
template <class T, std::size_t N, class A>
using Min = Combine<T, N, MinU<T>, IdU<T>, A>;

template <class T> struct MaxU {
  static T apply(T a, T b) { return a > b ? a : b; }
};
template <class T, std::size_t N, class A, class B>
using ElementwiseMax = BinaryOp<T, N, MaxU<T>, A, B>;
template <class T, std::size_t N, class A>
using Max = Combine<T, N, MaxU<T>, IdU<T>, A>;

// -----------------------FUNCTIONS----------------------

template <class A, class = typename rem_cvr<A>::thisisavvec>
typename A::value_type sum(const A &a) {
  return Sum<typename A::value_type, A::size, A>(a).evaluate();
}

template <class A, class = typename rem_cvr<A>::thisisavvec>
typename A::value_type product(const A &a) {
  return Product<typename A::value_type, A::size, A>(a).evaluate();
}

template <class A, class = typename rem_cvr<A>::thisisavvec>
typename A::value_type min(const A &a) {
  return Min<typename A::value_type, A::size, A>(a).evaluate();
}

template <class A, class = typename rem_cvr<A>::thisisavvec>
typename A::value_type max(const A &a) {
  return Max<typename A::value_type, A::size, A>(a).evaluate();
}

template <class A, class B, class = typename rem_cvr<A>::thisisavvec,
          class = typename rem_cvr<B>::thisisavvec,
          class = typename std::enable_if<
              std::is_same<typename rem_cvr<A>::value_type,
                           typename rem_cvr<B>::value_type>::value &&
              rem_cvr<A>::size == rem_cvr<B>::size>::type>
ElementwiseMin<typename A::value_type, A::size, A, B>
elementwiseMin(const A &a, const B &b) {
  return {a, b};
}

template <class A, class B, class = typename rem_cvr<A>::thisisavvec,
          class = typename rem_cvr<B>::thisisavvec,
          class = typename std::enable_if<
              std::is_same<typename rem_cvr<A>::value_type,
                           typename rem_cvr<B>::value_type>::value &&
              rem_cvr<A>::size == rem_cvr<B>::size>::type>
ElementwiseMax<typename A::value_type, A::size, A, B>
elementwiseMax(const A &a, const B &b) {
  return {a, b};
}

template <class A, class B, class = typename rem_cvr<A>::thisisavvec,
          class = typename rem_cvr<B>::thisisavvec,
          class = typename std::enable_if<
              std::is_same<typename rem_cvr<A>::value_type,
                           typename rem_cvr<B>::value_type>::value &&
              rem_cvr<A>::size == rem_cvr<B>::size>::type>
typename A::value_type dot(const A &a, const B &b) {
  return sum(a * b);
}

template <class A, class = typename rem_cvr<A>::thisisavvec>
typename A::value_type norm2(const A &a) {
  return sum(a * a);
}

template <class A, class B, class = typename rem_cvr<A>::thisisavvec,
          class = typename rem_cvr<B>::thisisavvec,
          class = typename std::enable_if<
              std::is_same<typename rem_cvr<A>::value_type,
                           typename rem_cvr<B>::value_type>::value &&
              rem_cvr<A>::size == rem_cvr<B>::size>::type>
typename A::value_type dist2(const A &a, const B &b) {
  return norm2(a - b);
}

template <class A, class = typename rem_cvr<A>::thisisavvec>
Abs<typename A::value_type, A::size, A> vabs(const A &a) {
  return {a};
}

// ----------------------VEC DEFINITION-----------------

template <class T, std::size_t N> class Vec {

  template <class U, std::size_t M> struct VecAssign {

    void operator()(T *v, const U &other) {
      v[M - 1] = other[M - 1];
      VecAssign<U, M - 1>{}(v, other);
    }
  };

  template <class U> struct VecAssign<U, 1> {

    void operator()(T *v, const U &other) { v[0] = other[0]; }
  };

  void assign(const T *other) { VecAssign<const T *, N>{}(data, other); }

  template <class U, class = typename rem_cvr<U>::thisisavvec,
            class = typename std::enable_if<
                std::is_same<T, typename rem_cvr<U>::value_type>::value &&
                N == rem_cvr<U>::size>::type>
  void assign(const U &other) {
    VecAssign<U, N>{}(data, other);
  }

public:
  typedef Vec<T, N> smaller;
  typedef void thisisavvec;
  typedef T value_type;
  static const std::size_t size = N;

  T data[N];

  Vec() {}

  Vec(const std::initializer_list<T> &inil) {
    if (inil.size() != N) {
      return;
    }
    std::size_t i = 0;
    for (const T &val : inil) {
      data[i++] = val;
    }
  }

  template <class U, class = typename rem_cvr<U>::thisisavvec,
            class = typename std::enable_if<
                std::is_same<T, typename rem_cvr<U>::value_type>::value &&
                N == rem_cvr<U>::size>::type>
  Vec(const U &other) {
    assign<U>(other);
  }

  explicit Vec(const T *other) { assign(other); }

  template <class U, class = typename rem_cvr<U>::thisisavvec,
            class = typename std::enable_if<
                std::is_same<T, typename rem_cvr<U>::value_type>::value &&
                N == rem_cvr<U>::size>::type>
  Vec<T, N> &operator=(const U &other) {
    assign<U>(other);
    return *this;
  }

  T &operator[](std::size_t i) { return data[i]; }
  const T &operator[](std::size_t i) const { return data[i]; }

  void copyFrom(const T *ptr) { assign(ptr); }

  template <class U, class = typename rem_cvr<U>::thisisavvec,
            class = typename std::enable_if<
                std::is_same<T, typename rem_cvr<U>::value_type>::value &&
                N == rem_cvr<U>::size>::type>
  Vec<T, N> &operator+=(const U &other) {
    return operator=(*this + other);
  }
  template <class U, class = typename rem_cvr<U>::thisisavvec,
            class = typename std::enable_if<
                std::is_same<T, typename rem_cvr<U>::value_type>::value &&
                N == rem_cvr<U>::size>::type>
  Vec<T, N> &operator-=(const U &other) {
    return operator=(*this - other);
  }
  template <class U, class = typename rem_cvr<U>::thisisavvec,
            class = typename std::enable_if<
                std::is_same<T, typename rem_cvr<U>::value_type>::value &&
                N == rem_cvr<U>::size>::type>
  Vec<T, N> &operator*=(const U &other) {
    return operator=(*this *other);
  }
  template <class U, class = typename rem_cvr<U>::thisisavvec,
            class = typename std::enable_if<
                std::is_same<T, typename rem_cvr<U>::value_type>::value &&
                N == rem_cvr<U>::size>::type>
  Vec<T, N> &operator/=(const U &other) {
    return operator=(*this / other);
  }

  Vec<T, N> &operator+=(T val) { return operator=(*this + val); }
  Vec<T, N> &operator-=(T val) { return operator=(*this - val); }
  Vec<T, N> &operator*=(T val) { return operator=(*this *val); }
  Vec<T, N> &operator/=(T val) { return operator=(*this / val); }
};

template <class T, T value, std::size_t N> struct ConstantVec {

  typedef ConstantVec<T, value, N - 1> smaller;
  typedef void thisisavvec;
  typedef T value_type;
  static const std::size_t size = N;

  constexpr T operator[](std::size_t) const { return value; }
};

// -----------------------VEC HELPERS-------------------

template <std::size_t N> using IVec = Vec<std::int32_t, N>;
template <std::size_t N> using FVec = Vec<float, N>;
template <std::size_t N> using DVec = Vec<double, N>;

// Source: https://en.wikipedia.org/wiki/MurmurHash
// I modified the ending of the algo so that it can
// run faster.

template <std::size_t N> struct IVecHash {

protected:
  std::uint32_t murmurscram(std::uint32_t val) const noexcept {
    val *= 0xcc9e2d51;
    val = (val << 15) | (val >> 17);
    val *= 0x1b873593;
    return val;
  }

  std::uint32_t murmurstep(std::uint32_t h, std::uint32_t val) const noexcept {
    h ^= murmurscram(val);
    h = (h << 13) | (h >> 19);
    return h * 5 + 0xe6546b64;
  }

public:
  std::uint32_t impl(const std::int32_t *argp) const noexcept {
    return murmurstep(IVecHash<N - 1>().impl(argp), argp[N - 1]);
  }

  std::size_t operator()(const IVec<N> &arg) const noexcept {
    return IVecHash<N - 1>().impl(&arg[0]) ^ murmurscram(arg[N - 1]);
  }
};

template <> struct IVecHash<1> : protected IVecHash<2> {

  std::uint32_t impl(const std::int32_t *argp) const noexcept {
    return murmurstep(0 /* seed */, argp[0]);
  }

  std::size_t operator()(const IVec<1> &arg) const noexcept {
    return murmurscram(arg[0]);
  }
};

} // namespace v

// -----------------------GENERIC OPERATORS-------------

template <class A, class B, class = typename rem_cvr<A>::thisisavvec,
          class = typename rem_cvr<B>::thisisavvec,
          class = typename std::enable_if<
              std::is_same<typename rem_cvr<A>::value_type,
                           typename rem_cvr<B>::value_type>::value &&
              rem_cvr<A>::size == rem_cvr<B>::size>::type>
v::Add<typename A::value_type, A::size, A, B> operator+(const A &a,
                                                        const B &b) {
  return {a, b};
}

template <class A, class B, class = typename rem_cvr<A>::thisisavvec,
          class = typename rem_cvr<B>::thisisavvec,
          class = typename std::enable_if<
              std::is_same<typename rem_cvr<A>::value_type,
                           typename rem_cvr<B>::value_type>::value &&
              rem_cvr<A>::size == rem_cvr<B>::size>::type>
v::Sub<typename A::value_type, A::size, A, B> operator-(const A &a,
                                                        const B &b) {
  return {a, b};
}

template <class A, class B, class = typename rem_cvr<A>::thisisavvec,
          class = typename rem_cvr<B>::thisisavvec,
          class = typename std::enable_if<
              std::is_same<typename rem_cvr<A>::value_type,
                           typename rem_cvr<B>::value_type>::value &&
              rem_cvr<A>::size == rem_cvr<B>::size>::type>
v::Mult<typename A::value_type, A::size, A, B> operator*(const A &a,
                                                         const B &b) {
  return {a, b};
}

template <class A, class B, class = typename rem_cvr<A>::thisisavvec,
          class = typename rem_cvr<B>::thisisavvec,
          class = typename std::enable_if<
              std::is_same<typename rem_cvr<A>::value_type,
                           typename rem_cvr<B>::value_type>::value &&
              rem_cvr<A>::size == rem_cvr<B>::size>::type>
v::Div<typename A::value_type, A::size, A, B> operator/(const A &a,
                                                        const B &b) {
  return {a, b};
}

template <class A, class = typename rem_cvr<A>::thisisavvec>
v::UnaryNeg<typename A::value_type, A::size, A> operator-(const A &a) {
  return {a};
}

template <class A, class = typename rem_cvr<A>::thisisavvec>
v::SAdd<typename A::value_type, A::size, A>
operator+(const A &a, typename A::value_type s) {
  return {a, s};
}

template <class A, class = typename rem_cvr<A>::thisisavvec>
v::SAdd<typename A::value_type, A::size, A> operator+(typename A::value_type s,
                                                      const A &a) {
  return {a, s};
}

template <class A, class = typename rem_cvr<A>::thisisavvec>
v::SSub<typename A::value_type, A::size, A>
operator-(const A &a, typename A::value_type s) {
  return {a, s};
}

template <class A, class = typename rem_cvr<A>::thisisavvec>
v::ReverseSSub<typename A::value_type, A::size, A>
operator-(typename A::value_type s, const A &a) {
  return {a, s};
}

template <class A, class = typename rem_cvr<A>::thisisavvec>
v::SMult<typename A::value_type, A::size, A>
operator*(const A &a, typename A::value_type s) {
  return {a, s};
}

template <class A, class = typename rem_cvr<A>::thisisavvec>
v::SMult<typename A::value_type, A::size, A> operator*(typename A::value_type s,
                                                       const A &a) {
  return {a, s};
}

template <class A, class = typename rem_cvr<A>::thisisavvec>
v::SDiv<typename A::value_type, A::size, A>
operator/(const A &a, typename A::value_type s) {
  return {a, s};
}

template <class A, class = typename rem_cvr<A>::thisisavvec>
v::ReverseSDiv<typename A::value_type, A::size, A>
operator/(typename A::value_type s, const A &a) {
  return {a, s};
}

template <class A, class B, class = typename rem_cvr<A>::thisisavvec,
          class = typename rem_cvr<B>::thisisavvec,
          class = typename std::enable_if<
              std::is_same<typename rem_cvr<A>::value_type,
                           typename rem_cvr<B>::value_type>::value &&
              rem_cvr<A>::size == rem_cvr<B>::size>::type>
bool operator==(const A &a, const B &b) {
  return v::EqualFunctor<A, B>{}(a, b);
}

template <class A, class B, class = typename rem_cvr<A>::thisisavvec,
          class = typename rem_cvr<B>::thisisavvec,
          class = typename std::enable_if<
              std::is_same<typename rem_cvr<A>::value_type,
                           typename rem_cvr<B>::value_type>::value &&
              rem_cvr<A>::size == rem_cvr<B>::size>::type>
bool operator!=(const A &a, const B &b) {
  return !(a == b);
}

template <class A, class = typename rem_cvr<A>::thisisavvec>
std::ostream &operator<<(std::ostream &os, const A &a) {
  os << "v[";
  for (std::size_t i = 0; i < A::size - 1; i++) {
    os << a[i] << " ";
  }
  os << a[A::size - 1] << "] ";
  return os;
}

#endif // V_VECTOR_HPP_
