#ifndef ROLLER_PHYSICS_OBB_HPP_
#define ROLLER_PHYSICS_OBB_HPP_

#include <cmath>

#include "aabb.hpp"
#include "screw.hpp"
#include "util.hpp"

namespace roller {

AABB getArcAABB(const v::DVec<3> &p, const ScrewM &sm, double tol) {
  double rnorm = std::sqrt(v::norm2(p - sm.center));
  double alpha = std::acos(rnorm / (tol + rnorm));
  double smrad = std::sqrt(v::norm2(sm.omega));
  double alphat1 = alpha / smrad;
  double alphat2 = std::tan(alpha) / alpha;
  v::Dvec<3> nomega = sm.omega / smrad;
  double theta = 0;
  v::DVec<3> currp = p;
  AABB toret = {p, p};
  while (theta + 2 * alpha < smrad) {
    v::DVec<3> rvec = currp - sm.center - sm.velo * (theta / smrad);
    v::DVec<3> outp =
        currp + sm.velo * alphat1 + v::cross3(sm.omega, rvec) * alphat2;
    toret = toret.combine({outp, outp});
    theta += 2 * alpha;
    currp = p;
    apply2Vec(currp, mult(theta / smrad, sm));
  }
  alpha = (smrad - theta) / 2;
  alphat1 = alpha / smrad;
  alphat2 = std::tan(alpha) / alpha;
  v::DVec<3> rvec = currp - sm.center - sm.velo * (theta / smrad);
  v::DVec<3> outp =
      currp + sm.velo * alphat1 + v::cross3(sm.omega, rvec) * alphat2;
  toret = toret.combine({outp, outp});
  currp = p;
  apply2Vec(currp, sm);
  return toret.combine({currp, currp});
}

struct OBB {
  v::DVec<3> b, s, x, y, z; /// basically b+gx+hy+iz where 0<=g,h,i<=s
  v::DVec<3> a, c;

  double tol = 1e-2;

  OBB() {}
  OBB(const v::DVec<3> &b, const v::DVec<3> &s, const v::DVec<3> &x,
      const v::DVec<3> &y, const v::DVec<3> &z)
      : b(b), s(s), x(x), y(y),
        z(z), a{v::dot(b, x), v::dot(b, y), v::dot(b, z)}, c(a + s) {}

  void properify() {
    a = {v::dot(b, x), v::dot(b, y), v::dot(b, z)};
    c = a + s;
  }

  bool isIn(const v::DVec<3> &p) const {
    v::DVec<3> pd = {v::dot(p, x), v::dot(p, y), v::dot(p, z)};
    return a[0] <= pd[0] && a[1] <= pd[1] && a[2] <= pd[2] && pd[0] <= c[0] &&
           pd[1] <= c[1] && pd[2] <= c[2];
  }

  v::DVec<3> center() const {
    return b + (s[0] / 2) * x + (s[1] / 2) * y + (s[2] / 2) * z;
  }
  /// minimum and maximum of u dot v, where v is a point in the OBB
  v::DVec<2> extrema(const v::DVec<3> &u) const {
    double o = v::dot(u, b);
    double lo = o, hi = o;
    double t;
    if ((t = v::dot(u, x)) > 0) {
      hi += s[0] * t;
    } else {
      lo += s[0] * t;
    }
    if ((t = v::dot(u, y)) > 0) {
      hi += s[1] * t;
    } else {
      lo += s[1] * t;
    }
    if ((t = v::dot(u, z)) > 0) {
      hi += s[2] * t;
    } else {
      lo += s[2] * t;
    }
    return {lo, hi};
  }
  /// finds the value of v such that u dot v is maximum (a.k.a. support)
  v::DVec<3> maximize(const v::DVec<3> &u) const {
    v::DVec<3> ret = b;
    if (v::dot(u, x) > 0) {
      ret += s[0] * x;
    }
    if (v::dot(u, y) > 0) {
      ret += s[1] * y;
    }
    if (v::dot(u, z) > 0) {
      ret += s[2] * z;
    }
    return ret;
  }

  /// for getAABB
  v::DVec<2> extremaAxis(int dim) const {
    double o = b[dim];
    double lo = o, hi = o;
    if (x[dim] > 0) {
      hi += s[0] * x[dim];
    } else {
      lo += s[0] * x[dim];
    }
    if (y[dim] > 0) {
      hi += s[1] * y[dim];
    } else {
      lo += s[1] * y[dim];
    }
    if (z[dim] > 0) {
      hi += s[2] * z[dim];
    } else {
      lo += s[2] * z[dim];
    }
    return {lo, hi};
  }
  AABB getAABB() const {
    v::DVec<2> xe = extremaAxis(0);
    v::DVec<2> ye = extremaAxis(1);
    v::DVec<2> ze = extremaAxis(2);
    return {{{xe[0], ye[0], ze[0]}, {xe[1], ye[1], ze[1]}}};
  }
  AABB getAABB(const ScrewM &sm) const {
    // XXX: proper limits maybe?
    if (v::norm2(sm.velo) < tol / 2 && v::norm2(sm.omega) < tol / 100) {
      return getAABB();
    }
    return getArcAABB(b, sm, tol)
        .combine(
            getArcAABB(b + s[0] * x, sm, tol)
                .combine(
                    getArcAABB(b + s[1] * y, sm, tol)
                        .combine(
                            getArcAABB(b + s[0] * x + s[1] * y, sm, tol)
                                .combine(
                                    getArcAABB(b + s[2] * z, sm, tol)
                                        .combine(
                                            getArcAABB(b + s[0] * x + s[2] * z,
                                                       sm, tol)
                                                .combine(
                                                    getArcAABB(b + s[1] * y +
                                                                   s[2] * z,
                                                               sm, tol)
                                                        .combine(getArcAABB(
                                                            b + s[0] * x +
                                                                s[1] * y +
                                                                s[2] * z,
                                                            sm, tol))))))));
  }

  /// constructs an obb with same x,y,z as *this but contains o
  OBB wrapOBB(const OBB &o) const {
    v::DVec<2> xe = o.extrema(x);
    v::DVec<2> ye = o.extrema(y);
    v::DVec<2> ze = o.extrema(z);

    OBB ret;
    ret.b = xe[0] * x + ye[0] * y + ze[0] * z;
    ret.x = x;
    ret.y = y;
    ret.z = z;
    ret.a = {xe[0], ye[0], ze[0]};
    ret.c = {xe[1], ye[1], ze[1]};
    ret.s = ret.c - ret.a;
    return ret;
  }
  /// assuming the same orientation, expands *this to contain o
  void fattenOBB(const OBB &o) {
    a = v::elementwiseMin(a, o.a);
    c = v::elementwiseMax(c, o.c);
    s = c - a;
    b = a[0] * x + a[1] * y + a[2] * z;
  }
  void inflate(double d) {
    a -= d;
    c += d;
    s = c - a;
    b = a[0] * x + a[1] * y + a[2] * z;
  }

  /// returns the lower corner of the edge
  v::DVec<3> maxEdge(const v::DVec<3> &u, int i) const {
    v::DVec<3> ret = b;
    if (i != 0 && v::dot(u, x) > 0) {
      ret += s[0] * x;
    }
    if (i != 1 && v::dot(u, y) > 0) {
      ret += s[1] * y;
    }
    if (i != 2 && v::dot(u, z) > 0) {
      ret += s[2] * z;
    }
    return ret;
  }

  /// returns a point on the line segment close to or in *this
  v::DVec<3> getAlongSeg(const v::DVec<3> &p, const v::DVec<3> &q) const {
    v::DVec<3> d = q - p;
    v::DVec<3> dd = {v::dot(d, x), v::dot(d, y), v::dot(d, z)};
    v::DVec<3> pp = {v::dot(p, x), v::dot(p, y), v::dot(p, z)};
    double lo = 0, hi = 1;
    for (int ind = 0; ind < 3; ind++) {
      if (-1e-6 < dd[ind] && dd[ind] < 1e-6) {
        // if ((a[ind] <= pp[ind]) ^ (pp[ind] <= c[ind])) {
        //   return false;
        // }
        continue;
      }
      double tlo = (a[ind] - pp[ind]) / dd[ind];
      double thi = (c[ind] - pp[ind]) / dd[ind];
      if (tlo > thi) {
        double tmp = tlo;
        tlo = thi;
        thi = tmp;
      }
      if (lo < tlo) {
        lo = tlo;
      }
      if (hi > thi) {
        hi = thi;
      }
    }
    return p + 0.5 * (lo + hi) * d;
  }
};

} // namespace roller

#endif // ROLLER_PHYSICS_OBB_HPP_
