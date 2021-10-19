#ifndef ROLLER_PHYSICS_OBB_HPP_
#define ROLLER_PHYSICS_OBB_HPP_

#include "aabb.hpp"
#include "vector.hpp"

namespace roller {

struct OBB {
  v::DVec<3> b, s, x, y, z; /// basically b+gx+hy+iz where 0<=g,h,i<=s
  v::DVec<3> a, c;

  OBB() {}
  OBB(const v::DVec<3> &b, const v::DVec<3> &s, const v::DVec<3> &x,
      const v::DVec<3> &y, const v::DVec<3> &z)
      : b(b), s(s), x(x), y(y),
        z(z), a{v::dot(b, x), v::dot(b, y), v::dot(b, z)}, c(a + s) {}

  bool isIn(const v::DVec<3> &p) const {
    v::DVec<3> pd = {v::dot(p, x), v::dot(p, y), v::dot(p, z)};
    return a[0] <= pd[0] && a[1] <= pd[1] && a[2] <= pd[2] && pd[0] <= c[0] &&
           pd[1] <= c[1] && pd[2] <= c[2];
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
      if (-1e-3 < dd[ind] && dd[ind] < 1e-3) {
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
};

} // namespace roller

#endif // ROLLER_PHYSICS_OBB_HPP_

