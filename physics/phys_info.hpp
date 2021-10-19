#ifndef ROLLER_PHYSICS_PHYS_INFO_HPP_
#define ROLLER_PHYSICS_PHYS_INFO_HPP_

#include <algorithm>
#include <cmath>
#include <unordered_set>

#include "pose.hpp"

namespace roller {

struct AuxPhysInfo {
  v::DVec<3> velo; /// velocity

  DMat3x3 rot;      /// rotation matrix
  DMat3x3 riner;    /// I^{-1}(t)
  v::DVec<3> omega; /// angular velocity
};

struct PhysInfo {
  Pose pose;
  v::DVec<3> lm{0, 0, 0}, am{0, 0, 0}; /// linear momentum, angular momentum

  double mass = 0;
  double massi = 0; /// 1 / mass
  DMat3x3 iner{{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
  DMat3x3 ineri{{0, 0, 0}, {0, 0, 0}, {0, 0, 0}}; /// iner^{-1}

  PhysInfo() {}
  PhysInfo(double massi, const DMat3x3 &ineri)
      : mass(massi < 1e-12 ? 0 : 1 / massi), massi(massi),
        iner(massi < 1e-12 ? DMat3x3{{0, 0, 0}, {0, 0, 0}, {0, 0, 0}}
                           : ineri.inverse()),
        ineri(ineri) {}

  void updateAuxRotInfo(AuxPhysInfo &aux) const {
    aux.rot = pose.toRotationMatrix();
    if (mass) {
      aux.riner = aux.rot * ineri * aux.rot.transpose();
      aux.omega = aux.riner * am;
    } else {
      aux.riner = ineri;
      aux.omega = {0, 0, 0};
    }
  }
  AuxPhysInfo getAuxInfo() const {
    AuxPhysInfo ret;
    if (mass) {
      ret.velo = massi * lm;
    } else {
      ret.velo = {0, 0, 0};
    }
    updateAuxRotInfo(ret);
    return ret;
  }

  /// velocity at point p
  v::DVec<3> getVelocity(const AuxPhysInfo &aux, const v::DVec<3> &p) const {
    return aux.velo + cross3(aux.omega, p - pose.p);
  }

  /// Uses Euler's method. Does not update velocity.
  template <bool updateAux>
  void stepTime(AuxPhysInfo &aux, double dt, double g) {
    if (!mass) {
      return;
    }
    pose.p += dt * aux.velo;
    pose.q = applyDAngVelo(pose.q, 0.5 * dt * aux.omega);
    if (updateAux) {
      updateAuxRotInfo(aux);
    }
  }
  template <bool updateAux>
  void stepVelo(AuxPhysInfo &aux, double dt, double g) {
    double imp = dt * g;
    lm[2] -= imp * mass;
    if (updateAux) {
      aux.velo[2] -= imp;
    }
  }
};

struct AABB {
  v::DVec<3> m[2]; /// min, max

  bool intersects(const AABB &c) const {
    return (m[0][0] < c.m[1][0]) && (m[1][0] > c.m[0][0]) &&
           (m[0][1] < c.m[1][1]) && (m[1][1] > c.m[0][1]) &&
           (m[0][2] < c.m[1][2]) && (m[1][2] > c.m[0][2]);
  }
};

struct Contact {
  double dist;
  v::DVec<3> p, n;

  Contact lo(const Contact &c) const {
    if (dist >= 0) {
      return dist < c.dist ? *this : c;
    }
    if (c.dist >= 0) {
      return *this;
    }
    return dist < c.dist ? c : *this;
  }
  Contact hi(const Contact &c) const { return dist < c.dist ? *this : c; }
};

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

  /// XXX: TRIGGER WARNING: NOT OPTIMIZED
  AABB getAABB() const {
    v::DVec<3> c0 = b;
    v::DVec<3> c1 = b + s[0] * x;
    v::DVec<3> c2 = b + s[1] * y;
    v::DVec<3> c3 = b + s[2] * z;
    v::DVec<3> c4 = c1 + s[1] * y;
    v::DVec<3> c5 = c2 + s[2] * z;
    v::DVec<3> c6 = c3 + s[0] * x;
    v::DVec<3> c7 = c4 + s[2] * z;
    return {
        elementwiseMin(
            elementwiseMin(elementwiseMin(c0, c1), elementwiseMin(c2, c3)),
            elementwiseMin(elementwiseMin(c4, c5), elementwiseMin(c6, c7))),
        elementwiseMax(
            elementwiseMax(elementwiseMax(c0, c1), elementwiseMax(c2, c3)),
            elementwiseMax(elementwiseMax(c4, c5), elementwiseMax(c6, c7)))};
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
  // finds the value of v such that u dot v is maximum (a.k.a. support)
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
  /// returns whether *this intersects the specified line
  /// segment, and outputs a point on both into `out`
  bool getAlongSeg(const v::DVec<3> &p, const v::DVec<3> &q,
                   v::DVec<3> &out) const {
    v::DVec<3> d = q - p;
    v::DVec<3> dd = {v::dot(d, x), v::dot(d, y), v::dot(d, z)};
    v::DVec<3> pp = {v::dot(p, x), v::dot(p, y), v::dot(p, z)};
    double lo = 0, hi = 1;
    for (int ind = 0; ind < 3; ind++) {
      if (-1e-3 < dd[ind] && dd[ind] < 1e-3) {
        if ((a[ind] <= pp[ind]) ^ (pp[ind] <= c[ind])) {
          return false;
        }
        continue;
      }
      double tlo = (a[ind] - pp[ind]) / dd[ind];
      double thi = (c[ind] - pp[ind]) / dd[ind];
      if (tlo > thi) {
        lo = std::fmax(lo, thi);
        hi = std::fmin(hi, tlo);
      } else {
        lo = std::fmax(lo, tlo);
        hi = std::fmin(hi, thi);
      }
    }
    if (lo > hi) {
      return false;
    }
    out = p + 0.5 * (lo + hi) * d;
    return true;
  }
  v::DVec<3> dirD(int i) const {
    switch (i) {
    case 0:
      return x;
    case 1:
      return y;
    case 2:
      return z;
    }
    return {0, 0, 0};
  }
  /// tries to find a point on o's wireframe in *this
  bool getInInter(const OBB &o, v::DVec<3> &out) const {
    for (int d1 = 0; d1 < 3; d1++) {
      int d2 = d1 >= 2 ? 0 : d1 + 1;
      int d3 = d2 >= 2 ? 0 : d2 + 1;
      for (int m1 = 0; m1 < 2; m1++) {
        for (int m2 = 0; m2 < 2; m2++) {
          v::DVec<3> pp =
              o.b + m1 * o.s[d2] * o.dirD(d2) + m2 * o.s[d3] * o.dirD(d3);
          if (getAlongSeg(pp, pp + o.s[d1] * o.dirD(d1), out)) {
            return true;
          }
        }
      }
    }
    return false;
  }
};

/// XXX: TRIGGER WARNING: NOT OPTIMIZED
struct OBBIntersector {

  OBB p, q;

  OBBIntersector(const OBB &p, const OBB &q) { initInt(p, q); }

  /// XXX: evaluate orientation-related constants for optimization
  void initInt(const OBB &pp, const OBB &qq) {
    p = pp;
    q = qq;
  }

  v::DVec<3> ee(const v::DVec<3> &a, const v::DVec<3> &b) {
    v::DVec<3> cros = cross3(a, b);
    double nrm = v::norm2(cros);
    if (nrm < 1e-12) {
      return p.x;
    }
    return cros / std::sqrt(nrm);
  }
  /// assuming u is unit vector
  bool vc(const v::DVec<3> &u) {
    v::DVec<2> a = p.extrema(u);
    v::DVec<2> b = q.extrema(u);
    return a[0] > b[1] || a[1] < b[0];
  }
  Contact getInts() {
    v::DVec<3> axes[] = {p.x,          p.y,          p.z, // p vert-face
                         q.x,          q.y,          q.z, // q vert-face
                         ee(p.x, q.x), ee(p.x, q.y), ee(p.x, q.z),
                         ee(p.y, q.x), ee(p.y, q.y), ee(p.y, q.z),
                         ee(p.z, q.x), ee(p.z, q.y), ee(p.z, q.z)};
    Contact ret;
    ret.dist = 100;
    for (const v::DVec<3> &axis : axes) {
      v::DVec<2> a = p.extrema(axis);
      v::DVec<2> b = q.extrema(axis);
      double forD = a[0] - b[1];
      double bakD = b[0] - a[1];
      if (forD > 0 || bakD > 0) {
        ret.dist = std::fmax(forD, bakD);
        return ret;
      }
      double tmpd;
      v::DVec<3> tmpn;
      if (forD > bakD) {
        tmpd = forD;
        tmpn = axis;
      } else {
        tmpd = bakD;
        tmpn = -axis;
      }
      if (ret.dist > 0 || ret.dist < tmpd) {
        ret.dist = tmpd;
        ret.n = tmpn;
      }
    }
    if (!p.getInInter(q, ret.p)) {
      if (!q.getInInter(p, ret.p)) {
        std::cout << "I'm gonna BAIL!!!!!!!" << std::endl;
        ret.dist = 100;
      }
    }

    std::cout << "OBBIntersector.getInts: " << ret.dist << " " << ret.p << " "
              << ret.n << std::endl;
    std::cout << "with " << p.b << p.s << p.x << p.y << p.z << std::endl;
    std::cout << "and " << q.b << q.s << q.x << q.y << q.z << std::endl;
    return ret;
  }
};

} // namespace roller

#endif // ROLLER_PHYSICS_PHYS_INFO_HPP_

