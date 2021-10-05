#ifndef ROLLER_PHYSICS_PHYS_INFO_HPP_
#define ROLLER_PHYSICS_PHYS_INFO_HPP_

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
  v::DVec<3> lm, am; /// linear momentum, angular momentum

  double massi;  /// 1 / mass
  DMat3x3 ineri; /// I^{-1}

  AuxPhysInfo getAuxInfo() const {
    AuxPhysInfo ret;
    ret.velo = massi * lm;
    ret.rot = pose.toRotationMatrix();
    ret.riner = ret.rot * ineri * ret.rot.transpose();
    ret.omega = ret.riner * am;
    return ret;
  }
};

struct AABB {
  v::DVec<3> a, b; /// min, max

  bool intersects(const AABB &c) const {
    return (a[0] < c.b[0]) && (b[0] > c.a[0]) && (a[1] < c.b[1]) &&
           (b[1] > c.a[1]) && (a[2] < c.b[2]) && (b[2] > c.a[2]);
  }
};

struct OBB {
  v::DVec<3> b, x, y, z; /// basically b+gx+hy+iz where 0<=g,h,i<=1
  v::DVec<3> a, c;

  OBB(const v::DVec<3> &b, const v::DVec<3> &x, const v::DVec<3> &y,
      const v::DVec<3> &z)
      : b(b), x(x), y(y), z(z), a{v::dot(b, x), v::dot(b, y), v::dot(b, z)},
        c(a + 1) {}

  AABB getAABB() const {
    // could use extrema below, but unnecessary multiplications
    AABB ret{b, b};
    if (x[0] > 0) {
      ret.b[0] += x[0];
    } else {
      ret.a[0] += x[0];
    }
    if (y[0] > 0) {
      ret.b[0] += y[0];
    } else {
      ret.a[0] += y[0];
    }
    if (z[0] > 0) {
      ret.b[0] += z[0];
    } else {
      ret.a[0] += z[0];
    }
    return ret;
  }

  /// minimum and maximum of u dot v, where v is a point in the OBB
  v::DVec<2> extrema(const v::DVec<3> &u) const {
    double s = v::dot(u, b);
    double a = 0, b = 0;
    double t;
    if ((t = v::dot(u, x)) > 0) {
      b += t;
    } else {
      a += t;
    }
    if ((t = v::dot(u, y)) > 0) {
      b += t;
    } else {
      a += t;
    }
    if ((t = v::dot(u, z)) > 0) {
      b += t;
    } else {
      a += t;
    }
    return {a, b};
  }

  bool intersects0(const OBB &o) const {
    v::DVec<2> tx = o.extrema(x);
    if (tx[0] > c[0] || tx[1] < a[0]) {
      return false;
    }
    v::DVec<2> ty = o.extrema(y);
    if (ty[0] > c[1] || ty[1] < a[1]) {
      return false;
    }
    v::DVec<2> tz = o.extrema(z);
    if (tz[0] > c[2] || tz[1] < a[2]) {
      return false;
    }
    return true;
  }
  bool intersects(const OBB &o) const {
    return intersects0(o) || o.intersects0(*this);
  }
};

} // namespace roller

#endif // ROLLER_PHYSICS_PHYS_INFO_HPP_

