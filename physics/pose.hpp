#ifndef ROLLER_PHYSICS_POSE_HPP_
#define ROLLER_PHYSICS_POSE_HPP_

#include <cmath>

#include "util.hpp"

namespace roller {

v::DVec<4> quaternionConj(const v::DVec<4> &a) {
  return {a[0], -a[1], -a[2], -a[3]};
}
v::DVec<4> quaternionMult(const v::DVec<4> &a, const v::DVec<4> &b) {
  return {a[0] * b[0] - a[1] * b[1] - a[2] * b[2] - a[3] * b[3],
          a[0] * b[1] + a[1] * b[0] + a[2] * b[3] - a[3] * b[2],
          a[0] * b[2] - a[1] * b[3] + a[2] * b[0] + a[3] * b[1],
          a[0] * b[3] + a[1] * b[2] - a[2] * b[1] + a[3] * b[0]};
}
v::DVec<3> applyQuaternion(const v::DVec<3> &p, const v::DVec<4> &q) {
  v::DVec<4> ch{0, p[0], p[1], p[2]};
  v::DVec<4> ro = quaternionMult(q, quaternionMult(ch, quaternionConj(q)));
  return {ro[1], ro[2], ro[3]};
}
// Source: https://stackoverflow.com/a/12934750
v::DVec<4> normalizeQuaternion(const v::DVec<4> &q) {
  double sqmag = v::norm2(q);
  if ((1 - 2.10735e-8) < sqmag && sqmag < (1 + 2.10735e-8)) {
    return (2 / (1 + sqmag)) * q;
  }
  return fastInvSqrt(sqmag) * q;
}
// yields a quaternion that represents constant angular rotation by w * 2
v::DVec<4> getRotQuaternion(const v::DVec<3> &w) {
  double wnorm2 = v::norm2(w);
  if (wnorm2 < 1e-12) {
    return {1, 0, 0, 0};
  }
  // XXX: FIND WHAT THIS NUMBER SHOULD ACTUALLY BE
  if (wnorm2 < 0.05) {
    // first-order Taylor approximation
    return {1, w[0], w[1], w[2]};
  }
  // XXX: USE FAST ALGORITHMS FOR THESE SQRT AND SINC
  double wnorm = std::sqrt(wnorm2);
  v::DVec<3> ww = (std::sin(wnorm) / wnorm) * w;
  return {std::cos(wnorm), ww[0], ww[1], ww[2]};
}
// applies getRotQuaternion(w) to q
v::DVec<4> applyDAngVelo(const v::DVec<4> &q, const v::DVec<3> &w) {
  return normalizeQuaternion(quaternionMult(getRotQuaternion(w), q));
}

struct Pose {

  v::DVec<3> p{0, 0, 0};    /// translation
  v::DVec<4> q{1, 0, 0, 0}; /// rotation

  DMat3x3 toRotationMatrix() const {
    double x2 = 2 * q[1] * q[1];
    double y2 = 2 * q[2] * q[2];
    double z2 = 2 * q[3] * q[3];
    double xy = q[1] * q[2];
    double xz = q[1] * q[3];
    double yz = q[2] * q[3];
    double sx = q[0] * q[1];
    double sy = q[0] * q[2];
    double sz = q[0] * q[3];
    return {{1 - y2 - z2, 2 * (xy - sz), 2 * (xz + sy)},
            {2 * (xy + sz), 1 - x2 - z2, 2 * (yz - sx)},
            {2 * (xz - sy), 2 * (yz + sx), 1 - x2 - y2}};
  }
  v::DVec<16> toProjMat() const {
    DMat3x3 r = toRotationMatrix();
    return {r.a[0], r.a[1], r.a[2], p[0], r.b[0], r.b[1], r.b[2], p[1],
            r.c[0], r.c[1], r.c[2], p[2], 0,      0,      0,      1};
  }
  v::DVec<16> toProjMat(const SliceDirs &sd) const {
    v::DVec<16> s = genProjMat(sd);
    v::DVec<16> t = toProjMat();
    return matMult4x4(s, t);
  }

  v::DVec<3> toShiftWorldCoords(const v::DVec<3> &c) const {
    return applyQuaternion(c, q);
  }
  v::DVec<3> toWorldCoords(const v::DVec<3> &c) const {
    return p + toShiftWorldCoords(c);
  }
  v::DVec<3> fromShiftWorldCoords(const v::DVec<3> &sc) const {
    v::DVec<4> sh{0, sc[0], sc[1], sc[2]};
    v::DVec<4> rh = quaternionMult(quaternionConj(q), quaternionMult(sh, q));
    return {rh[1], rh[2], rh[3]};
  }
  v::DVec<3> fromWorldCoords(const v::DVec<3> &c) const {
    return fromShiftWorldCoords(c - p);
  }
};

} // namespace roller

#endif // ROLLER_PHYSICS_POSE_HPP_

