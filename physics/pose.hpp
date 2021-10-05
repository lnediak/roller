#ifndef ROLLER_PHYSICS_POSE_HPP_
#define ROLLER_PHYSICS_POSE_HPP_

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

struct Pose {

  v::DVec<3> p; /// translation
  v::DVec<4> q; /// rotation

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

  v::DVec<3> toWorldCoords(const v::DVec<3> &c) const {
    v::DVec<4> ch{1, c[0], c[1], c[2]};
    v::DVec<4> ro = quaternionMult(q, quaternionMult(ch, quaternionConj(q)));
    return p + v::DVec<3>{ro[1], ro[2], ro[3]};
  }
  v::DVec<3> fromWorldCoords(const v::DVec<3> &c) const {
    v::DVec<3> sc = c - p;
    v::DVec<4> sh{1, sc[0], sc[1], sc[2]};
    v::DVec<4> rh = quaternionMult(quaternionConj(q), quaternionMult(sh, q));
    return {rh[1], rh[2], rh[3]};
  }
};

} // namespace roller

#endif // ROLLER_PHYSICS_POSE_HPP_

