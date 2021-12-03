#ifndef ROLLER_PHYSICS_SCREW_HPP_
#define ROLLER_PHYSICS_SCREW_HPP_

#include "pose.hpp"
#include "vector.hpp"

namespace roller {

struct ScrewM {
  /// velocity
  v::DVec<3> velo;
  /// angular velocity
  v::DVec<3> omega;
  /// center of rotation (which moves with velo, btw)
  v::DVec<3> center;
};

ScrewM mult(double t, const ScrewM &pm) {
  return {t * pm.velo, t * pm.omega, pm.center};
}

/// straight-up gives you a quaternion usable for doing the rotational part
v::DVec<4> getRotQuat(const ScrewM &pm) {
  return getRotQuaternion(pm.omega / 2);
}

void apply2Vec(v::DVec<3> &p, const v::DVec<4> &orot) {
  p = applyQuaternion(p, orot);
}
void apply2Vec(v::DVec<3> &p, const ScrewM &pm) {
  apply2Vec(p, getRotQuat(pm));
}

void apply2Pos(v::DVec<3> &p, const v::DVec<3> &velo, const v::DVec<4> &orot,
               const v::DVec<3> &center) {
  p = velo + center + applyQuaternion(p - center, orot);
}
void apply2Pos(v::DVec<3> &p, const ScrewM &pm) {
  apply2Pos(p, pm.velo, getRotQuat(pm), pm.center);
}

void apply2Quat(v::DVec<4> &q, const v::DVec<4> &orot) {
  q = normalizeQuaternion(quaternionMult(orot, q));
}
void apply2Quat(v::DVec<4> &q, const ScrewM &pm) {
  apply2Quat(q, getRotQuat(pm));
}

void apply2Pose(Pose &p, const ScrewM &pm) {
  v::DVec<4> orot = getRotQuat(pm);
  apply2Pos(p.p, pm.velo, orot, pm.center);
  apply2Quat(p.q, orot);
}

} // namespace roller

#endif // ROLLER_PHYSICS_SCREW_HPP_
