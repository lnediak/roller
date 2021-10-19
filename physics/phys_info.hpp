#ifndef ROLLER_PHYSICS_PHYS_INFO_HPP_
#define ROLLER_PHYSICS_PHYS_INFO_HPP_

#include <algorithm>
#include <cmath>
#include <unordered_set>

#include "aabb.hpp"
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

} // namespace roller

#endif // ROLLER_PHYSICS_PHYS_INFO_HPP_

