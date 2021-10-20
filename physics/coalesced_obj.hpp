#ifndef ROLLER_PHYSICS_COALESCED_OBJ_HPP_
#define ROLLER_PHYSICS_COALESCED_OBJ_HPP_

#include <vector>

#include "phys_info.hpp"

namespace roller {

/// Obj as defined in solver.hpp.
template <class Obj> class CoalescedObj {

  std::vector<Obj *> objs;
  std::vector<Pose> poses;

  PhysInfo pi;

  /// df is the new center relative to the center of mass
  DMat3x3 adjustCenter(const DMat3x3 &iner, double mass, const v::DVec<3> &df) {
    v::DVec<3> mm = mass * df;
    double xx = mm[0] * df[0];
    double xy = mm[0] * df[1];
    double xz = mm[0] * df[2];
    double yy = mm[1] * df[1];
    double yz = mm[1] * df[2];
    double zz = mm[2] * df[2];
    DMat3x3 dblCrossMat = {
        {yy + zz, -xy, -xz}, {-xy, xx + zz, -yz}, {-xz, -yz, xx + yy}};
    return iner + dblCrossMat;
  }

  void evalPhysInfo() const {
    double tm = 0;     // total mass
    v::DVec<3> cm = 0; // center of mass
    bool fat = false;
    for (const Obj &o : objs) {
      PhysInfo opi = o->physInfo();
      double m = opi.mass;
      v::DVec<3> c = opi.pose.p;
      if (!m) {
        tm = 0;
        cm = c;
        fat = true;
        break;
      }
      tm += m;
      cm += m * c;
    }
    pi.mass = tm;
    double tmi = 0;
    if (fat) {
      pi.massi = tmi;
      pi.pose.p = cm;
    } else {
      tmi = 1 / tm;
      cm *= tmi;
      pi.massi = tmi;
      pi.pose.p = cm;
    }
    v::DVec<3> tlm = {0, 0, 0};
    DMat3x3 iner = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
    v::DVec<3> tam = {0, 0, 0};
    for (const Obj &o : objs) {
      PhysInfo opi = o->physInfo();
      poses.push_back(opi.pose);
      v::DVec<3> relp = poses.back().p -= cm;
      if (!fat) {
        iner += adjustCenter(opi.getAuxInfo().riner, opi.mass, -relp);
        tlm += opi.lm;
        tam += opi.am + cross3(relp, opi.lm);
      }
    }
    pi.lm = tlm;
    pi.am = tam;
    pi.iner = iner;
    if (fat) {
      pi.ineri = iner;
    } else {
      pi.ineri = iner.inverse();
    }
  }

public:
  CoalescedObj(const std::vector<Obj *> &objs) : objs(objs) { evalPhysInfo(); }

  PhysInfo getPhysInfo() const { return pi; }
  void updatePhysInfo0(PhysInfo &opi, const Pose &p) const {
    opi.pose.p = p.p + pi.pose.p;
    opi.pose.q = quaternionMult(pi.pose.q, p.q);
  }
  void updatePhysInfo(std::size_t i) {
    PhysInfo opi = objs[i].physInfo();
    updatePhysInfo0(opi, poses[i]);
    objs[i].setPhysInfo(opi);
  }
  void setPhysInfo(const PhysInfo &ppi) {
    if (v::norm2(pi.pose.p - ppi.pose.p) + v::norm2(pi.pose.q - ppi.pose.q) >=
        1e-16) {
      pi = ppi;
      for (std::size_t i = 0, sz = objs.size(); i < sz; i++) {
        updatePhysInfo(i);
      }
      return;
    }
    pi = ppi;
  }
  // note: no doCCD, as that basically bypasses broad-phase, which is not good

  void addObj(Obj *obj) {
    PhysInfo opi = obj->physInfo();
    objs.push_back(obj);
    poses.push_back(opi.pose);
    if (!pi.mass) {
      poses.back().p -= pi.pose.p;
      return;
    }
    if (!opi.mass) {
      pi.mass = pi.massi = 0;
      pi.pose.p = opi.pose.p;
      poses.back().p = {0, 0, 0};
      return;
    }
    pi.lm += opi.lm;
    v::DVec<3> origP = pi.pose.p;
    double origMass = pi.mass;
    pi.mass += opi.mass;
    pi.massi = 1 / pi.mass;

    pi.pose.p = (origMass * pi.pose.p + opi.mass * opi.pose.p) * pi.massi;
    v::DVec<3> relp = opi.pose.p - pi.pose.p;
    pi.am += cross3(origP - pi.pose.p, pi.lm) + opi.am + cross3(relp, opi.lm);
    relp = poses.back().p = pi.pose.fromShiftWorldCoords(relp);
    // this line is for the getAuxInfo to evaluate in the object space of *this
    opi.pose.q = poses.back().q =
        quaternionMult(quaternionConj(pi.pose.q), opi.pose.q);

    pi.iner = adjustCenter(pi.iner, origMass,
                           pi.pose.fromShiftWorldCoords(pi.pose.p - origP)) +
              adjustCenter(opi.getAuxInfo().riner, opi.mass, -relp);
    pi.ineri = pi.iner.inverse();
  }

  void removeObj(std::size_t i) {
    PhysInfo opi = objs[i].physInfo();
    objs.erase(objs.begin() + i);
    if (opi.massi < 1e-8) {
      // in this case the stuff is unstable
      evalPhysInfo();
      return;
    }
    if (!pi.mass) {
      return;
    }
    pi.lm -= opi.lm;
    v::DVec<3> origP = pi.pose.p;
    double fatMass = pi.mass;
    pi.mass -= opi.mass;
    pi.massi = 1 / pi.mass;

    pi.pose.p = (fatMass * pi.pose.p - opi.mass * opi.pose.p) * pi.massi;
    pi.am += cross3(origP - pi.pose.p, pi.lm) - opi.am -
             cross3(pi.pose.toShiftWorldCoords(opi.pose.p), opi.lm);

    pi.iner = adjustCenter(
        pi.iner - adjustCenter(opi.getAuxInfo().riner, opi.mass, -opi.pose.p),
        pi.mass, pi.pose.fromShiftWorldCoords(pi.pose.p - origP));
    pi.ineri = pi.iner.inverse();
  }
  void removeObj(Obj *obj) { removeObj(objs.find(obj) - objs.begin()); }
};

} // namespace roller

#endif // ROLLER_PHYSICS_COALESCED_OBJ_HPP_

