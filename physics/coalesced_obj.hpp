#ifndef ROLLER_PHYSICS_COALESCED_OBJ_HPP_
#define ROLLER_PHYSICS_COALESCED_OBJ_HPP_

#include <vector>

#include "phys_info.hpp"

namespace roller {

class CoalescedObj {

  std::vector<CPhysInfo *> objs;
  std::vector<Pose> poses;

  CPhysInfo cpi;

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

  /// evaluates and writes into pi and poses
  void evalPhysInfo() {
    double tm = 0;             // total mass
    v::DVec<3> cm = {0, 0, 0}; // center of mass
    bool fat = false;
    for (const CPhysInfo *o : objs) {
      double m = o->pi.mass;
      v::DVec<3> c = o->pi.pose.p;
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
    for (const CPhysInfo *o : objs) {
      poses.push_back(o->pi.pose);
      v::DVec<3> relp = poses.back().p -= cm;
      if (!fat) {
        iner += adjustCenter(o->aux.riner, o->pi.mass, -relp);
        tlm += o->pi.lm;
        tam += o->pi.am + cross3(relp, o->pi.lm);
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
  CoalescedObj(const std::vector<CPhysInfo *> &objs) : objs(objs) {
    evalPhysInfo();
  }

  CPhysInfo physInfo() const { return pi; }
  void updatePhysInfo0(CPhysInfo *o, const Pose &p) const {
    o->pi.pose.p = p.p + pi.pose.p;
    o->pi.pose.q = quaternionMult(pi.pose.q, p.q);
    o->updateAux();
  }
  void setPhysInfo(const PhysInfo &ppi) {
    if (v::norm2(pi.pose.p - ppi.pose.p) + v::norm2(pi.pose.q - ppi.pose.q) >=
        1e-16) {
      pi = ppi;
      for (std::size_t i = 0, sz = objs.size(); i < sz; i++) {
        updatePhysInfo0(objs[i], poses[i]);
      }
    }
  }

  void addObj(CPhysInfo *o) {
    objs.push_back(o);
    poses.push_back(o->pi.pose);
    if (!pi.mass) {
      poses.back().p -= pi.pose.p;
      return;
    }
    if (!o->pi.mass) {
      pi.mass = pi.massi = 0;
      pi.pose.p = o->pi.pose.p;
      poses.back().p = {0, 0, 0};
      return;
    }
    v::DVec<3> origP = pi.pose.p;
    double origMass = pi.mass;
    pi.mass += o->pi.mass;

    pi.pose.p = (origMass * pi.pose.p + o->pi.mass * o->pi.pose.p) * pi.massi;
    v::DVec<3> relOrigP = pi.pose.p - origP;
    for (Pose &pose : poses) {
      pose.p += relOrigP;
    }

    v::DVec<3> relp = o->pi.pose.p - pi.pose.p;
    pi.am +=
        cross3(origP - pi.pose.p, pi.lm) + o->pi.am + cross3(relp, o->pi.lm);
    pi.lm += o->pi.lm;
    // these update poses.back() to actually store relative pose in object space
    relp = poses.back().p = pi.pose.fromShiftWorldCoords(relp);
    poses.back().q = quaternionMult(quaternionConj(pi.pose.q), o->pi.pose.q);

    DMat3x3 rotm = poses.back().toRotationMatrix();
    if (pi.mass) {
      pi.iner =
          adjustCenter(pi.iner, origMass,
                       pi.pose.fromShiftWorldCoords(pi.pose.p - origP)) +
          adjustCenter(rotm * o->pi.iner * rotm.transpose(), o->pi.mass, -relp);
      pi.ineri = pi.iner.inverse();
    } else {
      // ineri needs to be 0, but iner simply useless in this case
      pi.iner = pi.ineri = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
    }
  }
  void addObjs(CoalescedObj *o) {
    objs.insert(objs.end(), o->objs.begin(), o->objs.end());
    v::DVec<3> origP = pi.pose.p;
    double origMass = pi.mass;
    pi.mass += o->pi.mass;
    pi.massi = (origMass && o->pi.mass) ? 1 / pi.mass : 0;

    if (origMass) {
      if (o->pi.mass) {
        pi.pose.p =
            (origMass * pi.pose.p + o->pi.mass * o->pi.pose.p) * pi.massi;
      } else {
        pi.pose.p = o->pi.pose.p;
      }
    }
    v::DVec<3> relOrigP = origP - pi.pose.p;
    if (origMass) {
      for (Pose &pose : poses) {
        pose.p += relOrigP;
      }
    }

    v::DVec<3> relp = o->pi.pose.p - pi.pose.p;
    v::DVec<3> shiftQuat =
        quaternionMult(quaternionConj(pi.pose.q), o->pi.pose.q);
    for (Pose pose : o->poses) {
      pose.p = applyQuaternion(pose.p, shiftQuat) + relp;
      pose.q = quaternionMult(shiftQuat, pose.q);
      poses.push_back(pose);
    }
    pi.am += cross3(relOrigP, pi.lm) + o->pi.am + cross3(relp, o->pi.lm);
    pi.lm += o->pi.lm;

    DMat3x3 rotm = quaternionToRotMat(shiftQuat);
    if (pi.mass) {
      pi.iner = adjustCenter(pi.iner, origMass,
                             -pi.pose.fromShiftWorldCoords(relOrigP)) +
                adjustCenter(rotm * o->pi.iner * rotm.transpose(), o->pi.mass,
                             -pi.pose.fromShiftWorldCoords(relp));
      pi.ineri = pi.iner.inverse();
    } else {
      pi.iner = pi.ineri = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
    }
  }
};

} // namespace roller

#endif // ROLLER_PHYSICS_COALESCED_OBJ_HPP_
