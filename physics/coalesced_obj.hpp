#ifndef ROLLER_PHYSICS_COALESCED_OBJ_HPP_
#define ROLLER_PHYSICS_COALESCED_OBJ_HPP_

#include <vector>

#include "phys_info.hpp"

namespace roller {

template <class Obj> class CoalescedObj {

  std::vector<Obj> objs;

  PhysInfo pi;

  /// df is the new center relative to the center of mass
  DMat3x3 adjustCenter(const DMat3x3 &iner, double mass, const v::DVec<3> &df) {
    double mm = mass * df;
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

public:
  CoalescedObj(const std::vector<Obj> &objs) : objs(objs) {
    double tm = 0;     // total mass
    v::DVec<3> cm = 0; // center of mass
    bool fat = false;
    for (const Obj &o : objs) {
      PhysInfo &opi = o.getPhysInfo();
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
    if (!fat) {
      for (const Obj &o : objs) {
        PhysInfo &opi = o.getPhysInfo();
        tlm += opi.lm;
        DMat3x3 niner = adjustCenter(opi.iner, opi.mass, cm - opi.pose.p);
        tam -= cross3(cm - opi.pose.p, opi.lm).
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

  //
};

} // namespace roller

#endif // ROLLER_PHYSICS_COALESCED_OBJ_HPP_

