#ifndef ROLLER_PHYSICS_ROT_CCD_SOLVER_HPP_
#define ROLLER_PHYSICS_ROT_CCD_SOLVER_HPP_

#include <cmath>
#include <unordered_map>
#include <vector>

#include "coalesced_obj.hpp"
#include "phys_info.hpp"

namespace roller {

bool doCollisionResponse(const Contact &con, CPhysInfo &api, CPhysInfo &bpi,
                         const v::DVec<3> &surfaceDetail) {
  PhysInfo aopi = api.pi;
  PhysInfo bopi = bpi.pi;
  api.stepTime<false>(con.t);
  bpi.stepTime<false>(con.t);
  v::DVec<3> pos = con.p;
  v::DVec<3> normal = con.n;
  v::DVec<3> urel = api.getVelocity(aaux, pos) - bpi.getVelocity(baux, pos);
  double ureln = v::dot(urel, normal);
  if (ureln > 0) {
    api.pi = aopi;
    bpi.pi = bopi;
    return false;
  }
  api.updateAux();
  bpi.updateAux();

  double el = surfaceDetail[0]; // elasticity (coefficient of restitution)
  v::DVec<3> lhs = -el * ureln * normal - urel;
  double mi = api.pi.massi + bpi.pi.massi;
  v::DVec<3> r1 = pos - api.pi.pose.p;
  DMat3x3 rcr1 = crossMat(r1);
  v::DVec<3> r2 = pos - bpi.pi.pose.p;
  DMat3x3 rcr2 = crossMat(r2);
  DMat3x3 kt = diagMat(mi, mi, mi) + rcr1.transpose() * api.aux.riner * rcr1 +
               rcr2.transpose() * bpi.aux.riner * rcr2;
  v::DVec<3> imp = kt.solve(lhs);
  double sf = surfaceDetail[1]; // static friction
  double jdn = v::dot(imp, normal);
  if (v::norm2(imp - jdn * normal) > jdn * jdn * sf * sf) {
    // not in friction cone; use regular response with kinetic friction
    double kf = surfaceDetail[2]; // kinetic friction
    v::DVec<3> utann = urel / std::sqrt(v::norm2(urel));
    jdn = -(el + 1) * ureln / v::dot(normal, kt * (normal - kf * utann));
    imp = jdn * normal;
  }
  if (con.depth > 1e-2) {
    // penalty, should never happen
    imp += 1e-3 * normal;
  }
  api.pi.lm += imp;
  bpi.pi.m -= imp;
  api.pi.am += cross3(r1, imp);
  bpi.pi.am -= cross3(r2, imp);
  api.updateAux();
  bpi.updateAux();
  return true;
}

/**
  Primary structure for using the physics stuff here. Stores objects (you can
  always make a wrapper class if you want to store the object elsewhere).

  Obj template: Refers to an object. Must have the following:

  v::DVec<3> getSurfaceDetail(const Obj &o, const v::DVec<3> &p) const
    - returns {elasticity, static friction, kinetic friction}
  PhysInfo obj.physInfo() const
    - returns this object's physInfo
  void obj.setPhysInfo(const PhysInfo &pi)
    - sets this object's physInfo

  Contact doCCD(v::DVec<3> velo1, v::DVec<3> omega1, v::DVec<3> c1, const Obj
  &other, v::DVec<3> velo2, v::DVec<3> omega2, v::DVec<3> c2) const
    - returns a Contact where 0<=c.t<=1, and for no contact, c.t can be any
      value >1. *this is considered to be travelling in screw motion with
      primary velocity velo1, and with circular portion around c1+t*velo1 at
      constant angular velocity omega1, while other going velo2 and is around
      c2+t*velo2 at angular velocity omega2
  AABB obj.getAABB(v::DVec<3> velo, v::DVec<3> omega, v::Dvec<3> c) const
    - gets an AABB for the object in screw motion (as above)
*/
template <class Obj> class RotCCDSolver {

  std::vector<Obj> objs;
  UnionFind coalesces;

  double g; // acceleration of gravity

  struct AABBProvider {
    RotCCDSolver &s;
    double dt;
    std::size_t numObjs() const { return s.objs.size(); }
    AABB getAABB(std::size_t i) const {
      PhysInfo pi = s.getBasePhysInfo(i);
      AuxPhysInfo aux = pi.getAuxInfo();
      return objs[i].getAABB(dt * aux.velo, dt * aux.omega, pi.pose.p);
    }
  };
  AABBProvider bProvider;
  BroadPhaseAABB<AABBProvider &> broadPhase;

  struct Collision {
    Contact c;
    int i, j;
  };

  // temporaries:
  std::unordered_map<int, CoalescedObj> cobjs;
  // indices of objs to collisions, representing earliest-time collision
  std::unordered_map<int, Collision> cols;
  // same as coalesces, but with the temporary coalescing
  UnionFind tmpCoalesces;

public:
  RotCCDSolver(std::vector<Obj> &&objs, double g)
      : objs(objs), coalesces(objs.size()),
        g(g), bProvider{*this, 1. /*arbitrary*/}, broadPhase(bProvider) {}

  /// for arbitrary i, call on tmpCoalesces.find(i)
  CPhysInfo &basePhysInfo(int i) {
    return tmpCoalesces.getSz(i) > 1 ? cpis[i] : objs[i];
  }

  /// this handles the possibility of collision between objs[i] and objs[j]
  void checkCollision(double time, double dt, bool doC, int i, int j) {
    int irep = tmpCoalesces.find(i);
    int jrep = tmpCoalesces.find(j);
    if (irep == jrep) {
      // part of same coalesced object
      return;
    }
    double idt = dt - time, jdt = dt - time;
    CPhysInfo cpii, cpij;
    if (tmpCoalesces.getSz(irep) > 1) {
      auto it = cols.find(irep);
      if (it != cols.end()) {
        idt = it->c.t;
      }
      cpii = cobjs[irep].pi;
    } else {
      // well, i == irep anyways in this case...
      cpii = objs[irep].pi;
    }
    if (!cpii.mass) {
      return;
    }
    if (tmpCoalesces.getSz(jrep) > 1) {
      auto it = cols.find(jrep);
      if (it != cols.end()) {
        jdt = it->c.t;
      }
      cpij = cobjs[jrep].pi;
    } else {
      cpij = objs[jrep].pi;
    }
    if (!cpij.mass) {
      return;
    }
    dt = idt < jdt ? idt : jdt;
    Contact c = objs[i].doCCD(dt * cpii.aux.velo, dt * cpii.aux.omega,
                              cpii.pi.pose.p, objs[j], dt * cpij.aux.velo,
                              dt * cpii.aux.omega, cpii.pi.pose.p);
    // TODO: NEED TO CONSIDER WHETHER THE COLLISION IS LEAVING
    if (c.t > 1) {
      return;
    }
    c.t *= dt;
    cols[i] = {c, i, j};
    cols[j] = {c, i, j};
  }

  /// returns true when no collisions
  bool processCollision(double &time) {
    for (auto iter = cols.cbegin(), iter_end = cols.cend(); iter != iter_end;
         ++iter) {
      const Collision &col = *iter;
      int i = col.i;
      int j = col.j;
      PhysInfo pii = getBasePhysInfo(i);
      PhysInfo pij = getBasePhysInfo(j);
      AuxPhysInfo auxi = pii.getAuxInfo(), auxj = pij.getAuxInfo();
      if (doCollisionResponse(col.c, pii, pij, auxi, auxj,
                              objs[i].getSurfaceDetail(objs[j], col.c.p))) {
        time += col.c.t;
        for (std::size_t ii = 0, sz = objs.size(); ii < sz; ii++) {
          if (ii != i && ii != j && tmpObjs[ii] < 0) {
            PhysInfo tmp = objs[ii].physInfo();
            AuxPhysInfo taux = tmp.getAuxInfo();
            tmp.stepTime<false>(taux, col.c.t);
            objs[ii].setPhysInfo(tmp);
          }
        }
        int iin = tmpObjs[i];
        int jjn = tmpObjs[j];
        for (std::size_t ii = 0, sz = cobjs.size(); ii < sz; ii++) {
          if (ii != iin && ii != jjn) {
            PhysInfo tmp = cobjs[ii].physInfo();
            AuxPhysInfo taux = tmp.getAuxInfo();
            tmp.stepTime<false>(taux, col.c.t);
            cobjs[ii].setPhysInfo(tmp);
          }
        }
        // XXX: do this more efficiently
        std::vector<Obj *> ccobjs; // for passing to CoalescedObj
        if (tmpObjs[i] < 0) {
          ccobjs.push_back(&objs[i]);
        } else {
          for (Obj *ptr : cobjs[tmpObjs[i]].listObjs()) {
            ccobjs.push_back(ptr);
          }
        }
        if (tmpObjs[j] < 0) {
          ccobjs.push_back(&objs[j]);
        } else {
          for (Obj *ptr : cobjs[tmpObjs[j]].listObjs()) {
            ccobjs.push_back(ptr);
          }
        }
        std::size_t ind = cobjs.size();
        cobjs.emplace_back(ccobjs);
        tmpObjs[i] = ind;
        tmpObjs[j] = ind;
        return false;
      }
    }
    return true;
  }

  void step(double dt) {
    cols.clear();
    tmpObjs = trueObjs;
    double time = 0;
    while (time < dt) {
      bProvider.dt = dt - time;
      broadPhase.exportAllAABBInts(
          [this](int i, int j) -> void { addCollision(i, j); });
      if (processCollisions(time)) {
        break;
      }
    }
    if (time < dt) {
      for (std::size_t i = 0, sz = objs.size(); i < sz; i++) {
        if (tmpObjs[i] < 0) {
          PhysInfo tmp = objs[i].physInfo();
          AuxPhysInfo taux = tmp.getAuxInfo();
          tmp.stepTime<false>(taux, dt - time);
          objs[i].setPhysInfo(tmp);
        }
      }
      for (std::size_t i = 0, sz = cobjs.size(); i < sz; i++) {
        PhysInfo tmp = cobjs[i].physInfo();
        AuxPhysInfo taux = tmp.getAuxInfo();
        tmp.stepTime<false>(taux, dt - time);
        cobjs[i].setPhysInfo(tmp);
      }
    }
    cobjs.resize(cobjsSz);
  }
};

} // namespace roller

#endif // ROLLER_PHYSICS_ROT_CCD_SOLVER_HPP_

