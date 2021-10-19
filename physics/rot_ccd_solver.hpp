#ifndef ROLLER_PHYSICS_ROT_CCD_SOLVER_HPP_
#define ROLLER_PHYSICS_ROT_CCD_SOLVER_HPP_

#include <cmath>
#include <vector>

#include "coalesced_obj.hpp"
#include "phys_info.hpp"

namespace roller {

bool doCollisionResponse(const Contact &con, PhysInfo &api, PhysInfo &bpi,
                         AuxPhysInfo &aaux, AuxPhysInfo &baux,
                         const v::DVec<3> &surfaceDetail) {
  PhysInfo aopi = api;
  PhysInfo bopi = bpi;
  api.stepTime<false>(aaux, con.t);
  bpi.stepTime<false>(baux, con.t);
  v::DVec<3> pos = con.p;
  v::DVec<3> normal = con.n;
  v::DVec<3> urel = api.getVelocity(aaux, pos) - bpi.getVelocity(baux, pos);
  double ureln = v::dot(urel, normal);
  if (ureln > 0) {
    api = aopi;
    bpi = bopi;
    return false;
  }
  aaux = api.getAuxInfo();
  baux = bpi.getAuxInfo();

  double el = surfaceDetail[0]; // elasticity (coefficient of restitution)
  v::DVec<3> lhs = -el * ureln * normal - urel;
  double mi = api.massi + bpi.massi;
  v::DVec<3> r1 = pos - api.pose.p;
  DMat3x3 rcr1 = crossMat(r1);
  v::DVec<3> r2 = pos - bpi.pose.p;
  DMat3x3 rcr2 = crossMat(r2);
  DMat3x3 kt = diagMat(mi, mi, mi) + rcr1.transpose() * aaux.riner * rcr1 +
               rcr2.transpose() * baux.riner * rcr2;
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
    imp += 1e-3 * normal;
  }
  api.lm += imp;
  bpi.lm -= imp;
  api.am += cross3(r1, imp);
  bpi.am -= cross3(r2, imp);
  aaux = api.getAuxInfo();
  baux = bpi.getAuxInfo();
  return true;
}

/**
  Primary structure for using the physics stuff here. Stores objects (you can
  always make a wrapper class if you want to store the object elsewhere).

  Obj template: Refers to an object. Must have the following:

  AABB obj.getAABB() const
    - gets an AABB for the object
  v::DVec<3> getSurfaceDetail(const Obj &o, const v::DVec<3> &p) const
    - returns {elasticity, static friction, kinetic friction}
  PhysInfo obj.physInfo() const
    - returns this object's physInfo
  void obj.setPhysInfo(const PhysInfo &pi)
    - sets this object's physInfo

  and the hardest one (needs to handle screw motion):

  template <class Fun1, class Fun2>
  Contact doCCD(double velo1, double omega1, v::DVec<3> c1, const Obj &other,
  double velo2, double omega2, v::DVec<3> c2)
    - returns a Contact where 0<=c.t<=1, and for no contact, c.t can be any
      value >1. *this is considered to be travelling in screw motion with
      primary velocity velo1, and with circular portion around
      pose.p+toR1+t*velo1 at constant angular velocity omega1, while other going
      velo2 and is around other.pose.p+toR2+t*velo2 at angular velocity omega2
*/
template <class Obj> class RotCCDSolver {

  std::vector<Obj> objs;
  std::vector<CoalescedObj<Obj>> cobjs;
  /// has same sz as objs, if positive is index in cobjs
  std::vector<int> trueObjs;

  double g; // acceleration of gravity

  struct Collision {
    Contact c;
    int i, j;
  };

  // temporaries:
  // same as trueObjs, but with the temporary coalescing
  std::vector<int> tmpObjs;
  std::vector<PhysInfo> p;
  std::vector<AuxPhysInfo> x;
  std::unordered_set<Collision> cols; // collisions

public:
  RotCCDSolver(W &&w, double g) : w(std::forward<W>(w)), g(g) {}

  void checkCollision(double dt, int i, int j) {
    if (!p[i].mass && !p[j].mass) {
      return;
    }
    auto &o1 = objs[i];
    auto &o2 = objs[j];
    Contact c = o1.doCCD(o2, dt);
    if (c.t > dt) {
      return;
    }
    cols.push_back({c, i, j});
  }

  /// returns true when no collisions
  bool processCollisions(double &time) {
    std::vector<decltype(cols.cend())> erasureList;
    for (auto iter = cols.cbegin(), iter_end = cols.cend(); iter != iter_end;
         ++iter) {
      const Collision &col = *iter;
      int i = col.i;
      int j = col.j;
      if (doCollisionResponse(col.c, p[i], p[j], x[i], x[j],
                              objs[i].getSurfaceDetail(objs[j]))) {
        time += col.c.t;
        for (std::size_t ii = 0, sz = objs.size(); ii < sz; ii++) {
          if (ii != i && ii != j) {
            p[ii].stepTime<true>(x[ii], col.c.t);
          }
        }
        return false;
      } else {
        erasureList.push_back(iter);
      }
    }
    for (const auto &ii : erasureList) {
      cols.erase(ii);
    }
    return true;
  }

  // FIXME: REWRITE
  void step(double dt) {
    cols.clear();
    // collision resolution passes
    for (int spam = 0; spam < 5; spam++) {
      for (std::size_t i = 0; i < numObjs; i++) {
        p[i].pose = porig[i].pose;
        x[i] = p[i].getAuxInfo();
        p[i].stepTime(x[i], dt, g);
        x[i] = p[i].getAuxInfo();
        w.getObj(i).setPhysInfo(p[i]);
      }
      w.exportAllAABBInts([this](int i, int j) -> void { addCollision(i, j); });
      if (cols.empty()) {
        break;
      }
      for (const Collision &col : cols) {
        processCollision(col, 1);
      }
    }
    for (std::size_t i = 0; i < numObjs; i++) {
      p[i].pose = porig[i].pose;
      p[i].stepVelo(dt, g);
    }
    // contact resolution passes
    for (int spam = 0; spam < 5; spam++) {
      for (std::size_t i = 0; i < numObjs; i++) {
        p[i].pose = porig[i].pose;
        x[i] = p[i].getAuxInfo();
        p[i].stepTime(x[i], dt, 0);
        x[i] = p[i].getAuxInfo();
        w.getObj(i).setPhysInfo(p[i]);
      }
      w.exportAllAABBInts([this](int i, int j) -> void { addCollision(i, j); });
      if (cols.empty()) {
        break;
      }
      for (const Collision &col : cols) {
        std::cout << "COLLISION?????" << std::endl;
        processCollision(col, 0);
      }
    }
    for (std::size_t i = 0; i < numObjs; i++) {
      p[i].pose = porig[i].pose;
      x[i] = p[i].getAuxInfo();
      p[i].stepTime(x[i], dt, 0);
      w.getObj(i).setPhysInfo(p[i]);
    }
  }
};

} // namespace roller

#endif // ROLLER_PHYSICS_ROT_CCD_SOLVER_HPP_

