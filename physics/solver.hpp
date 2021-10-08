#ifndef ROLLER_PHYSICS_SOLVER_HPP_
#define ROLLER_PHYSICS_SOLVER_HPP_

#include <cmath>
#include <utility>

#include "phys_info.hpp"

namespace roller {

/**
  W template: All this needs two functions:
  std::size_t numObjs() const, and
  Obj &getObj(std::size_t).

  Obj template: This needs a little more:
  PhysInfo getPhysInfo() const,
  void setPhysInfo(const PhysInfo &pi),
  AABB getAABB() const, and
  Contact getContact(Obj other) const.
*/
template <class W> class Solver {

  World<W> w;
  double g; // acceleration of gravity

  struct Collision {
    Contact c;
    int i, j;
  };

  // temporaries:
  std::vector<PhysInfo> p;
  std::vector<AuxPhysInfo> x;
  std::vector<Collision> cols; // collisions

  double eps; // collision threshold

public:
  Solver(W &&w, double g) : w(std::forward<W>(w)) {}

  void addCollision(int i, int j) {
    auto &o1 = w.getObj(i);
    auto &o2 = w.getObj(j);
    Contact c = o1.getContact(o2);
    if (c.dist > 0) {
      return;
    }
    v::DVec<3> v1 = p[i].getVelocity(x[i], c.p);
    v::DVec<3> v2 = p[j].getVelocity(x[j], c.p);
    double ureln = v::dot(v1 - v2, c.n);
    if (ureln < eps) {
      cols.push_back({c, i, j, v1, v2, ureln});
    }
  }

  void processCollision(const Collision &col) {
    int i = col.i;
    int j = col.j;
    v::DVec<3> p = col.c.p;
    v::DVec<3> n = col.c.n;
    v::DVec<3> urel = p[i].getVelocity(x[i], c.p) - p[j].getVelocity(x[j], c.p);
    double ureln = v::dot(urel, n);
    if (ureln >= eps) {
      // other collisions changed this
      return;
    }
    double el = std::min(p[i].el, p[j].el);
    v::DVec<3> lhs = -el * ureln * n - urel;
    double mi = p[i].massi + p[j].massi;
    v::DVec<3> r1 = p - p[i].pose.p;
    DMat3x3 rcr1 = crossMat(r1);
    v::DVec<3> r2 = p - p[j].pose.p;
    DMat3x3 rcr2 = crossMat(r2);
    DMat3x3 kt = diagMat(mi) + rcr1.transpose() * p[i].ineri * rcr1 +
                 rcr2.transpose() * p[j].ineri * rcr2;
    v::DVec<3> j = kt.solve(lhs);
    // test if in friction cone
    // XXX: MODIFY STATIC AND KINETIC FRICTION FORMULAS BELOW
    double jdn = v::dot(j, n);
    double sf = (p[i].sf + p[j].sf) / 2.;
    if (v::norm2(j - jdn * n) <= jdn * jdn * sf * sf) {
      // in cone; update velocities
      p[i].lm += j;
      p[j].lm -= j;
      p[i].am += cross3(r1, j);
      p[j].am -= cross3(r2, j);
      x[i] = p[i].getAuxInfo();
      x[j] = p[j].getAuxInfo();
      w.getObj(i).setPhysInfo(p[i]);
      w.getObj(j).setPhysInfo(p[j]);
      return;
    }
    double kf = (p[i].kf + p[j].kf) / 2.;
    v::DVec<3> utann = urel / std::sqrt(v::norm2(urel));
    jdn = -(el + 1) * ureln / v::dot(n, kt * (n - kf * utann));
    j = jdn * n;

    p[i].lm += j;
    p[j].lm -= j;
    p[i].am += cross3(r1, j);
    p[j].am -= cross3(r2, j);
    x[i] = p[i].getAuxInfo();
    x[j] = p[j].getAuxInfo();
    w.getObj(i).setPhysInfo(p[i]);
    w.getObj(j).setPhysInfo(p[j]);
  }

  // XXX: Adaptive timestep duration
  void step(double dt) {
    cols.clear();
    std::size_t numObjs = w.numObjs();
    p.resize(numObjs);
    for (std::size_t i = 0; i < numObjs; i++) {
      p[i] = w.getObj(i).getPhysInfo();
      x[i] = p[i].getAuxInfo();
      p[i].stepTime(x[i], dt);
      x[i] = p[i].getAuxInfo();
    }
    // collision resolution passes
    for (int spam = 0; spam < 10; spam++) {
      w.exportAllAABBInts(
          [this](int i, int j) -> void { addCollision(i, j, -1e-3); });
      if (cols.empty()) {
        break;
      }
      for (const Collision &col : cols) {
        processCollision(col);
      }
    }
    // TODO: CONTACT HANDLING
    for (std::size_t i = 0; i < numObjs; i++) {
      p[i].stepTime(x[i], dt);
    }
  }
};

} // namespace roller

#endif // ROLLER_PHYSICS_SOLVER_HPP_

