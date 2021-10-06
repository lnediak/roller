#ifndef ROLLER_PHYSICS_SOLVER_HPP_
#define ROLLER_PHYSICS_SOLVER_HPP_

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
  std::vector<Collision> c; // collisions

public:
  Solver(W &&w, double g) : w(std::forward<W>(w)) {}

  void addCollision(int i, int j) {
    auto &o1 = w.getObj(i);
    auto &o2 = w.getObj(j);
    Contact c = o1.getContact(o2);
    v::DVec<3> v1 = p[i].getVelocity(x[i], c.p);
    v::DVec<3> v2 = p[j].getVelocity(x[j], c.p);
    if (v::dot(v1 - v2, c.n) < 1e-5) {
      // TODO: we have a collision
    }
  }

  // XXX: Adaptive timestep duration
  void step(double dt) {
    std::size_t numObjs = w.numObjs();
    p.resize(numObjs);
    for (std::size_t i = 0; i < numObjs; i++) {
      p[i] = w.getObj(i).getPhysInfo();
      x[i] = p[i].getAuxInfo();
      p[i].stepTime(x[i], dt);
      x[i] = p[i].getAuxInfo();
    }
  }
};

} // namespace roller

#endif // ROLLER_PHYSICS_SOLVER_HPP_

