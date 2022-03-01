#ifndef GENERIC_OBJ_HPP_
#define GENERIC_OBJ_HPP_

#include "world.hpp"

namespace roller {

/// requires Prim::setPose and Prim::getPhysInfo
template <class Prim> struct GenericPrimWrapper {
  Prim prim;
  Pose relP;

  void setPose(const Pose &pose) { prim.setPose(pose); }
  PhysInfo getPhysInfo() const { return prim.getPhysInfo(); }

  Pose getPose() const { return prim.getPose(); }
  AABB getAABB(const ScrewM &sm) const { return prim.getAABB(sm); }
  Contact doCCD(double t1, double t2, const ScrewM &sm1,
                const GenericPrimWrapper<Prim> &prim2,
                const ScrewM &sm2) const {
    return prim.doCCD(t1, t2, sm1, prim2.prim, sm2);
  }
};

struct GenericObj {
  template <class Iterator>
  GenericObj(const Iterator &primIter, const Iterator &primIterEnd, CPhysInfo &ocpi) {
    CPhysInfo tmp;
    auto lambda = [&tmp](const Iterator &iter) -> CPhysInfo & {
      tmp = iter->getPhysInfo();
      return tmp;
    };
    ocpi = combinedCPI(IteratorWrapper<Iterator, decltype(lambda)>(
                                     primIter, std::move(lambda)),
                                 IteratorWrapper<Iterator, decltype(lambda)>(
                                     primIterEnd, std::move(lambda)));
    Iterator iter = primIter;
    do {
      iter->relP = ocpi.pi.pose.p.fromWorldCoords(iter->getPose());
    } while (++iter != primIterEnd);
  }

  template <class Iterator>
  void setPose(const Pose &pose, const Iterator &primIter) const {
    Iterator iter = primIter;
    do {
      iter->setPose(pose.toWorldCoords(iter->relP));
    } while (++iter != primIter);
  }

  template <class Iterator> void addPrims(const Iterator &primIter) const {
    // XXX: implement changeable objects at some point plz
  }
};

} // namespace roller

#endif // GENERIC_OBJ_HPP_
