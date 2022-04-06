#ifndef GENERIC_OBJ_HPP_
#define GENERIC_OBJ_HPP_

#include "world.hpp"

namespace roller {

/*
  Although world.hpp doesn't require these methods, this header expects Prim to
  additionally have the setPose, getPhysInfo, and getSurfaceDetail methods.
  btw, getSurfaceDetail returns v::DVec<3> with {coefficient of restitution,
  static friction, kinetic friction}.
*/

template <class Prim> struct GenericPrimWrapper : public Prim {
  Pose relP;
  GenericPrimWrapper() {}
  GenericPrimWrapper(Prim &&prim) : Prim(std::move(prim)) {}
};

struct GenericObj {
  GenericObj() {}

  template <class Iterator>
  GenericObj(const Iterator &primIter, const Iterator &primIterEnd,
             CPhysInfo &ocpi) {
    CPhysInfo tmp;
    auto lambda =
        [&tmp](const typename std::iterator_traits<Iterator>::value_type &prim)
        -> CPhysInfo & {
      tmp = prim.getPhysInfo();
      return tmp;
    };
    ocpi = combinedCPI(IteratorWrapper<Iterator, decltype(lambda)>(
                           primIter, std::move(lambda)),
                       IteratorWrapper<Iterator, decltype(lambda)>(
                           primIterEnd, std::move(lambda)));
    Iterator iter = primIter;
    do {
      iter->relP = ocpi.pi.pose.fromWorldCoords(iter->getPose());
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

/// XXX: make a more interesting calculation
template <class Prim>
v::DVec<3> calcSurfaceDetail(const Prim &ap, const Prim &bp) {
  v::DVec<3> sd1 = ap.getSurfaceDetail();
  v::DVec<3> sd2 = bp.getSurfaceDetail();
  return v::elementwiseMax(sd1, sd2);
}

template <class Prim>
v::DVec<3> genericImpulseCalc(const ccd::Contact &con, const Prim &ap,
                              const CPhysInfo &api, const Prim &bp,
                              const CPhysInfo &bpi) {
  v::DVec<3> surfaceDetail = calcSurfaceDetail(ap, bp);
  v::DVec<3> pos = con.p;
  v::DVec<3> normal = con.n;
  v::DVec<3> urel = api.getVelocity(pos) - bpi.getVelocity(pos);
  double ureln = v::dot(urel, normal);
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
  if (con.d > 1e-2) {
    // penalty, should never happen
    imp += 2e-2 * normal;
  }
  return imp;
}

} // namespace roller

#endif // GENERIC_OBJ_HPP_
