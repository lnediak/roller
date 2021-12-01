#ifndef ROLLER_PHYSICS_SCREW_HPP_
#define ROLLER_PHYSICS_SCREW_HPP_

#include "vector.hpp"

struct ScrewM {
  /// velocity
  v::DVec<3> velo;
  /// angular velocity
  v::DVec<3> omega;
  /// center of rotation (which moves with velo, btw)
  v::DVec<3> center;
};

#endif // ROLLER_PHYSICS_SCREW_HPP_
