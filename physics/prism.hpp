#ifndef ROLLER_PHYSICS_PRISM_HPP_
#define ROLLER_PHYSICS_PRISM_HPP_

#include "pose.hpp"

template <class Tag> struct Prism {

  v::IVec<3> p; /// lowest corner
  v::IVec<3> s; /// sidelengths

  Pose pose;
};

#endif // ROLLER_PHYSICS_PRISM_HPP_

