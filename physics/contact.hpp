#ifndef ROLLER_PHYSICS_CONTACT_HPP_
#define ROLLER_PHYSICS_CONTACT_HPP_

#include "vector.hpp"

namespace roller {

namespace ccd {

struct Contact {
  double t;     /// time of intersection
  double d;     /// depth of penetration
  v::DVec<3> p; /// contact point
  v::DVec<3> n; /// contact normal (pointing towards first object)
};

} // namespace ccd

} // namespace roller

#endif // ROLLER_PHYSICS_CONTACT_HPP_

