#ifndef ROLLER_PHYSICS_AABB_HPP_
#define ROLLER_PHYSICS_AABB_HPP_

#include "vector.hpp"

namespace roller {

struct AABB {
  v::DVec<3> m[2]; /// min, max

  bool intersects(const AABB &c) const {
    return (m[0][0] < c.m[1][0]) && (m[1][0] > c.m[0][0]) &&
           (m[0][1] < c.m[1][1]) && (m[1][1] > c.m[0][1]) &&
           (m[0][2] < c.m[1][2]) && (m[1][2] > c.m[0][2]);
  }
};

} // namespace roller

#endif // ROLLER_PHYSICS_AABB_HPP_

