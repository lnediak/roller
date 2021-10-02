#ifndef ROLLER_BOOL_TERRAIN_GENERATOR_HPP_
#define ROLLER_BOOL_TERRAIN_GENERATOR_HPP_

#include "vector.hpp"

namespace roller {

struct BoolBlockdata {
  bool b;
  bool isAir() const { return !b; }
  bool isOpaque() const { return b; }
  unsigned getColor(int bb) const {
    if (!b) {
      return 0;
    }
    const unsigned results[] = {0x222222FF, 0x444444FF, 0x666666FF, 0x888888FF,
                                0xAAAAAAFF, 0xCCCCCCFF, 0xFFFFFFFF};
    return results[bb];
  }
};

template <class FTerGen> struct BoolTerrainGenerator {
  typedef BoolBlockdata value_type;

  FTerGen fTerGen;

  value_type operator()(const v::IVec<3> &v) { return {fTerGen(v) > 0}; }
};

} // namespace roller

#endif // ROLLER_BOOL_TERRAIN_GENERATOR_HPP_

