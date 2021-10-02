#ifndef ROLLER_BOOL_TERRAIN_GENERATOR_HPP_
#define ROLLER_BOOL_TERRAIN_GENERATOR_HPP_

#include "util.hpp"
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
  void operator()(v::IVec<3> sc, value_type *data, int sl) {
    v::DVec<8> vecpow;
    vecpow[0] = fTerGen(sc);
    sc[0] += sl;
    vecpow[1] = fTerGen(sc);
    sc[1] += sl;
    vecpow[3] = fTerGen(sc);
    sc[0] -= sl;
    vecpow[2] = fTerGen(sc);
    sc[2] += sl;
    vecpow[6] = fTerGen(sc);
    sc[0] += sl;
    vecpow[7] = fTerGen(sc);
    sc[1] -= sl;
    vecpow[5] = fTerGen(sc);
    sc[0] -= sl;
    vecpow[4] = fTerGen(sc);
    sc[2] -= sl;

    double s = 1. / sl;
    v::DVec<3> v;
    std::size_t ind = 0;
    for (v[0] = 0; v[0] < 0.9999; v[0] += s) {
      for (v[1] = 0; v[1] < 0.9999; v[1] += s) {
        for (v[2] = 0; v[2] < 0.9999; v[2] += s, ind++) {
          v::DVec<3> lerp = v * v * (3 - 2 * v);
          data[ind] = {vlerpAll(vecpow, lerp) > 0};
        }
      }
    }
  }
};

} // namespace roller

#endif // ROLLER_BOOL_TERRAIN_GENERATOR_HPP_

