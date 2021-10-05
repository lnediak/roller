#ifndef ROLLER_BOOL_TERRAIN_GENERATOR_HPP_
#define ROLLER_BOOL_TERRAIN_GENERATOR_HPP_

#include "util.hpp"
#include "vector.hpp"

namespace roller {

struct BoolBlockdata {
  bool b;
  bool isAir() const { return !b; }
  bool isOpaque() const { return b; }
  Color getColor(int bb) const {
    if (!b) {
      return {{0, 0, 0, 0}};
    }
    const Color results[] = {
        {0x00, 0x00, 0x00, 0xFF}, {0x22, 0x22, 0x22, 0xFF},
        {0x44, 0x44, 0x44, 0xFF}, {0x66, 0x66, 0x66, 0xFF},
        {0x88, 0x88, 0x88, 0xFF}, {0xAA, 0xAA, 0xAA, 0xFF},
        {0xCC, 0xCC, 0xCC, 0xFF}, {0xFF, 0xFF, 0xFF, 0xFF}};
    return results[bb];
    return {0xFF, 0xFF, 0, 0};
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
    v::DVec<3> lerp;
    for (v[0] = 0; v[0] < 0.9999; v[0] += s) {
      lerp[0] = v[0] * v[0] * (3 - 2 * v[0]);
      for (v[1] = 0; v[1] < 0.9999; v[1] += s) {
        lerp[1] = v[1] * v[1] * (3 - 2 * v[1]);
        for (v[2] = 0; v[2] < 0.9999; v[2] += s, ind++) {
          lerp[2] = v[2] * v[2] * (3 - 2 * v[2]);
          data[ind] = {vlerpAll(vecpow, lerp) > 0};
        }
      }
    }
  }
};

} // namespace roller

#endif // ROLLER_BOOL_TERRAIN_GENERATOR_HPP_

