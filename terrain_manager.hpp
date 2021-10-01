#ifndef ROLLER_TERRAIN_MANAGER_HPP_
#define ROLLER_TERRAIN_MANAGER_HPP_

#include <memory>
#include <unordered_map>

namespace roller {

template <class TerGen> struct TerrainManager {

  typedef TerGen::value_type value_type;

  TerGen terGen;

  static const std::size_t CHNK_BT = 4;
  static const std::size_t CHNK_SL = 16;
  static const std::size_t CHNK_VL = 4096;
  class Chunk {
    v::IVec<3> coord;  /// in chunk coordinates
    v::IVec<3> scoord; /// in world coordinates (divisible by CHNK_SL)
    std::unique_ptr<value_type[]> data;
    std::unique_ptr<int[]> vInds; /// visible indices
    std::size_t viSz;             /// length of above

    bool marked = false;

  public:
    Chunk(const v::IVec<3> &coord)
        : coord(coord), scoord(CHNK_SL * coord), data(new value_type[CHNK_VL]) {
    }

    bool isInChunk(const v::IVec<3> &v) const { return v / CHNK_SL == coord; }

    v::IVec<3> chunkLocal(const v::IVec<3> &v) const { return v - scoord; }
    std::size_t getIndL(const v::IVec<3> &v) const {
      return v[2] + CHNK_SL * v[1] + CHNK_SL * CHNK_SL * v[0];
    }
    value_type getValueL(const v::IVec<3> &v) const { return data[getIndL(v)]; }
    std::size_t getInd(const v::IVec<3> &v) const {
      return getIndL(chunkLocal(v));
    }
    value_type getValue(const v::IVec<3> &v) const { return data[getInd(v)]; }

    void simpleGen(TerGen &terGen) {
      v::IVec<3> scm = scoord + CHNK_SL;
      v::IVec<3> v;
      v[0] = scoord[0];
      std::size_t ind = 0;
      for (; v[0] < scm[0]; v[0]++) {
        v[1] = scoord[1];
        for (; v[1] < scm[1]; v[1]++) {
          v[2] = scoord[2];
          for (; v[2] < scm[2]; v[2]++, ind++) {
            data[ind] = terGen(v);
          }
        }
      }
    }

    // TODO: IMPLEMENT WEIRD FLOOD FILL
  };

  typedef std::unordered_map<v::IVec<3>, Chunk, v::IVecHash<3>> cmap;
  cmap chunks;
};

} // namespace roller

#endif // ROLLER_TERRAIN_MANAGER_HPP_

