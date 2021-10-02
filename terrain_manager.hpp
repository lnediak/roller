#ifndef ROLLER_TERRAIN_MANAGER_HPP_
#define ROLLER_TERRAIN_MANAGER_HPP_

#include <memory>
#include <unordered_map>
#include <vector>

#include "hashtable_cache.hpp"
#include "util.hpp"
#include "vector.hpp"

namespace roller {

/**
  TerGen must be a functor that takes v::IVec<3> as input, and produces
  TerGen::value_type, which must be a valid blockdata type, which means that it
  must contain functions isOpaque(), isAir(), and getColor(int), which
  takes a value 0..7 where the bits represent which corner of the cube we're
  getting the color of, and coordinate 0 is the least significant digit while
  coordinate 2 is the most significant digit.
*/
template <class TerGen> struct TerrainManager {

  typedef typename TerGen::value_type blockdata;

  TerGen terGen;

  static const std::size_t CHNK_BT = 4;
  static const std::size_t CHNK_SL = 16;
  static const std::size_t CHNK_VL = 4096;
  class Chunk {
    v::IVec<3> coord;  /// in chunk coordinates
    v::IVec<3> scoord; /// in world coordinates (divisible by CHNK_SL)
    std::unique_ptr<blockdata[]> data;

    /// format: x, y, z, c; stride is 4 floats (c is uint)
    /// will backface cull, so all will be ccw.
    std::vector<float> triangles;
    bool meshed = false;

  public:
    Chunk(const v::IVec<3> &coord, TerGen &terGen)
        : coord(coord), scoord(CHNK_SL * coord), data(new blockdata[CHNK_VL]) {
      simpleGen(terGen);
    }
    Chunk() {}

    bool isInChunk(const v::IVec<3> &v) const { return v / CHNK_SL == coord; }

    v::IVec<3> chunkLocal(const v::IVec<3> &v) const { return v - scoord; }
    std::size_t getIndL(const v::IVec<3> &v) const {
      return v[2] + CHNK_SL * v[1] + CHNK_SL * CHNK_SL * v[0];
    }
    blockdata getValueL(const v::IVec<3> &v) const { return data[getIndL(v)]; }
    std::size_t getInd(const v::IVec<3> &v) const {
      return getIndL(chunkLocal(v));
    }
    blockdata getValue(const v::IVec<3> &v) const { return data[getInd(v)]; }

    void simpleGen(TerGen &terGen) {
      v::IVec<3> scm = scoord + CHNK_SL;
      v::IVec<3> v;
      std::size_t ind = 0;
      for (v[0] = scoord[0]; v[0] < scm[0]; v[0]++) {
        for (v[1] = scoord[1]; v[1] < scm[1]; v[1]++) {
          for (v[2] = scoord[2]; v[2] < scm[2]; v[2]++, ind++) {
            data[ind] = terGen(v);
          }
        }
      }
    }

  private:
    void triPush(const v::IVec<3> &v, unsigned color) {
      triangles.push_back(v[0]);
      triangles.push_back(v[1]);
      triangles.push_back(v[2]);
      triangles.push_back(*((float *)&color));
    }

  public:
    /**
      dim can be 1, 2, 3, -1, -2, or -3, for +x, +y, +z, -x, -y, -z,
      respectively
    */
    void addTriangle(v::IVec<3> v, blockdata bd, int dim) {
      int d0 = dim > 0 ? dim - 1 : -dim - 1;
      int d1 = d0 >= 2 ? 0 : d0 + 1;
      int d2 = d1 >= 2 ? 0 : d1 + 1;

      int bb = 0;
      // the key here is that v_d1 X v_d2 = looking flat at triangle
      if (dim > 0) {
        int tmp = d1;
        d1 = d2;
        d2 = tmp;
        v[d0]++;
        bb += 1 << d0;
      }
      v[d1]++;
      bb += 1 << d1;
      triPush(v, bd.getColor(bb));
      v[d1]--;
      bb -= 1 << d1;
      triPush(v, bd.getColor(bb));
      v[d2]++;
      bb += 1 << d2;
      triPush(v, bd.getColor(bb));

      triPush(v, bd.getColor(bb));
      v[d1]++;
      bb += 1 << d1;
      triPush(v, bd.getColor(bb));
      v[d2]--;
      bb -= 1 << d2;
      triPush(v, bd.getColor(bb));
    }

    /**
      Only does anything if this chunk is not already meshed.
      dim can be 0, 1, or 2.
    */
    void meshAllIfNo(const Chunk &lo, const Chunk &hi, int dim) {
      if (meshed) {
        return;
      }
      triangles.clear();
      int d0 = dim;
      int d1 = d0 >= 2 ? 0 : d0 + 1;
      int d2 = d1 >= 2 ? 0 : d1 + 1;

      v::IVec<3> scm = scoord + CHNK_SL;
      v::IVec<3> v;
      for (v[d1] = scoord[d1]; v[d1] < scm[d1]; v[d1]++) {
        for (v[d2] = scoord[d2]; v[d2] < scm[d2]; v[d2]++) {
          v[d0] = scoord[d0] - 1;
          blockdata lobd = lo->getValue(v);
          v[d0]++;
          blockdata hibd = getValue(v);
          if (lobd.isAir() && !hibd.isAir()) {
            addTriangle(v, -1 - d0);
          }
          lobd = hibd;
          v[d0]++;
          for (; v[d0] < scm[d0]; v[d0]++) {
            hibd = getValue(v);
            if (lobd.isAir() && !hibd.isAir()) {
              addTriangle(v, -1 - d0);
            } else if (!lobd.isAir() && hibd.isAir()) {
              v[d0]--;
              addTriangle(v, d0 + 1);
              v[d0]++;
            }
            lobd = hibd;
          }
          hibd = hi->getValue(v);
          if (!lobd.isAir() && hibd.isAir()) {
            v[d0]--;
            addTriangle(v, d0 + 1);
            v[d0]++;
          }
        }
      }
      meshed = true;
    }

    template <class Fun> void exportTriangles(Fun &&fun) {
      fun(&triangles[0], triangles.size());
    }
  };

  HashtableCache<v::IVec<3>, Chunk, v::IVecHash<3>,
                 v::EqualFunctor<v::IVec<3>, v::IVec<3>>>
      ccache;

  template <class Fun> void exportAllTriangles(const SliceDirs &sd, Fun &&fun) {
    // XXX: PROPER CHUNK-IN-VIEW CALCULATIONS
    std::unordered_map<v::IVec<3>, Chunk *> chunks(4096);
    v::IVec<3> mins = sdMins(sd) / CHNK_SL - 1;
    v::IVec<3> maxs = (sdMaxs(sd) + 1) / CHNK_SL + 1;
    v::IVec<3> cc;
    for (cc[0] = mins[0]; cc[0] <= maxs[0]; cc[0]++) {
      for (cc[1] = mins[1]; cc[1] <= maxs[1]; cc[1]++) {
        for (cc[2] = mins[2]; cc[2] <= maxs[2]; cc[2]++) {
          chunks[cc] = &ccache.applyIfNew(cc, [this, &cc](Chunk &c) -> void {
            c.~Chunk();
            new (&c) Chunk(cc, terGen);
          });
        }
      }
    }
    mins += 1;
    maxs -= 1;
    for (int d0 = 0; d0 < 3; d0++) {
      int d1 = d0 >= 2 ? 0 : d0 + 1;
      int d2 = d1 >= 2 ? 0 : d1 + 1;
      for (cc[d0] = mins[d0]; cc[d0] <= maxs[d0]; cc[d0]++) {
        for (cc[d1] = mins[d1]; cc[d1] <= maxs[d1]; cc[d1]++) {
          for (cc[d2] = mins[d2]; cc[d2] <= maxs[d2]; cc[d2]++) {
            cc[d2]--;
            const Chunk &lo = chunks[cc];
            cc[d2] += 2;
            const Chunk &hi = chunks[cc];
            cc[d2]--;
            chunks[cc].meshAllIfNo(lo, hi, d2);
            chunks[cc].exportTriangles(fun);
          }
        }
      }
    }
  }
};

} // namespace roller

#endif // ROLLER_TERRAIN_MANAGER_HPP_

