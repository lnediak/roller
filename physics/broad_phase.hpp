#ifndef ROLLER_PHYSICS_BROAD_PHASE_HPP_
#define ROLLER_PHYSICS_BROAD_PHASE_HPP_

#include <algorithm>
#include <unordered_set>
#include <vector>

#include "aabb.hpp"

namespace roller {

/**
  W template: All this needs is two functions:
  std::size_t numObjs() const, and
  AABB &getAABB(std::size_t).
*/
template <class W> class BroadPhaseAABB {

  struct IP {
    int i, j; /// note: unordered
  };
  struct IPHash {
    std::size_t operator()(const IP &a) const noexcept {
      unsigned mix = ((a.i + a.j) << 16) + a.i * a.j;
      // part of murmurhash's mixing
      mix *= 0xcc9e2d51;
      // mix = (mix << 15) | (mix >> 17);
      // mix *= 0x1b873593;
      return mix;
    }
  };
  struct IPEq {
    bool operator()(const IP &a, const IP &b) const noexcept {
      return (a.i == b.i && a.j == b.j) || (a.i == b.j && a.j == b.i);
    }
  };

  W w; /// access point for all info about world

  // AABB stuff
  std::size_t nObjs;
  std::vector<int> xsort, ysort, zsort; /// each is 2 * obj_ind + (is AABB max)
  std::unordered_set<IP, IPHash, IPEq> aabbInts; /// pairs of intersecting AABBs

  void initAABBStuff() {
    nObjs = w.numObjs();
    xsort.resize(nObjs * 2);
    ysort.resize(nObjs * 2);
    zsort.resize(nObjs * 2);
    std::vector<double> cached[3];
    cached[0].resize(nObjs * 2);
    cached[1].resize(nObjs * 2);
    cached[2].resize(nObjs * 2);
    for (std::size_t in = 0; in < nObjs; in++) {
      AABB aabb = w.getAABB(in);
      std::size_t i = 2 * in;
      xsort[i] = i;
      cached[0][i] = aabb.m[0][0];
      ysort[i] = i;
      cached[1][i] = aabb.m[0][1];
      zsort[i] = i;
      cached[2][i] = aabb.m[0][2];
      i++;
      xsort[i] = i;
      cached[0][i] = aabb.m[1][0];
      ysort[i] = i;
      cached[1][i] = aabb.m[1][1];
      zsort[i] = i;
      cached[2][i] = aabb.m[1][2];
    }
    int xyz;
    auto cmp = [this, &cached, &xyz](int a, int b) -> bool {
      return cached[xyz][a] < cached[xyz][b];
    };
    // sort...
    xyz = 0;
    std::sort(xsort.begin(), xsort.end(), cmp);
    xyz = 1;
    std::sort(ysort.begin(), ysort.end(), cmp);
    xyz = 2;
    std::sort(zsort.begin(), zsort.end(), cmp);
    // and sweep! (not most efficient lol)
    std::size_t bloat = 3 * nObjs / 2;
    decltype(aabbInts) tmpx(bloat), tmpy(bloat);
    aabbInts.reserve(bloat);
    std::unordered_set<int> active(bloat);
    for (int i : xsort) {
      if (i & 1) {
        active.erase(i / 2); // note: assuming no degenerate AABBs
      } else {
        int ii = i / 2;
        for (int j : active) {
          tmpx.insert({ii, j});
        }
        active.insert(ii);
      }
    }
    // active should be empty right now
    for (int i : ysort) {
      if (i & 1) {
        active.erase(i / 2);
      } else {
        int ii = i / 2;
        for (int j : active) {
          if (tmpx.count({ii, j})) {
            tmpy.insert({ii, j});
          }
        }
        active.insert(ii);
      }
    }
    for (int i : zsort) {
      if (i & 1) {
        active.erase(i / 2);
      } else {
        int ii = i / 2;
        for (int j : active) {
          if (tmpy.count({ii, j})) {
            aabbInts.insert({ii, j});
          }
        }
        active.insert(ii);
      }
    }
  }

public:
  BroadPhase(W &&w) : w(w) { initAABBStuff(); }

private:
  std::vector<AABB> tmp;
  void updateAABBStuff0(std::vector<int> &dsort, std::size_t d) {
    int rpi = dsort[0]; // real previous index
    double rp = tmp[rpi / 2].m[rpi & 1][d];
    for (std::size_t ind = 1, sz = dsort.size(); ind < sz; ind++) {
      std::size_t cind = ind;
      int tci = dsort[cind];
      int tpi = rpi;
      double tc = tmp[tci / 2].m[tci & 1][d];
      double tp = rp;
      if (tp <= tc) {
        rpi = tci;
        rp = tc;
        continue;
      }
      int rci = tci;
      int rc = tc;
      do {
        if ((tci & 1) != (tpi & 1)) {
          if (tmp[tpi / 2].intersects(tmp[tci / 2])) {
            aabbInts.insert({tpi / 2, tci / 2});
          } else {
            aabbInts.erase({tpi / 2, tci / 2});
          }
        }
        dsort[cind] = tpi;
        dsort[--cind] = tci;

        if (!cind) {
          break;
        }
        tci = tpi;
        tpi = dsort[cind - 1];
        tc = tp;
        tp = tmp[tpi / 2].m[tpi & 1][d];
        if (tp <= tc) {
          break;
        }
      } while (true);
      rpi = rci;
      rp = rc;
    }
  }

public:
  void updateAABBStuff() {
    if (nObjs != w.numObjs()) {
      // maybe can do an incremental update that's faster
      initAABBStuff();
      return;
    }
    tmp.resize(nObjs);
    for (std::size_t i = 0; i < nObjs; i++) {
      tmp[i] = w.getAABB(i);
    }
    updateAABBStuff0(xsort, 0);
    updateAABBStuff0(ysort, 1);
    updateAABBStuff0(zsort, 2);
  }

  /// calls fun(i, j) for all i,j indices of intersecting object AABBs
  template <class Fun> void exportInts(Fun &&fun) {
    updateAABBStuff();
    for (IP ip : aabbInts) {
      fun(ip.i, ip.j);
    }
  }
};

} // namespace roller

#endif // ROLLER_PHYSICS_BROAD_PHASE_HPP_

