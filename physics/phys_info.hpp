#ifndef ROLLER_PHYSICS_PHYS_INFO_HPP_
#define ROLLER_PHYSICS_PHYS_INFO_HPP_

#include <unordered_set>

#include "pose.hpp"

namespace roller {

struct AuxPhysInfo {
  v::DVec<3> velo; /// velocity

  DMat3x3 rot;      /// rotation matrix
  DMat3x3 riner;    /// I^{-1}(t)
  v::DVec<3> omega; /// angular velocity
};

struct PhysInfo {
  Pose pose;
  v::DVec<3> lm, am; /// linear momentum, angular momentum

  double massi;  /// 1 / mass
  DMat3x3 ineri; /// I^{-1}

  AuxPhysInfo getAuxInfo() const {
    AuxPhysInfo ret;
    ret.velo = massi * lm;
    ret.rot = pose.toRotationMatrix();
    ret.riner = ret.rot * ineri * ret.rot.transpose();
    ret.omega = ret.riner * am;
    return ret;
  }
};

struct AABB {
  v::DVec<3> m[2]; /// min, max

  bool intersects(const AABB &c) const {
    return (m[0][0] < c.m[1][0]) && (m[1][0] > c.m[0][0]) &&
           (m[0][1] < c.m[1][1]) && (m[1][1] > c.m[0][1]) &&
           (m[0][2] < c.m[1][2]) && (m[1][2] > c.m[0][2]);
  }
};

struct OBB {
  v::DVec<3> b, x, y, z; /// basically b+gx+hy+iz where 0<=g,h,i<=1
  v::DVec<3> a, c;

  OBB(const v::DVec<3> &b, const v::DVec<3> &x, const v::DVec<3> &y,
      const v::DVec<3> &z)
      : b(b), x(x), y(y), z(z), a{v::dot(b, x), v::dot(b, y), v::dot(b, z)},
        c(a + 1) {}

  AABB getAABB() const {
    // could use extrema below, but unnecessary multiplications
    AABB ret{{b, b}};
    if (x[0] > 0) {
      ret.m[1][0] += x[0];
    } else {
      ret.m[0][0] += x[0];
    }
    if (y[0] > 0) {
      ret.m[1][0] += y[0];
    } else {
      ret.m[0][0] += y[0];
    }
    if (z[0] > 0) {
      ret.m[1][0] += z[0];
    } else {
      ret.m[0][0] += z[0];
    }
    return ret;
  }

  /// minimum and maximum of u dot v, where v is a point in the OBB
  v::DVec<2> extrema(const v::DVec<3> &u) const {
    double s = v::dot(u, b);
    double a = 0, b = 0;
    double t;
    if ((t = v::dot(u, x)) > 0) {
      b += t;
    } else {
      a += t;
    }
    if ((t = v::dot(u, y)) > 0) {
      b += t;
    } else {
      a += t;
    }
    if ((t = v::dot(u, z)) > 0) {
      b += t;
    } else {
      a += t;
    }
    return {a, b};
  }

  bool intersects0(const OBB &o) const {
    v::DVec<2> tx = o.extrema(x);
    if (tx[0] > c[0] || tx[1] < a[0]) {
      return false;
    }
    v::DVec<2> ty = o.extrema(y);
    if (ty[0] > c[1] || ty[1] < a[1]) {
      return false;
    }
    v::DVec<2> tz = o.extrema(z);
    if (tz[0] > c[2] || tz[1] < a[2]) {
      return false;
    }
    return true;
  }
  bool intersects(const OBB &o) const {
    return intersects0(o) || o.intersects0(*this);
  }
};

template <class W> class World {

  struct IP {
    unsigned i, j; /// note: unordered
  };
  struct IPHash {
    std::size_t operator()(const IP &a) const noexcept {
      unsigned mix = ((a.i + a.j) << 16) + a.i * a.j;
      // murmurhash's mixing
      mix *= 0xcc9e2d51;
      mix = (mix << 15) | (mix >> 17);
      mix *= 0x1b873593;
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
  std::size_t numObjs;
  std::vector<int> xsort, ysort, zsort; /// each is 2 * obj_ind + (is AABB max)
  std::unordered_set<IP, IPHash, IPEq> aabbCol; /// collisions of AABBs

  void initAABBStuff() {
    numObjs = w.numObjects();
    xsort.resize(numObjs * 2);
    ysort.resize(numObjs * 2);
    zsort.resize(numObjs * 2);
    cached.resize(numObjs);
    std::vector<double> cached[3](numObjs * 2);
    std::size_t i = 0;
    for (auto &obj : w) {
      AABB aabb = obj.getAABB();
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
      i++;
    }
    int xyz;
    auto cmp = [this, &cached, &xyz](int a, int b) -> bool {
      return cached[xyz][a] < cached[xyz][b];
    };
    // sort...
    xyz = 0;
    std::sort(xsort.beg(), xsort.end(), cmp);
    xyz = 1;
    std::sort(ysort.beg(), ysort.end(), cmp);
    xyz = 2;
    std::sort(zsort.beg(), zsort.end(), cmp);
    // and sweep! (not most efficient lol)
    std::size_t bloat = 3 * numObjs / 2;
    decltype(aabbCol) tmpx(bloat), tmpy(bloat);
    aabbCol.reserve(bloat);
    std::unordered_set<int> active(bloat);
    for (int i : xsort) {
      if (i & 1) {
        active.erase(i / 2); // note: assuming no degenerate AABBs
      } else {
        for (int j : active) {
          tmpx.insert({i, j});
        }
        active.insert(i / 2);
      }
    }
    // active should be empty right now
    for (int i : ysort) {
      if (i & 1) {
        active.erase(i / 2);
      } else {
        int ii = i / 2;
        for (int j : active) {
          if (tmpx.count({i, j})) {
            tmpy.insert({i, j});
          }
        }
        active.insert(i / 2);
      }
    }
    for (int i : zsort) {
      if (i & 1) {
        active.erase(i / 2);
      } else {
        int ii = i / 2;
        for (int j : active) {
          if (tmpy.count({i, j})) {
            aabbCol.insert({i, j});
          }
        }
        active.insert(i / 2);
      }
    }
  }

  World(W &&w) : w(w) { initAABBStuff(); }

  // TODO: PROOFREAD PLZ
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
            aabbCols.insert({tpi / 2, tci / 2});
          } else {
            aabbCols.erase({tpi / 2, tci / 2});
          }
        }
        dsort[cind] = tpi;
        dsort[--cind] = tci;

        if (!(--cind)) {
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
  void updateAABBStuff() {
    if (numObjs != w.numObjs()) {
      // maybe can do an incremental update that's faster
      initAABBStuff();
      return;
    }
    tmp.resize(numObjs);
    for (std::size_t i = 0; i < numObjs; i++) {
      tmp[i] = w.getObj(i).getAABB();
    }
    updateAABBStuff0(xsort, 0);
    updateAABBStuff0(ysort, 1);
    updateAABBStuff0(zsort, 2);
  }
};

} // namespace roller

#endif // ROLLER_PHYSICS_PHYS_INFO_HPP_

