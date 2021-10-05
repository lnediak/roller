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
  v::DVec<3> lm{0, 0, 0}, am{0, 0, 0}; /// linear momentum, angular momentum

  double massi = 0;                               /// 1 / mass
  DMat3x3 ineri{{0, 0, 0}, {0, 0, 0}, {0, 0, 0}}; /// I^{-1}

  PhysInfo(double massi, const DMat3x3 &ineri) : massi(massi), ineri(ineri) {}

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
  v::DVec<3> b, s, x, y, z; /// basically b+gx+hy+iz where 0<=g,h,i<=s
  v::DVec<3> a, c;

  OBB(const v::DVec<3> &b, const v::DVec<3> &s, const v::DVec<3> &x,
      const v::DVec<3> &y, const v::DVec<3> &z)
      : b(b), s(s), x(x), y(y),
        z(z), a{v::dot(b, x), v::dot(b, y), v::dot(b, z)}, c(a + s) {}

  AABB getAABB() const {
    // could use extrema below, but unnecessary multiplications
    AABB ret{{b, b}};
    if (x[0] > 0) {
      ret.m[1][0] += s[0] * x[0];
    } else {
      ret.m[0][0] += s[0] * x[0];
    }
    if (y[0] > 0) {
      ret.m[1][0] += s[1] * y[0];
    } else {
      ret.m[0][0] += s[1] * y[0];
    }
    if (z[0] > 0) {
      ret.m[1][0] += s[2] * z[0];
    } else {
      ret.m[0][0] += s[2] * z[0];
    }
    return ret;
  }
};

struct Contact {
  double dist;
  v::DVec<3> p, n;
};

struct OBBIntersector {
  v::DVec<3> wor[2][3];
  v::DVec<3> rel[2][3];

  v::DVec<3> min[2], max[2];

  OBBIntersector(const OBB &p, const OBB &q)
      : wor{{p.x, p.y, p.z}, {q.x, q.y, q.z}} {
    DMat3x3 relm = DMat3x3{p.x, p.y, p.z} * DMat3x3{q.x, q.y, q.z}.transpose();
    rel[1][0] = relm.a;
    rel[1][1] = relm.b;
    rel[1][2] = relm.c;
    relm = relm.transpose();
    rel[0][0] = relm.a;
    rel[0][1] = relm.b;
    rel[0][2] = relm.c;
  }

  void initInt(const OBB &p, const OBB &q) {
    min[0] = p.a;
    min[1] = q.a;
    max[0] = p.c;
    max[1] = q.c;
  }

  void upDown0(double &lo, double &hi, double a) const {
    if (a > 0) {
      hi += a;
    } else {
      lo += a;
    }
  }
  void upDown03(double &lo, double &hi, double a, double b, double c) const {
    upDown0(lo, hi, a);
    upDown0(lo, hi, b);
    upDown0(lo, hi, c);
  }
  void upv0(v::DVec<3> &hi, int w, int d, double a) const {
    if (a > 0) {
      hi += max[w][d] * wor[w][d];
    } else {
      hi += min[w][d] * wor[w][d];
    }
  }
  void upv03(v::DVec<3> &hi, int w, double a, double b, double c) const {
    hi = 0;
    upv0(hi, w, 0, a);
    upv0(hi, w, 1, b);
    upv0(hi, w, 2, c);
  }
  void downv0(v::DVec<3> &lo, int w, int d, double a) const {
    if (a > 0) {
      lo += min[w][d] * wor[w][d];
    } else {
      lo += max[w][d] * wor[w][d];
    }
  }
  void downv03(v::DVec<3> &lo, int w, double a, double b, double c) const {
    lo = 0;
    downv0(lo, w, 0, a);
    downv0(lo, w, 1, b);
    downv0(lo, w, 2, c);
  }

  /// w: 0 or 1; d: 0, 1, or 2
  Contact distByFace(int w, int d) const {
    double melo = m[w][d];
    double mehi = x[w][d];
    double yulo, yuhi;
    double a = rel[w][0][d], b = rel[w][1][d], c = rel[w][2][d];
    upDown03(yulo, yuhi, a, b, c);

    double forP = yuhi - melo;
    double bakP = yulo - mehi;
    if (forP < 0 && bakP < 0) {
      if (forP > bakP) {
        v::DVec<3> hi;
        upv03(hi, w, a, b, c);
        return {forP, hi, wor[w][d]};
      } else {
        v::DVec<3> lo;
        downv03(lo, w, a, b, c);
        return {bakP, lo, -wor[w][d]};
      }
    }
    Contact c;
    c.dist = forP > bakP ? bakP : forP;
    return c;
  }
  /// d0, d1: 0, 1, or 2
  Contact distByEE(int d0, int d1) const {
    //
  }
};

/// basically handles broad-phase
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
  std::unordered_set<IP, IPHash, IPEq> aabbInts; /// pairs of intersecting AABBs

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
    decltype(aabbInts) tmpx(bloat), tmpy(bloat);
    aabbInts.reserve(bloat);
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
            aabbInts.insert({i, j});
          }
        }
        active.insert(i / 2);
      }
    }
  }

  World(W &&w) : w(w) { initAABBStuff(); }

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

  /// calls fun(obj1, obj2) for all obj1 and obj2 with intersecting AABBs
  template <class Fun> void exportAllAABBInts(Fun &&fun) {
    updateAABBStuff();
    for (IP ip : aabbInts) {
      fun(w.getObj(ip.i), w.getObj(ip.j));
    }
  }
};

} // namespace roller

#endif // ROLLER_PHYSICS_PHYS_INFO_HPP_

