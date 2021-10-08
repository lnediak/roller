#ifndef ROLLER_PHYSICS_PHYS_INFO_HPP_
#define ROLLER_PHYSICS_PHYS_INFO_HPP_

#include <cmath>
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
  double el;                           /// coefficient of restitution
  double sf, kf;                       /// friction

  double massi = 0;                               /// 1 / mass
  DMat3x3 ineri{{0, 0, 0}, {0, 0, 0}, {0, 0, 0}}; /// I^{-1}

  PhysInfo(double massi, const DMat3x3 &ineri, double sf, double kf)
      : massi(massi), ineri(ineri), sf(sf), kf(kf) {}

  AuxPhysInfo getAuxInfo() const {
    AuxPhysInfo ret;
    ret.velo = massi * lm;
    ret.rot = pose.toRotationMatrix();
    ret.riner = ret.rot * ineri * ret.rot.transpose();
    ret.omega = ret.riner * am;
    return ret;
  }

  /// velocity at point p
  v::DVec<3> getVelocity(const AuxPhysInfo &aux, const v::DVec<3> &p) const {
    return velo + cross3(aux.omega, p - pose.p);
  }

  /// Euler's method, minus updating momentum
  void stepTime(const AuxPhysInfo &aux, double dt) {
    if (massi < 1e-12) {
      return;
    }
    pose.p += dt * aux.velo;
    pose.q += dt * 0.5 * aux.omega * pose.q;
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

  /// minimum and maximum of u dot v, where v is a point in the OBB
  v::DVec<2> extrema(const v::DVec<3> &u) const {
    double o = v::dot(u, b);
    double lo = o, hi = o;
    out = b;
    double t;
    if ((t = v::dot(u, x)) > 0) {
      hi += s[0] * t;
    } else {
      lo += s[0] * t;
    }
    if ((t = v::dot(u, y)) > 0) {
      hi += s[0] * t;
    } else {
      lo += s[0] * t;
    }
    if ((t = v::dot(u, z)) > 0) {
      hi += s[0] * t;
    } else {
      lo += s[0] * t;
    }
    return {lo, hi};
  }
  // finds the value of v such that u dot v is maximum (like extrema)
  v::DVec<3> maximize(const v::DVec<3> &u) const {
    v::DVec<3> ret = b;
    if (v::dot(u, x) > 0) {
      ret += s[0] * x;
    }
    if (v::dot(u, y) > 0) {
      ret += s[1] * y;
    }
    if (v::dot(u, z) > 0) {
      ret += s[2] * z;
    }
    return ret;
  }
};

struct Contact {
  double dist;
  v::DVec<3> p, n;

  Contact lo(const Contact &c) const {
    if (dist < c.dist) {
      return *this;
    }
    return c;
  }
};

/// XXX: TRIGGER WARNING: NOT OPTIMIZED
struct OBBIntersector {

  OBB p, q;

  /// XXX: evaluate orientation-related constants for optimization
  OBBIntersector(const OBB &p, const OBB &q) { initInt(p, q); }

  void initInt(const OBB &pp, const OBB &qq) {
    p = pp;
    q = qq;
  }

  v::DVec<3> ee(const v::DVec<3> &a, const v::DVec<3> &b) {
    v::DVec<3> cros = cross3(a, b);
    double nrm = v::norm2(cros);
    if (nrm < 1e-12) {
      return p.x;
    }
    return cros / std::sqrt(nrm); // I told this you this is not optimized
  }
  /// assuming u is unit vector
  Contact vc(const v::DVec<3> &u) {
    v::DVec<2> a = p.extrema(u);
    v::Dvec<2> b = q.extrema(u);
    double forD = a[0] - b[1];
    double bakD = b[0] - a[1];
    if (forD < 0 && bakD < 0) {
      if (forD > bakD) {
        return {forD, q.maximize(u), u};
      }
      return {bakD, q.maximize(-u), -u};
    }
    Contact c;
    c.dist = forD > bakD ? bakD : forD;
    return c;
  }
  Contact getInts() {
    cons.clear();
    Contact ret = vc(p.x).lo(vc(p.y).lo(vc(p.z)));
    ret = ret.lo(vc(q.x).lo(vc(q.y).lo(vc(q.z))));
    ret = ret.lo(vc(ee(p.x, q.x)).lo(vc(ee(p.x, q.y)).lo(vc(ee(p.x, q.z)))));
    ret = ret.lo(vc(ee(p.y, q.x)).lo(vc(ee(p.y, q.y)).lo(vc(ee(p.y, q.z)))));
    return ret.lo(vc(ee(p.z, q.x)).lo(vc(ee(p.z, q.y)).lo(vc(ee(p.z, q.z)))));
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
    numObjs = w.numObjs();
    xsort.resize(numObjs * 2);
    ysort.resize(numObjs * 2);
    zsort.resize(numObjs * 2);
    cached.resize(numObjs);
    std::vector<double> cached[3](numObjs * 2);
    for (std::size_t in = 0; in < numObjs; in++) {
      AABB aabb = w.getObj(i).getAABB();
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

  /// calls fun(i, j) for all i,j indices of intersecting object AABBs
  template <class Fun> void exportAllAABBInts(Fun &&fun) {
    updateAABBStuff();
    for (IP ip : aabbInts) {
      fun(ip.i, ip.j);
    }
  }
};

} // namespace roller

#endif // ROLLER_PHYSICS_PHYS_INFO_HPP_

