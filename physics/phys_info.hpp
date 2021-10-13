#ifndef ROLLER_PHYSICS_PHYS_INFO_HPP_
#define ROLLER_PHYSICS_PHYS_INFO_HPP_

#include <algorithm>
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
  double el = 1;                       /// coefficient of restitution
  double sf = 0.9, kf = 0.3;           /// friction

  double massi = 0;                               /// 1 / mass
  DMat3x3 ineri{{0, 0, 0}, {0, 0, 0}, {0, 0, 0}}; /// I^{-1}

  PhysInfo() {}
  PhysInfo(double massi, const DMat3x3 &ineri, double sf, double kf)
      : sf(sf), kf(kf), massi(massi), ineri(ineri) {}

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
    return aux.velo + cross3(aux.omega, p - pose.p);
  }

  /// Euler's method, btw
  void stepTime(const AuxPhysInfo &aux, double dt, double g, bool updateVelo) {
    if (massi < 1e-12) {
      return;
    }
    v::DVec<3> nvelo = aux.velo;
    nvelo[2] -= dt * g;
    if (updateVelo) {
      lm[2] -= dt * g / massi;
    }
    pose.p += dt * nvelo;
    // XXX: MORE INTELLIGENT QUATERNION DRIFT HANDLING
    pose.q +=
        dt * 0.5 *
        quaternionMult(v::DVec<4>{0, aux.omega[0], aux.omega[1], aux.omega[2]},
                       pose.q);
    pose.q /= std::sqrt(v::norm2(pose.q));
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

struct Contact {
  double dist;
  v::DVec<3> p, n;

  Contact lo(const Contact &c) const {
    if (dist >= 0) {
      return dist < c.dist ? *this : c;
    }
    if (c.dist >= 0) {
      return *this;
    }
    return dist < c.dist ? c : *this;
  }
  Contact hi(const Contact &c) const { return dist < c.dist ? *this : c; }
};

struct OBB {
  v::DVec<3> b, s, x, y, z; /// basically b+gx+hy+iz where 0<=g,h,i<=s
  v::DVec<3> a, c;

  OBB() {}
  OBB(const v::DVec<3> &b, const v::DVec<3> &s, const v::DVec<3> &x,
      const v::DVec<3> &y, const v::DVec<3> &z)
      : b(b), s(s), x(x), y(y),
        z(z), a{v::dot(b, x), v::dot(b, y), v::dot(b, z)}, c(a + s) {}

  bool isIn(const v::DVec<3> &p) const {
    v::DVec<3> pd = {v::dot(p, x), v::dot(p, y), v::dot(p, z)};
    return a[0] <= pd[0] && a[1] <= pd[1] && a[2] <= pd[2] && pd[0] <= c[0] &&
           pd[1] <= c[1] && pd[2] <= c[2];
  }

  /// XXX: TRIGGER WARNING: NOT OPTIMIZED
  AABB getAABB() const {
    v::DVec<3> c0 = b;
    v::DVec<3> c1 = b + s[0] * x;
    v::DVec<3> c2 = b + s[1] * y;
    v::DVec<3> c3 = b + s[2] * z;
    v::DVec<3> c4 = c1 + s[1] * y;
    v::DVec<3> c5 = c2 + s[2] * z;
    v::DVec<3> c6 = c3 + s[0] * x;
    v::DVec<3> c7 = c4 + s[2] * z;
    return {
        elementwiseMin(
            elementwiseMin(elementwiseMin(c0, c1), elementwiseMin(c2, c3)),
            elementwiseMin(elementwiseMin(c4, c5), elementwiseMin(c6, c7))),
        elementwiseMax(
            elementwiseMax(elementwiseMax(c0, c1), elementwiseMax(c2, c3)),
            elementwiseMax(elementwiseMax(c4, c5), elementwiseMax(c6, c7)))};
  }

  /// minimum and maximum of u dot v, where v is a point in the OBB
  v::DVec<2> extrema(const v::DVec<3> &u) const {
    double o = v::dot(u, b);
    double lo = o, hi = o;
    double t;
    if ((t = v::dot(u, x)) > 0) {
      hi += s[0] * t;
    } else {
      lo += s[0] * t;
    }
    if ((t = v::dot(u, y)) > 0) {
      hi += s[1] * t;
    } else {
      lo += s[1] * t;
    }
    if ((t = v::dot(u, z)) > 0) {
      hi += s[2] * t;
    } else {
      lo += s[2] * t;
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
  // I'll implement this efficiently later. Right now, I'm jjust
  // a little sick of this and want it to WORKKKKKKK already
  // XXX: Use GJK/EPA or something instead of my hideous garbage below
  Contact deepest(const v::DVec<3> &u, const OBB &o, bool flipN) const {
    Contact ret;
    ret.dist = 100;
    for (int d1 = 0; d1 < 3; d1++) {
      int d2 = d1 >= 2 ? 0 : d1 + 1;
      int d3 = d2 >= 2 ? 0 : d2 + 1;
      for (int m1 = 0; m1 < 2; m1++) {
        for (int m2 = 0; m2 < 2; m2++) {
          v::DVec<3> pp =
              o.b + m1 * o.s[d2] * o.dirD(d2) + m2 * o.s[d3] * o.dirD(d3);
          Contact tmp = deepestAlongSeg(u, pp, pp + o.s[d1] * o.dirD(d1));
          ret = ret.hi(tmp);
          std::cout << "deepest. tmp: " << tmp.dist << std::endl;
          if (tmp.dist < 0) {
            std::cout << tmp.p << tmp.n << std::endl;
          }
        }
      }
    }
    if (flipN) {
      ret.n = -ret.n;
    }
    std::cout << "returning: " << ret.dist << std::endl;
    if (ret.dist < 0) {
      std::cout << ret.p << ret.n << std::endl;
    }
    std::cout << std::endl;
    return ret;
  }
  v::DVec<3> dirD(int i) const {
    switch (i) {
    case 0:
      return x;
    case 1:
      return y;
    case 2:
      return z;
    }
    return {0, 0, 0};
  }
  // normal will point away from *this
  Contact deepestAlongSeg(const v::DVec<3> &u, v::DVec<3> p,
                          v::DVec<3> q) const {
    std::cout << "deepestAlongSeg: " << u << p << q << std::endl;
    if (v::dot(p, u) < v::dot(q, u)) {
      v::DVec<3> tmp = p;
      p = q;
      q = tmp;
    }
    v::DVec<3> d = q - p;
    v::DVec<3> dd = {v::dot(d, x), v::dot(d, y), v::dot(d, z)};
    v::DVec<3> pp = {v::dot(p, x), v::dot(p, y), v::dot(p, z)};
    std::cout << "a, c, dd, pp: " << a << c << dd << pp << std::endl;
    double lo = 0, hi = 1;
    for (int ind = 0; ind < 3; ind++) {
      if (-1e-3 < dd[ind] && dd[ind] < 1e-3) {
        if ((a[ind] < pp[ind]) ^ (pp[ind] < c[ind])) {
          Contact ret;
          ret.dist = 100;
          return ret;
        }
      }
      double tlo = (a[ind] - pp[ind]) / dd[ind];
      double thi = (c[ind] - pp[ind]) / dd[ind];
      if (tlo > thi) {
        lo = std::fmax(lo, thi);
        hi = std::fmin(hi, tlo);
      } else {
        lo = std::fmax(lo, tlo);
        hi = std::fmin(hi, thi);
      }
    }
    std::cout << "lo, hi: " << lo << " " << hi << std::endl;
    if (hi >= 1 || hi >= 1 - lo) {
      return {(1 - lo) * v::dot(d, u), p + hi * q, u};
    }
    return {hi * v::dot(d, u), p + lo * d, -u};
  }
};

/// XXX: TRIGGER WARNING: NOT OPTIMIZED
struct OBBIntersector {

  OBB p, q;

  OBBIntersector(const OBB &p, const OBB &q) { initInt(p, q); }

  /// XXX: evaluate orientation-related constants for optimization
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
    return cros / std::sqrt(nrm);
  }
  /// assuming u is unit vector
  bool vc(const v::DVec<3> &u) {
    v::DVec<2> a = p.extrema(u);
    v::DVec<2> b = q.extrema(u);
    return a[0] > b[1] || a[1] < b[0];
  }
  Contact getInts() {
    if (vc(p.x) || vc(p.y) || vc(p.z) || vc(q.x) || vc(q.y) || vc(q.z) ||
        vc(ee(p.x, q.x)) || vc(ee(p.x, q.y)) || vc(ee(p.x, q.z)) ||
        vc(ee(p.y, q.x)) || vc(ee(p.y, q.y)) || vc(ee(p.y, q.z)) ||
        vc(ee(p.z, q.x)) || vc(ee(p.z, q.y)) || vc(ee(p.z, q.z))) {
      Contact ret;
      ret.dist = 10; // junk
      return ret;
    }
    // I'm not even trying to make this even remotely optimized
    Contact ret;
    ret.dist = 100;
    ret = ret.lo(p.deepest(q.x, q, true));
    ret = ret.lo(p.deepest(q.y, q, true));
    ret = ret.lo(p.deepest(q.z, q, true));
    ret = ret.lo(q.deepest(p.x, p, false));
    ret = ret.lo(q.deepest(p.y, p, false));
    ret = ret.lo(q.deepest(p.z, p, false));

    std::cout << "after vertices. ret: " << ret.dist << std::endl;
    if (ret.dist < 0) {
      std::cout << ret.p << ret.n << std::endl;
    }

    ret = ret.lo(q.deepest(ee(p.x, q.x), p, false));
    ret = ret.lo(q.deepest(ee(p.x, q.y), p, false));
    ret = ret.lo(q.deepest(ee(p.x, q.z), p, false));
    ret = ret.lo(q.deepest(ee(p.y, q.x), p, false));
    ret = ret.lo(q.deepest(ee(p.y, q.y), p, false));
    ret = ret.lo(q.deepest(ee(p.y, q.z), p, false));
    ret = ret.lo(q.deepest(ee(p.z, q.x), p, false));
    ret = ret.lo(q.deepest(ee(p.z, q.y), p, false));
    ret = ret.lo(q.deepest(ee(p.z, q.z), p, false));

    std::cout << "OBBIntersector.getInts: " << ret.dist << " " << ret.p << " "
              << ret.n << std::endl;
    std::cout << "with " << p.b << p.s << p.x << p.y << p.z << std::endl;
    std::cout << "and " << q.b << q.s << q.x << q.y << q.z << std::endl;
    return ret;
  }
};

/// basically handles broad-phase
template <class W> class World {

  struct IP {
    int i, j; /// note: unordered
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
      AABB aabb = w.getObj(in).getAABB();
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
  World(W &&w) : w(w) { initAABBStuff(); }

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
            std::cout << "interesting..." << std::endl;
            std::cout << tmp[tpi / 2].m[0] << " " << tmp[tpi / 2].m[1]
                      << std::endl;
            std::cout << tmp[tci / 2].m[0] << " " << tmp[tci / 2].m[1]
                      << std::endl;
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
      tmp[i] = w.getObj(i).getAABB();
    }
    updateAABBStuff0(xsort, 0);
    updateAABBStuff0(ysort, 1);
    updateAABBStuff0(zsort, 2);
  }

  std::size_t numObjs() const { return nObjs; }
  decltype(w.getObj(0)) &getObj(std::size_t i) { return w.getObj(i); }

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

