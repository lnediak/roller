#ifndef ROLLER_PHYSICS_CCD_NARROW_HPP_
#define ROLLER_PHYSICS_CCD_NARROW_HPP_

#include "obb.hpp"
#include "util.hpp"

namespace roller {

struct Contact {
  double t;     /// time of intersection
  v::DVec<3> p; /// contact point
  v::DVec<3> n; /// contact normal (pointing towards first object)
};

/// can intersect many pairs, but orientation must remain constant
struct CCDOBBIntersector {

  OBB p, q;

  CCDOBBIntersector(const OBB &p, const OBB &q) { initInt(p, q); }

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
  static double dmin(double a, double b) { return a > b ? b : a; }
  static double dmax(double a, double b) { return a < b ? b : a; }

  /**
    Returns [t0,t1] where t0 is intersection time (clamped to 0) and t1 is
    exit time. If no intersection with t>=0, returns an interval where t0>t1.
    Guaranteed is no intersection from 0<=t<=t0, between [a0,a1] and
    [b0+t*relVelo,b1+t*relVelo], and df1=a0-b1,df2=a1-b0
  */
  static v::DVec<2> intervalInts(double relVelo, double df1, double df2) {
    double min = 1, max = 0;
    // note that df2-df1=a1-a0+b1-b0>0 so df2>df1
    if (df1 <= 0 && df2 >= 0) {
      min = 0;
      if (relVelo > 1e-8) {
        max = df2 / relVelo;
      } else if (relVelo < -1e-8) {
        max = df1 / relVelo;
      } else {
        max = 100; // I consider an error of ~1e-6 to be acceptable
      }
    } else {
      if (df1 > 0 && relVelo > 1e-8) {
        min = df1 / relVelo;
        max = df2 / relVelo;
      } else if (df2 < 0 && relVelo < -1e-8) {
        min = df2 / relVelo;
        max = df1 / relVelo;
      } else {
        min = 1e9;
        max = 0;
      }
    }
    return {min, max};
  }

private:
  // just a helper for getInts
  struct AxisIDetail {
    double tmin, tmax;
    v::DVec<3> n;
    int type, index;
  };
  static bool updateAID(const AxisIDetail &tmp, AxisIDetail &main) {
    if (tmp.tmin > tmp.tmax) {
      main = tmp;
      return true;
    }
    if (main.tmin < tmp.tmin) {
      if (main.tmax > tmp.tmax) {
        main = tmp;
      } else {
        main.tmin = tmp.tmin;
        main.n = tmp.n;
        main.type = tmp.type;
        main.index = tmp.index;
      }
    } else if (main.tmin == tmp.tmin) {
      if (main.tmax > tmp.tmax) {
        main = tmp;
      }
    } else if (main.tmax > tmp.tmax) {
      main.tmax = tmp.tmax;
    }
    return false;
  }

public:
  // TODO: OPTIMIZE THE BELOW BY CHANGING BASES
  /// we suppose that p moves linearly by pt, and q by qt
  Contact getInts(const v::DVec<3> &pt, const v::DVec<3> &qt) {
    v::DVec<3> pA[] = {p.x, p.y, p.z}; // p vert-face
    v::DVec<3> qA[] = {q.x, q.y, q.z}; // q vert-face
    v::DVec<3> eA[] = {ee(p.x, q.x), ee(p.x, q.y), ee(p.x, q.z),
                       ee(p.y, q.x), ee(p.y, q.y), ee(p.y, q.z),
                       ee(p.z, q.x), ee(p.z, q.y), ee(p.z, q.z)};
    v::DVec<3> relVelo = qt - pt;

    Contact ret;
    AxisIDetail main{-1, 1e9, {0, 0, 0}, 0, 0};
    for (int i = 0; i < 3; i++) {
      v::DVec<2> a = p.extrema(pA[i]);
      v::DVec<2> b = q.extrema(pA[i]);
      double tmpRel = v::dot(pA[i], relVelo);
      v::DVec<2> mm = intervalInts(tmpRel, a[0] - b[1], a[1] - b[0]);
      // we don't use the index anyway in this case
      if (updateAID(
              {mm[0], mm[1], tmpRel > 0 ? pA[i] : v::DVec<3>(-pA[i]), 0, 0},
              main)) {
        break;
      }
      a = p.extrema(qA[i]);
      b = q.extrema(qA[i]);
      tmpRel = v::dot(qA[i], relVelo);
      mm = intervalInts(tmpRel, a[0] - b[1], a[1] - b[0]);
      if (updateAID(
              {mm[0], mm[1], tmpRel > 0 ? pA[i] : v::DVec<3>(-pA[i]), 1, 0},
              main)) {
        break;
      }
    }
    if (main.tmin > main.tmax) {
      ret.t = main.tmin;
      ret.n = main.n;
      return ret;
    }
    for (int i = 0; i < 9; i++) {
      v::DVec<2> a = p.extrema(eA[i]);
      v::DVec<2> b = q.extrema(eA[i]);
      double tmpRel = v::dot(eA[i], relVelo);
      v::DVec<2> mm = intervalInts(tmpRel, a[0] - b[1], a[1] - b[0]);
      if (updateAID(
              {mm[0], mm[1], tmpRel > 0 ? eA[i] : v::DVec<3>(-eA[i]), 2, i},
              main)) {
        break;
      }
    }
    ret.t = main.tmin;
    ret.n = main.n;
    if (main.tmin > main.tmax) {
      return ret;
    }

    OBB pp = p;
    pp.b += ret.t * pt;
    OBB qq = q;
    qq.b += ret.t * qt;
    switch (main.type) {
    case 0:
      ret.p = qq.maximize(ret.n);
      return ret;
    case 1:
      ret.p = pp.maximize(-ret.n);
      return ret;
    default:
      int qi = main.index % 3;
      v::DVec<3> ql = q.maxEdge(ret.n, qi);
      ret.p = p.getAlongSeg(ql, ql + q.s[qi] * qA[qi]);
      return ret;
    }
  }
  //
}; // namespace roller

} // namespace roller

#endif // ROLLER_PHYSICS_CCD_NARROW_HPP_

