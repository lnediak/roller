#ifndef ROLLER_PHYSICS_CCD_NARROW_HPP_
#define ROLLER_PHYSICS_CCD_NARROW_HPP_

#include "obb.hpp"
#include "pose.hpp"
#include "util.hpp"

namespace roller {

struct Contact {
  double t;     /// time of intersection
  double d;     /// depth of penetration
  v::DVec<3> p; /// contact point
  v::DVec<3> n; /// contact normal (pointing towards first object)
};

/// can intersect many pairs, but orientation must remain constant
struct CCDOBBIntersector {

  OBB p, q;

  CCDOBBIntersector() {}
  CCDOBBIntersector(const OBB &p, const OBB &q) { initInt(p, q); }

  /// XXX: evaluate orientation-related constants for optimization
  void initInt(const OBB &pp, const OBB &qq) {
    p = pp;
    q = qq;
  }

  /// only sets position and sidelength
  void setOBBs(OBB pp, OBB qq) {
    p.b = pp.b;
    p.s = pp.s;
    q.b = qq.b;
    q.s = qq.s;
    p.properify();
    q.properify();
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
    Returns [t0,t1,d] where t0 is intersection time (clamped to 0), t1 is exit
    time, and d is penetration depth. If no intersection with t>=0, returns an
    interval where t0>t1. Guaranteed is no intersection from 0<=t<t0, between
    [a0,a1] and [b0+t*relVelo,b1+t*relVelo], and df1=a0-b1,df2=a1-b0
  */
  static v::DVec<3> intervalInt(double relVelo, double df1, double df2) {
    double min = 1, max = 0, depth = 0;
    // note that df2-df1=a1-a0+b1-b0>0 so df2>df1
    if (df1 <= 0 && df2 >= 0) {
      min = 0;
      depth = dmin(-df1, df2);
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
    return {min, max, depth};
  }

private:
  // just a helper for getInt
  struct AxisIDetail {
    double tmin, tmax, depth; /// note: tmin and depth cannot both be nonzero
    v::DVec<3> n;
    int type, index;
  };
  static bool updateAID(const AxisIDetail &tmp, AxisIDetail &main) {
    std::cout << "tmp.tmin,tmp.tmax,tmp.depth,tmp.n: " << tmp.tmin << " "
              << tmp.tmax << " " << tmp.depth << " " << tmp.n << std::endl;
    std::cout << "main.tmin,main.tmax,main.depth,main.n: " << main.tmin << " "
              << main.tmax << " " << main.depth << " " << main.n << std::endl;
    if (tmp.tmin > tmp.tmax) {
      main = tmp;
      return true;
    }
    if (main.tmin < tmp.tmin) {
      if (main.tmax > tmp.tmax) {
        main = tmp;
      } else {
        main.tmin = tmp.tmin;
        main.depth = tmp.depth;
        main.n = tmp.n;
        main.type = tmp.type;
        main.index = tmp.index;
        if (main.tmin > main.tmax) {
          main.tmin = 1e9;
          return true;
        }
      }
    } else if (main.tmin == tmp.tmin) {
      if (main.tmax > tmp.tmax) {
        main = tmp;
      }
    } else if (main.tmax > tmp.tmax) {
      main.tmax = tmp.tmax;
      if (main.tmin > main.tmax) {
        main.tmin = 1e9;
        return true;
      }
    }
    return false;
  }

public:
  // TODO: OPTIMIZE THE BELOW BY CHANGING BASES
  /// we suppose that p moves linearly by pt, and q by qt
  Contact getInt(const v::DVec<3> &pt, const v::DVec<3> &qt) {
    std::cout << "ENTERING!!" << std::endl;
    v::DVec<3> pA[] = {p.x, p.y, p.z}; // p vert-face
    v::DVec<3> qA[] = {q.x, q.y, q.z}; // q vert-face
    v::DVec<3> eA[] = {ee(p.x, q.x), ee(p.x, q.y), ee(p.x, q.z),
                       ee(p.y, q.x), ee(p.y, q.y), ee(p.y, q.z),
                       ee(p.z, q.x), ee(p.z, q.y), ee(p.z, q.z)};
    v::DVec<3> relVelo = qt - pt;

    Contact ret;
    AxisIDetail main{-1, 1e9, 1e9, {0, 0, 0}, 0, 0};
    for (int i = 0; i < 3; i++) {
      v::DVec<2> a = p.extrema(pA[i]);
      v::DVec<2> b = q.extrema(pA[i]);
      double tmpRel = v::dot(pA[i], relVelo);
      v::DVec<3> mm = intervalInt(tmpRel, a[0] - b[1], a[1] - b[0]);
      std::cout << "a,b,tmpRel,mm: " << a << b << tmpRel << " " << mm
                << std::endl;
      // we don't use the index anyway in this case
      if (updateAID({mm[0], mm[1], mm[2],
                     tmpRel > 0 ? pA[i] : v::DVec<3>(-pA[i]), 0, 0},
                    main)) {
        break;
      }
      a = p.extrema(qA[i]);
      b = q.extrema(qA[i]);
      tmpRel = v::dot(qA[i], relVelo);
      mm = intervalInt(tmpRel, a[0] - b[1], a[1] - b[0]);
      if (updateAID({mm[0], mm[1], mm[2],
                     tmpRel > 0 ? qA[i] : v::DVec<3>(-qA[i]), 1, 0},
                    main)) {
        break;
      }
    }
    if (main.tmin > main.tmax) {
      ret.t = main.tmin;
      ret.n = main.n;
      std::cout << "kek1: " << ret.t << std::endl;
      return ret;
    }
    for (int i = 0; i < 9; i++) {
      v::DVec<2> a = p.extrema(eA[i]);
      v::DVec<2> b = q.extrema(eA[i]);
      double tmpRel = v::dot(eA[i], relVelo);
      v::DVec<3> mm = intervalInt(tmpRel, a[0] - b[1], a[1] - b[0]);
      if (updateAID({mm[0], mm[1], mm[2],
                     tmpRel > 0 ? eA[i] : v::DVec<3>(-eA[i]), 2, i},
                    main)) {
        break;
      }
    }
    ret.t = main.tmin;
    ret.n = main.n;
    if (main.tmin > main.tmax) {
      std::cout << "kek2: " << ret.t << std::endl;
      return ret;
    }
    ret.d = main.depth;

    OBB pp = p;
    pp.b += ret.t * pt;
    OBB qq = q;
    qq.b += ret.t * qt;
    switch (main.type) {
    case 0:
      std::cout << "HERE WE GO! " << qq.maximize(ret.n) << " " << qt
                << std::endl;
      ret.p = qq.maximize(ret.n) + ret.t * qt;
      return ret;
    case 1:
      std::cout << "maybe???! " << pp.maximize(-ret.n) << " " << pt
                << std::endl;
      ret.p = pp.maximize(-ret.n) + ret.t * pt;
      return ret;
    default:
      int qi = main.index % 3;
      v::DVec<3> ql = q.maxEdge(ret.n, qi);
      std::cout << "edge.......???! "
                << p.getAlongSeg(ql, ql + q.s[qi] * qA[qi]) << " " << pt
                << std::endl;
      ret.p = p.getAlongSeg(ql, ql + q.s[qi] * qA[qi]) + ret.t * pt;
      return ret;
    }
  }
};

// XXX: REUSE POSES AT ENDPOINTS OF ARCS
/// CCD for OBBs in screw motion
struct CCDRotOBBIntersector {

  v::DVec<3> ps, qs; // OBB sidelengths
  Pose p, q;
  v::DVec<3> w1, w2; // these are the angular velocities applied
  v::DVec<3> pc, qc; // centers of rotations

  double tol; /// permitted tolerance

  struct RotState {
    bool isInit = false;
    v::DVec<4> u1, u2;                 // quaternions for the below
    v::DVec<3> px, py, pz, qx, qy, qz; // x,y,z of the respective OBBs
  };

  // 1 means no subdividing, 2 means subdividing once in the middle, etc.
#define MAX_LEVELS 7
  // because we have to calculate using the orientation in the middle of the arc
#define MAX_SUBDIVISION (1 << MAX_LEVELS)
  RotState states[MAX_SUBDIVISION + 1];
  CCDOBBIntersector locals[MAX_SUBDIVISION];

  v::DVec<4> smolRot1[MAX_LEVELS];
  v::DVec<4> smolRot2[MAX_LEVELS];

  /// the pose.p is the OBB's b, not the center of mass
  CCDRotOBBIntersector(const Pose &p, v::DVec<3> ps, v::DVec<3> omega1,
                       v::DVec<3> pc, const Pose &q, v::DVec<3> qs,
                       v::DVec<3> omega2, v::DVec<3> qc, double tol)
      : ps(ps), qs(qs), p(p), q(q), w1(omega1 / 2), w2(omega2 / 2), pc(pc),
        qc(qc), tol(tol) {
    states[0].u1 = p.q;
    states[0].u2 = q.q;
    OBB tmp;
    setOBBOrientation(p, tmp);
    states[0].px = tmp.x;
    states[0].py = tmp.y;
    states[0].pz = tmp.z;
    setOBBOrientation(q, tmp);
    states[0].qx = tmp.x;
    states[0].qy = tmp.y;
    states[0].qz = tmp.z;
    states[0].isInit = true;
    for (int lvl = 0; lvl < MAX_LEVELS; lvl++) {
      smolRot1[lvl] =
          normalizeQuaternion(getRotQuaternion(w1 / (double)(1 << (lvl + 1))));
      smolRot2[lvl] =
          normalizeQuaternion(getRotQuaternion(w2 / (double)(1 << (lvl + 1))));
    }
  }

  // just the obb.b of both
  void updateOBBs(v::DVec<3> np, v::DVec<3> nps, v::DVec<3> nq,
                  v::DVec<3> nqs) {
    p.p = np;
    ps = nps;
    q.p = nq;
    qs = nqs;
  }

  static void setOBBOrientation(const Pose &pose, OBB &obb) {
    DMat3x3 ro = pose.toRotationMatrix().transpose();
    obb.x = ro.a;
    obb.y = ro.b;
    obb.z = ro.c;
  }
  static OBB poseToOBB(const Pose &pose, const v::DVec<3> &sl) {
    OBB ret;
    setOBBOrientation(pose, ret);
    ret.b = pose.p;
    ret.s = sl;
    ret.properify();
    return ret;
  }

  struct Subdiv {
    std::size_t ind, lvl;
    v::DVec<4> w1, w2; // quaternions to rotate half the arc
    Pose pp, qq;       // earlier end of arc
    OBB po, qo;        // obbs for pp and qq
  };
  // moves sd.pp,qq,po,and qo forward by half the arc, updates ind accordingly
  template <bool updateLocals> void modSubdiv(Subdiv &sd) {
    sd.pp.p = applyQuaternion(sd.pp.p, sd.w1);
    sd.qq.p = applyQuaternion(sd.qq.p, sd.w2);
    sd.po.b = sd.pp.p + pc;
    sd.qo.b = sd.qq.p + qc;
    std::size_t sii = sd.ind + (MAX_SUBDIVISION >> (sd.lvl + 1));
    sd.ind = sii;
    if (states[sii].isInit) {
      sd.pp.q = states[sii].u1;
      sd.qq.q = states[sii].u2;
      sd.po.x = states[sii].px;
      sd.po.y = states[sii].py;
      sd.po.z = states[sii].pz;
      sd.qo.x = states[sii].qx;
      sd.qo.y = states[sii].qy;
      sd.qo.z = states[sii].qz;
      if (updateLocals) {
        locals[sii].setOBBs(sd.po, sd.qo);
      }
    } else {
      sd.pp.q = normalizeQuaternion(quaternionMult(sd.w1, sd.pp.q));
      states[sii].u1 = sd.pp.q;
      sd.qq.q = normalizeQuaternion(quaternionMult(sd.w1, sd.qq.q));
      states[sii].u2 = sd.qq.q;
      setOBBOrientation(sd.pp, sd.po);
      states[sii].px = sd.po.x;
      states[sii].py = sd.po.y;
      states[sii].pz = sd.po.z;
      setOBBOrientation(sd.qq, sd.qo);
      states[sii].qx = sd.qo.x;
      states[sii].qy = sd.qo.y;
      states[sii].qz = sd.qo.z;
      if (updateLocals) {
        locals[sii].~CCDOBBIntersector();
        new (&locals[sii]) CCDOBBIntersector(sd.po, sd.qo);
      }
    }
  }
  /// pv and qv are velocities
  Contact getInt(v::DVec<3> pv, v::DVec<3> qv) {
    std::vector<Subdiv> stack;
    v::DVec<4> tmpw1 = smolRot1[0];
    v::DVec<4> tmpw2 = smolRot2[0];
    OBB tmppo, tmpqo;
    tmppo.x = states[0].px;
    tmppo.y = states[0].py;
    tmppo.z = states[0].pz;
    tmpqo.x = states[0].qx;
    tmpqo.y = states[0].qy;
    tmpqo.z = states[0].qz;
    tmppo = OBB(p.p, ps, states[0].px, states[0].py, states[0].pz);
    tmpqo = OBB(q.p, qs, states[0].qx, states[0].qy, states[0].qz);

    Pose tmpp = p;
    tmpp.p -= pc;
    Pose tmpq = q;
    tmpq.p -= qc;
    double radiusp = std::sqrt(v::norm2(tmpp.p));
    double radiusq = std::sqrt(v::norm2(tmpq.p));
    stack.push_back({0, 0, tmpw1, tmpw2, tmpp, tmpq, tmppo, tmpqo});
    do {
      Subdiv subdiv = stack.back();
      stack.pop_back();
      Subdiv subdivm = subdiv;
      modSubdiv<true>(subdivm);
      Subdiv subdivl = subdivm;
      modSubdiv<false>(subdivl);
      CCDOBBIntersector &inn = locals[subdivm.ind];

      OBB obb1 = subdivm.po.wrapOBB(subdiv.po);
      OBB obb2 = subdivm.qo.wrapOBB(subdiv.qo);
      OBB obb3 = subdivm.po.wrapOBB(subdivl.po);
      OBB obb4 = subdivm.qo.wrapOBB(subdivl.qo);
      v::DVec<3> pdisp = obb3.center() - obb1.center();
      obb3.b -= pdisp;
      obb1.fattenOBB(obb3);
      v::DVec<3> qdisp = obb4.center() - obb2.center();
      obb4.b -= qdisp;
      obb2.fattenOBB(obb4);

      v::DVec<3> pvelo = pv + pdisp;
      v::DVec<3> qvelo = qv + qdisp;

      double err2 = subdivm.w1[0] * subdivm.w1[0];
      double err2p = err2 * radiusp;
      double err2q = err2 * radiusq;
      // remember that theta^2/2>=1-cos(theta)
      double errp = err2p / 2;
      double errq = err2q / 2;
      obb1.b -= errp;
      obb1.s += err2p;
      obb2.b -= errq;
      obb2.s += err2q;

      Contact contact = inn.getInt(pvelo, qvelo);
      if (contact.t <= 1) {
        if ((errp < tol && errq < tol) || (subdivm.lvl + 1 >= MAX_LEVELS)) {
          return contact;
        }
        subdivm.lvl = ++subdiv.lvl;
        subdivm.w1 = subdiv.w1 = smolRot1[subdiv.lvl];
        subdivm.w2 = subdiv.w2 = smolRot2[subdiv.lvl];
        stack.push_back(subdivm);
        stack.push_back(subdiv);
      }
    } while (stack.size());
    Contact ret;
    // no collision
    ret.t = 2;
    return ret;
  }
};

} // namespace roller

#endif // ROLLER_PHYSICS_CCD_NARROW_HPP_

