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

  DMat3x3 qrelori; /// orientation of q relative to p

  CCDOBBIntersector() {}
  CCDOBBIntersector(const OBB &p, const OBB &q) { initInt(p, q); }

  void initInt(const OBB &pp, const OBB &qq) {
    p = pp;
    q = qq;
    qrelori = DMat3x3{p.x, p.y, p.z} * DMat3x3{q.x, q.y, q.z}.transpose();
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
        max = 100;
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
      double tmax = main.tmax > tmp.tmax ? tmp.tmax : main.tmax;
      if (main.depth > tmp.depth) {
        main = tmp;
      }
      main.tmax = tmax;
    } else if (main.tmax > tmp.tmax) {
      main.tmax = tmp.tmax;
      if (main.tmin > main.tmax) {
        main.tmin = 1e9;
        return true;
      }
    }
    return false;
  }

  // variation of OBB::extrema
  static v::DVec<2> dimExtrema(double base, const v::DVec<3> (&axes)[3],
                               const v::DVec<3> &sl, int dim) {
    double a = base, b = base;
    if (axes[0][dim] < 0) {
      a += sl[0] * axes[0][dim];
    } else {
      b += sl[0] * axes[0][dim];
    }
    if (axes[1][dim] < 0) {
      a += sl[1] * axes[1][dim];
    } else {
      b += sl[1] * axes[1][dim];
    }
    if (axes[2][dim] < 0) {
      a += sl[2] * axes[2][dim];
    } else {
      b += sl[2] * axes[2][dim];
    }
    return {a, b};
  }
  static v::DVec<2> eeExtrema(double base, const v::DVec<2> &axis,
                              const v::DVec<3> &sl, int dim1, int dim2) {
    double a = base, b = base;
    if (axis[0] < 0) {
      a += sl[dim1] * axis[0];
    } else {
      b += sl[dim1] * axis[0];
    }
    if (axis[1] < 0) {
      a += sl[dim2] * axis[1];
    } else {
      b += sl[dim2] * axis[1];
    }
    return {a, b};
  }

public:
  // Further work: Identify precise reason for the ridiculously rare slightly
  // off error we get in the test (that I don't care much about btw)

  /// we suppose that p moves linearly by pt, and q by qt
  Contact getInt(const v::DVec<3> &pt, const v::DVec<3> &qt) {
    v::DVec<3> pA[] = {p.x, p.y, p.z}; // p vert-face
    v::DVec<3> qA[] = {q.x, q.y, q.z}; // q vert-face
    v::DVec<3> relVelo = qt - pt;

    // face normals relative to the other OBB
    v::DVec<3> prA[] = {qrelori.a, qrelori.b, qrelori.c};
    DMat3x3 prelori = qrelori.transpose();
    v::DVec<3> qrA[] = {prelori.a, prelori.b, prelori.c};

    Contact ret;
    AxisIDetail main{-1, 1e9, 1e9, {0, 0, 0}, 0, 0};
    for (int i = 0; i < 3; i++) {
      v::DVec<2> a = {p.a[i], p.c[i]};
      v::DVec<2> b = dimExtrema(v::dot(pA[i], q.b), qrA, q.s, i);
      double tmpRel = v::dot(pA[i], relVelo);
      v::DVec<3> mm = intervalInt(tmpRel, a[0] - b[1], a[1] - b[0]);
      bool ncond;
      if (mm[0] > 0) {
        ncond = tmpRel > 0;
      } else {
        ncond = a[1] - b[0] > b[1] - a[0];
      }
      v::DVec<3> normal = ncond ? pA[i] : v::DVec<3>(-pA[i]);
      // we don't use the index anyway in this case
      if (updateAID({mm[0], mm[1], mm[2], normal, 0, 0}, main)) {
        ret.t = main.tmin;
        ret.n = main.n;
        return ret;
      }
      a = dimExtrema(v::dot(qA[i], p.b), prA, p.s, i);
      b = {q.a[i], q.c[i]};
      tmpRel = v::dot(qA[i], relVelo);
      mm = intervalInt(tmpRel, a[0] - b[1], a[1] - b[0]);
      if (mm[0] > 0) {
        ncond = tmpRel > 0;
      } else {
        ncond = a[1] - b[0] > b[1] - a[0];
      }
      normal = ncond ? qA[i] : v::DVec<3>(-qA[i]);
      if (updateAID({mm[0], mm[1], mm[2], normal, 1, 0}, main)) {
        ret.t = main.tmin;
        ret.n = main.n;
        return ret;
      }
    }
    int others[3][2] = {{1, 2}, {2, 0}, {0, 1}}; // for cross-product
    for (int indp = 0; indp < 3; indp++) {
      for (int indq = 0; indq < 3; indq++) {
        if (1 - qrA[indq][indp] < 1e-8) {
          // degerate cross product
          continue;
        }
        int o0p = others[indp][0];
        int o1p = others[indp][1];
        int o0q = others[indq][0];
        int o1q = others[indq][1];
        // lol this is not a unit vector
        v::DVec<3> cros = cross3(pA[indp], qA[indq]);
        v::DVec<2> a =
            eeExtrema(v::dot(cros, p.b), {-qrA[indq][o1p], qrA[indq][o0p]}, p.s,
                      o0p, o1p);
        v::DVec<2> b =
            eeExtrema(v::dot(cros, q.b), {prA[indp][o1q], -prA[indp][o0q]}, q.s,
                      o0q, o1q);
        double tmpRel = v::dot(cros, relVelo);
        v::DVec<3> mm = intervalInt(tmpRel, a[0] - b[1], a[1] - b[0]);
        bool ncond;
        if (mm[0] > 0) {
          ncond = tmpRel > 0;
        } else {
          ncond = a[1] - b[0] > b[1] - a[0];
        }
        v::DVec<3> normal = ncond ? cros : v::DVec<3>(-cros);
        // note: might remove indp from index, since I'm not using it anyways
        if (updateAID({mm[0], mm[1], mm[2], normal, 2, (indp << 2) | indq},
                      main)) {
          ret.t = main.tmin;
          ret.n = main.n;
          return ret;
        }
      }
    }
    ret.t = main.tmin;
    ret.n = main.n;
    ret.d = main.depth;

    // note: this calculation of contact point would be really bad if
    // this were not CCD, but based on testing, I think the degenerate
    // cases are really rare
    switch (main.type) {
    case 0:
      ret.p = q.maximize(ret.n) + ret.t * qt;
      return ret;
    case 1:
      ret.p = p.maximize(-ret.n) + ret.t * pt;
      return ret;
    default:
      int qi = main.index & 3;
      double invsqn = fastInvSqrt(v::norm2(ret.n));
      ret.n *= invsqn;
      // talk about trippy
      ret.d *= invsqn;
      v::DVec<3> ql = q.maxEdge(ret.n, qi) + ret.t * (qt - pt);
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
  v::DVec<3> w1, w2; // these are the angular velocities (/2 in constructor)
  double w1n, w2n;   // norms of w1 and w2
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

  /// the pose.p is the OBB's b, not the center of mass (tol>2e-5)
  CCDRotOBBIntersector(const Pose &p, v::DVec<3> ps, v::DVec<3> omega1,
                       v::DVec<3> pc, const Pose &q, v::DVec<3> qs,
                       v::DVec<3> omega2, v::DVec<3> qc, double tol)
      : ps(ps), qs(qs), p(p), q(q), w1(omega1 / 2), w2(omega2 / 2),
        w1n(std::sqrt(v::norm2(w1))), w2n(std::sqrt(v::norm2(w2))), pc(pc),
        qc(qc), tol(tol - 1e-5) {
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
    std::size_t lvll = MAX_SUBDIVISION >> (sd.lvl + 1);
    sd.ind += lvll;
    if (states[sd.ind].isInit) {
      sd.pp.q = states[sd.ind].u1;
      sd.qq.q = states[sd.ind].u2;
      sd.po.x = states[sd.ind].px;
      sd.po.y = states[sd.ind].py;
      sd.po.z = states[sd.ind].pz;
      sd.qo.x = states[sd.ind].qx;
      sd.qo.y = states[sd.ind].qy;
      sd.qo.z = states[sd.ind].qz;
    } else {
      sd.pp.q = normalizeQuaternion(quaternionMult(sd.w1, sd.pp.q));
      states[sd.ind].u1 = sd.pp.q;
      sd.qq.q = normalizeQuaternion(quaternionMult(sd.w2, sd.qq.q));
      states[sd.ind].u2 = sd.qq.q;
      setOBBOrientation(sd.pp, sd.po);
      states[sd.ind].px = sd.po.x;
      states[sd.ind].py = sd.po.y;
      states[sd.ind].pz = sd.po.z;
      setOBBOrientation(sd.qq, sd.qo);
      states[sd.ind].qx = sd.qo.x;
      states[sd.ind].qy = sd.qo.y;
      states[sd.ind].qz = sd.qo.z;
      states[sd.ind].isInit = true;
      if (updateLocals) {
        locals[sd.ind].~CCDOBBIntersector();
        new (&locals[sd.ind]) CCDOBBIntersector(sd.po, sd.qo);
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
      // XXX: CHECK TO SEE IF THIS IS REALLY NECESSARY
      OBB obb3 = subdivm.po.wrapOBB(subdivl.po);
      OBB obb4 = subdivm.qo.wrapOBB(subdivl.qo);
      v::DVec<3> pdisp = obb3.center() - obb1.center();
      obb3.b -= pdisp;
      obb3.properify();
      obb1.fattenOBB(obb3);
      v::DVec<3> qdisp = obb4.center() - obb2.center();
      obb4.b -= qdisp;
      obb4.properify();
      obb2.fattenOBB(obb4);

      double timedivf = 1 << subdiv.lvl;
      v::DVec<3> pvelo = pv / timedivf + pdisp;
      v::DVec<3> qvelo = qv / timedivf + qdisp;

      // recall that subdivm.w1/w2 only rotates by half the desired arc
      double cos1 = 2 * subdivm.w1[0] * subdivm.w1[0] - 1; // dbl angle formula
      double cos2 = 2 * subdivm.w2[0] * subdivm.w2[0] - 1;
      double serr1 = 1 - cos1;
      double serr2 = 1 - cos2;
      double errp = serr1 * radiusp;
      double errq = serr2 * radiusq;
      obb1.inflate(errp);
      obb2.inflate(errq);

      double indf = subdiv.ind;
      indf /= MAX_SUBDIVISION;
      obb1.b += pv * indf;
      obb2.b += qv * indf;

      inn.setOBBs(obb1, obb2);
      Contact contact = inn.getInt(pvelo, qvelo);
      if (contact.t <= 1) {
        double tmp1 = w1n / timedivf;
        double tmp2 = w2n / timedivf;
        // sum here (also the 2*) is fat overestimate, but why not
        double errpt = 2 * (tmp1 * v::sum(obb1.s) + errp);
        double errqt = 2 * (tmp2 * v::sum(obb2.s) + errq);
        if ((errpt < tol && errqt < tol) || (subdivm.lvl + 1 >= MAX_LEVELS)) {
          contact.t /= timedivf;
          contact.t += indf;
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

