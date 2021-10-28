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
      if (relVelo > 1e-12) {
        max = df2 / relVelo;
      } else if (relVelo < -1e-12) {
        max = df1 / relVelo;
      } else {
        max = 100;
      }
    } else {
      if (df1 > 0 && relVelo > 1e-12) {
        min = df1 / relVelo;
        max = df2 / relVelo;
      } else if (df2 < 0 && relVelo < -1e-12) {
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

public:
#ifdef TEST_BUG_DEBUG
  bool overlapKek = false;
#endif
  // TODO: OPTIMIZE THE BELOW BY CHANGING BASES
  /// we suppose that p moves linearly by pt, and q by qt
  Contact getInt(const v::DVec<3> &pt, const v::DVec<3> &qt) {
#ifdef TEST_BUG_DEBUG
    overlapKek = false;
#endif
    // std::cout << "ENTERING!!" << std::endl;
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
#ifdef TEST_BUG_DEBUG
      overlapKek = overlapKek || ((a[0] <= b[0] && a[1] >= b[1]) ||
                                  (a[0] >= b[0] && a[1] <= b[1]));
#endif
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
        break;
      }
      a = p.extrema(qA[i]);
      b = q.extrema(qA[i]);
#ifdef TEST_BUG_DEBUG
      overlapKek = overlapKek || ((a[0] <= b[0] && a[1] >= b[1]) ||
                                  (a[0] >= b[0] && a[1] <= b[1]));
#endif
      tmpRel = v::dot(qA[i], relVelo);
      mm = intervalInt(tmpRel, a[0] - b[1], a[1] - b[0]);
      if (mm[0] > 0) {
        ncond = tmpRel > 0;
      } else {
        ncond = a[1] - b[0] > b[1] - a[0];
      }
      normal = ncond ? qA[i] : v::DVec<3>(-qA[i]);
      if (updateAID({mm[0], mm[1], mm[2], normal, 1, 0}, main)) {
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
#ifdef TEST_BUG_DEBUG
      overlapKek = overlapKek || ((a[0] <= b[0] && a[1] >= b[1]) ||
                                  (a[0] >= b[0] && a[1] <= b[1]));
#endif
      double tmpRel = v::dot(eA[i], relVelo);
      v::DVec<3> mm = intervalInt(tmpRel, a[0] - b[1], a[1] - b[0]);
      bool ncond;
      if (mm[0] > 0) {
        ncond = tmpRel > 0;
      } else {
        ncond = a[1] - b[0] > b[1] - a[0];
      }
      v::DVec<3> normal = ncond ? eA[i] : v::DVec<3>(-eA[i]);
      if (updateAID({mm[0], mm[1], mm[2], normal, 2, i}, main)) {
        break;
      }
    }
    ret.t = main.tmin;
    ret.n = main.n;
    if (main.tmin > main.tmax) {
      return ret;
    }
    ret.d = main.depth;

    switch (main.type) {
    case 0:
      ret.p = q.maximize(ret.n) + ret.t * qt;
      return ret;
    case 1:
      ret.p = p.maximize(-ret.n) + ret.t * pt;
      return ret;
    default:
      int qi = main.index % 3;
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
    std::size_t lvll = MAX_SUBDIVISION >> (sd.lvl + 1);
    sd.ind += lvll;
    std::cout << "modSubdiv: " << lvll << " " << sd.ind << std::endl;
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
      sd.qq.q = normalizeQuaternion(quaternionMult(sd.w1, sd.qq.q));
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
        std::cout << "modSubdiv, initializing local " << sd.ind << std::endl;
        locals[sd.ind].~CCDOBBIntersector();
        new (&locals[sd.ind]) CCDOBBIntersector(sd.po, sd.qo);
      }
    }
  }
#ifdef TEST_BUG_DEBUG
  bool overlapKek = false;
#endif
  /// pv and qv are velocities
  Contact getInt(v::DVec<3> pv, v::DVec<3> qv) {
    std::cout << "getInt entry here-----------" << std::endl;

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

    std::cout << "tmppo: " << tmppo.b << tmppo.s << tmppo.x << tmppo.y
              << tmppo.z << std::endl;
    std::cout << "tmpqo: " << tmpqo.b << tmpqo.s << tmpqo.x << tmpqo.y
              << tmpqo.z << std::endl;

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
      std::cout << "do-while." << std::endl;

      std::cout << subdiv.ind << " " << subdiv.lvl << std::endl;
      std::cout << subdiv.w1 << " " << subdiv.w2 << std::endl;
      std::cout << subdiv.pp.p << subdiv.pp.q << subdiv.qq.p << subdiv.qq.q
                << std::endl;
      std::cout << subdiv.po.b << subdiv.po.s << subdiv.po.x << subdiv.po.y
                << subdiv.po.z << std::endl;
      std::cout << subdiv.qo.b << subdiv.qo.s << subdiv.qo.x << subdiv.qo.y
                << subdiv.qo.z << std::endl;
      std::cout << "subdivl: " << std::endl;
      std::cout << subdivl.ind << " " << subdivl.lvl << std::endl;
      std::cout << subdivl.w1 << " " << subdivl.w2 << std::endl;
      std::cout << subdivl.pp.p << subdivl.pp.q << subdivl.qq.p << subdivl.qq.q
                << std::endl;
      std::cout << subdivl.po.b << subdivl.po.s << subdivl.po.x << subdivl.po.y
                << subdivl.po.z << std::endl;
      std::cout << subdivl.qo.b << subdivl.qo.s << subdivl.qo.x << subdivl.qo.y
                << subdivl.qo.z << std::endl;
      /*std::cout << "adjusted?" << std::endl;
      OBB utterspam1 = subdivl.po;
      OBB utterspam2 = subdivl.qo;
      utterspam1.b += (((double)subdivl.ind) / MAX_SUBDIVISION) * pv;
      utterspam2.b += (((double)subdivl.ind) / MAX_SUBDIVISION) * qv;
      utterspam1.properify();
      utterspam2.properify();
      std::cout << utterspam1.b << utterspam1.s << utterspam1.x << utterspam1.y
                << utterspam1.z << utterspam1.a << utterspam1.c << std::endl;
      std::cout << utterspam2.b << utterspam2.s << utterspam2.x << utterspam2.y
                << utterspam2.z << utterspam2.a << utterspam2.c << std::endl;*/

      OBB obb1 = subdivm.po.wrapOBB(subdiv.po);
      OBB obb2 = subdivm.qo.wrapOBB(subdiv.qo);
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

      // HORRIBLE DEBUG HERE!!!!!!!!!
      /*
      double step = 1 / 50.;
      double mval = 0;
      for (double vv1 = 0; vv1 < 1; vv1 += step) {
        for (double vv2 = 0; vv2 < 1; vv2 += step) {
          for (int type1 = 0; type1 < 6; type1++) {
            double mmval = 10000000;
            for (double vv4 = 0; vv4 < 1; vv4 += step) {
              for (double vv5 = 0; vv5 < 1; vv5 += step) {
                for (int type2 = 0; type2 < 6; type2++) {
                  double v1, v2, v3, v4, v5, v6;
                  switch (type1) {
                  case 0:
                    v1 = vv1;
                    v2 = vv2;
                    v3 = 0;
                    break;
                  case 1:
                    v1 = vv1;
                    v2 = vv2;
                    v3 = 1;
                    break;
                  case 2:
                    v1 = vv1;
                    v2 = 0;
                    v3 = vv2;
                    break;
                  case 3:
                    v1 = vv1;
                    v2 = 1;
                    v3 = vv2;
                    break;
                  case 4:
                    v1 = 0;
                    v2 = vv1;
                    v3 = vv2;
                    break;
                  default:
                    v1 = 1;
                    v2 = vv1;
                    v3 = vv2;
                  }
                  switch (type2) {
                  case 0:
                    v4 = vv4;
                    v5 = vv5;
                    v6 = 0;
                    break;
                  case 1:
                    v4 = vv4;
                    v5 = vv5;
                    v6 = 1;
                    break;
                  case 2:
                    v4 = vv4;
                    v5 = 0;
                    v6 = vv5;
                    break;
                  case 3:
                    v4 = vv4;
                    v5 = 1;
                    v6 = vv5;
                    break;
                  case 4:
                    v4 = 0;
                    v5 = vv4;
                    v6 = vv5;
                    break;
                  default:
                    v4 = 1;
                    v5 = vv4;
                    v6 = vv5;
                  }
                  v::DVec<3> pnt1 = obb1.b + v1 * obb1.s[0] * obb1.x +
                                    v2 * obb1.s[1] * obb1.y +
                                    v3 * obb1.s[2] * obb1.z;
                  v::DVec<3> pnt2 = subdiv.po.b +
                                    v4 * subdiv.po.s[0] * subdiv.po.x +
                                    v5 * subdiv.po.s[1] * subdiv.po.y +
                                    v6 * subdiv.po.s[2] * subdiv.po.z;
                  double nval = v::norm2(pnt2 - pnt1);
                  if (nval < mmval) {
                    mmval = nval;
                  }
                }
              }
            }
            if (mmval > mval) {
              mval = mmval;
            }
          }
        }
      }
      std::cout << "HARD WORK HERE!!!!!!!!!" << std::sqrt(mval) << std::endl;
      */
      // END HORRIBLE DEBUG HERE!!!!!!!!!!!!!
      std::cout << "difference? " << obb1.s - subdivm.po.s << std::endl;

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

      // DEBUG HERE
      obb1.properify();
      obb2.properify();
      std::cout << "pvelo,qvelo: " << pvelo << " " << qvelo << std::endl;
      std::cout << "subdivm.w1[0],subdivm.w2[0]: " << subdivm.w1[0] << " "
                << subdivm.w2[0] << std::endl;
      std::cout << "errp,errq: " << errp << " " << errq << std::endl;
      std::cout << obb1.b << obb1.s << obb1.x << obb1.y << obb1.z << obb1.a
                << obb1.c << std::endl;
      std::cout << obb2.b << obb2.s << obb2.x << obb2.y << obb2.z << obb2.a
                << obb2.c << std::endl;

      inn.setOBBs(obb1, obb2);
      Contact contact = inn.getInt(pvelo, qvelo);
      if (contact.t <= 1) {
        double tmp1 = std::sin(std::sqrt(v::norm2(w1)) / timedivf);
        double tmp2 = std::sin(std::sqrt(v::norm2(w2)) / timedivf);
        // double tmp2 = std::sqrt(1 - cos2 * cos2);

        // sum here is fat overestimate, but why not
        double errpt = 2 * (tmp1 * v::sum(obb1.s) + errp);
        double errqt = 2 * (tmp2 * v::sum(obb2.s) + errq);
        std::cout << "errpt,errqt: " << errpt << " " << errqt << std::endl;
        if ((errpt < tol && errqt < tol) || (subdivm.lvl + 1 >= MAX_LEVELS)) {
#ifdef TEST_BUG_DEBUG
          overlapKek = inn.overlapKek;
#endif
          contact.t /= timedivf;
          contact.t += indf;
          std::cout << "found! " << contact.t << " " << subdiv.ind << " "
                    << subdiv.lvl << std::endl;
          return contact;
        }
        subdivm.lvl = ++subdiv.lvl;
        subdivm.w1 = subdiv.w1 = smolRot1[subdiv.lvl];
        subdivm.w2 = subdiv.w2 = smolRot2[subdiv.lvl];
        std::cout << "pushing: " << subdivm.ind << " " << subdivm.lvl
                  << std::endl;
        stack.push_back(subdivm);
        std::cout << "pushing: " << subdiv.ind << " " << subdiv.lvl
                  << std::endl;
        stack.push_back(subdiv);
      }
      std::cout << std::endl;
    } while (stack.size());
    std::cout << "no collision" << std::endl;
    Contact ret;
    // no collision
    ret.t = 2;
    return ret;
  }
};

} // namespace roller

#endif // ROLLER_PHYSICS_CCD_NARROW_HPP_

