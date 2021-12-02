#include <iostream>
#include <random>

#define TEST_BUG_DEBUG
#include "ccd_narrow.hpp"

// #include "old_ccd.hpp"

using namespace roller;

bool stupidOBBInt(const OBB &a, const OBB &b) {
  struct Plane {
    v::DVec<3> v;
    double p;

    v::DVec<3> get(const Plane &a, const Plane &b) const {
      if (v::norm2(v - a.v) < 1e-20 || v::norm2(v - b.v) < 1e-20 ||
          v::norm2(a.v - b.v) < 1e-20) {
        // all my tests will be within 100 of the origin
        return {1000, 1000, 1000};
      }
      return DMat3x3{v, a.v, b.v}.solve({p, a.p, b.p});
    }
  };
  v::DVec<3> aa = a.a + 1e-8;
  v::DVec<3> ac = a.c - 1e-8;
  v::DVec<3> ba = b.a + 1e-8;
  v::DVec<3> bc = b.c - 1e-8;
  Plane planes[] = {{a.x, aa[0]}, {a.x, ac[0]}, {a.y, aa[1]}, {a.y, ac[1]},
                    {a.z, aa[2]}, {a.z, ac[2]}, {b.x, ba[0]}, {b.x, bc[0]},
                    {b.y, ba[1]}, {b.y, bc[1]}, {b.z, ba[2]}, {b.z, bc[2]}};
  for (int i = 0; i < 10; i++) {
    for (int j = i + 1; j < 11; j++) {
      for (int k = j + 1; k < 12; k++) {
        v::DVec<3> pp = planes[i].get(planes[j], planes[k]);
        if (a.isIn(pp) && b.isIn(pp)) {
          // std::cout << "stupid: " << pp << std::endl;
          return true;
        }
      }
    }
  }
  return false;
}

void apply(OBB &p, v::DVec<3> trans, v::DVec<3> omega, v::DVec<3> c) {
  omega /= 2;
  v::DVec<4> orot = getRotQuaternion(omega);
  p.b = trans + c + applyQuaternion(p.b - c, orot);
  p.x = applyQuaternion(p.x, orot);
  p.y = applyQuaternion(p.y, orot);
  p.z = applyQuaternion(p.z, orot);
  p.properify();
}
void apply2Pose(Pose &p, v::DVec<3> trans, v::DVec<3> omega, v::DVec<3> c) {
  omega /= 2;
  v::DVec<4> orot = getRotQuaternion(omega);
  p.p = trans + c + applyQuaternion(p.p - c, orot);
  p.q = normalizeQuaternion(quaternionMult(orot, p.q));
}
void apply(OBB &p, const ScrewM &pm) { apply(p, pm.velo, pm.omega, pm.center); }
void apply2Pose(Pose &p, const ScrewM &pm) {
  apply2Pose(p, pm.velo, pm.omega, pm.center);
}

double getTime(const OBB &p, v::DVec<3> pv, v::DVec<3> omega1, v::DVec<3> pc,
               const OBB &q, v::DVec<3> qv, v::DVec<3> omega2, v::DVec<3> qc,
               double dt) {
  for (double t = 0; t < 1; t += dt) {
    OBB pp = p;
    apply(pp, t * pv, t * omega1, pc);
    OBB qq = q;
    apply(qq, t * qv, t * omega2, qc);
    if (stupidOBBInt(pp, qq)) {
      return t;
    }
  }
  return 1.0001;
}
double getTime(const OBB &p, const ScrewM &pm, const OBB &q, const ScrewM &qm,
               double dt) {
  return getTime(p, pm.velo, pm.omega, pm.center, q, qm.velo, qm.omega,
                 qm.center, dt);
}

int testRot(int numIters, int hpthresh, int stupidthresh) {
  std::uniform_real_distribution<> distro(-10, 10);
  std::uniform_real_distribution<> smol(-0.05, 0.05);
  std::uniform_real_distribution<> pdistro(0, 10);
  for (int spam = 0; spam < numIters; spam++) {
    std::mt19937 mtrand(spam);
    if (spam % 1000 == 0) {
      std::cout << "Iteration #" << spam << std::endl;
    }
    Pose pp;
    pp.p = {distro(mtrand), distro(mtrand), distro(mtrand)};
    pp.q = normalizeQuaternion(
        {distro(mtrand), distro(mtrand), distro(mtrand), distro(mtrand)});
    Pose qp;
    qp.p = {distro(mtrand), distro(mtrand), distro(mtrand)};
    qp.q = normalizeQuaternion(
        {distro(mtrand), distro(mtrand), distro(mtrand), distro(mtrand)});

    v::DVec<3> pv{distro(mtrand), distro(mtrand), distro(mtrand)};
    v::DVec<3> qv{distro(mtrand), distro(mtrand), distro(mtrand)};
    v::DVec<3> omega1{smol(mtrand), smol(mtrand), smol(mtrand)};
    v::DVec<3> omega2{smol(mtrand), smol(mtrand), smol(mtrand)};
    v::DVec<3> pc{distro(mtrand), distro(mtrand), distro(mtrand)};
    v::DVec<3> qc{distro(mtrand), distro(mtrand), distro(mtrand)};
    ScrewM pm = {pv, omega1, pc};
    ScrewM qm = {qv, omega2, qc};

    OBB p =
        ccd::poseToOBB(pp, {pdistro(mtrand), pdistro(mtrand), pdistro(mtrand)});
    OBB q =
        ccd::poseToOBB(qp, {pdistro(mtrand), pdistro(mtrand), pdistro(mtrand)});

    /*
    CCDRotOBBIntersector inn(pp, p.s, pm.omega, pm.center, qp, q.s, qm.omega,
                             qm.center, 1e-2);
    Contact rescc = inn.getInt(pm.velo, qm.velo);
    */

    ccd::Contact resc = ccd::obbRot(pp, p.s, pm, qp, q.s, qm, 1e-2);
    double res1 = resc.t;
    if (res1 > 1) {
      res1 = 1.0001;
    }
    double res2 = res1;
    double dt = 0.01;
    double dttol = 0.03;
    if (spam < hpthresh) {
      dt = 0.001;
      dttol = 0.02;
    }
    if (spam < stupidthresh) {
      res2 = getTime(p, pm, q, qm, dt);
    }

    /*
    OBB spamp = p;
    apply(spamp, pv, omega1, pc);
    OBB spamq = q;
    apply(spamq, qv, omega2, qc);
    std::cout << std::endl << "spamp, spamq:" << std::endl;
    std::cout << spamp.b << spamp.s << spamp.x << spamp.y << spamp.z
              << spamp.a << spamp.c << std::endl;
    std::cout << spamq.b << spamq.s << spamq.x << spamq.y << spamq.z
              << spamq.a << spamq.c << std::endl;
              */

    OBB op = p;
    apply(op, res1 * pv, res1 * omega1, pc);
    OBB oq = q;
    apply(oq, res1 * qv, res1 * omega2, qc);
    op.a -= 0.01;
    op.c += 0.01;
    oq.a -= 0.01;
    oq.c += 0.01;
    double df = res1 > res2 ? res1 - res2 : res2 - res1;
    bool kekis =
        res1 != 0 && res1 <= 1 && (!op.isIn(resc.p) || !oq.isIn(resc.p));
    if (df >= dttol || kekis) {
      if (res1 <= 1 && res2 >= 1 && !kekis) {
        std::cout << "false positive..." << std::endl;
        continue;
      }
      std::cout << op.b << op.s << op.x << op.y << op.z << op.a << op.c
                << std::endl;
      std::cout << oq.b << oq.s << oq.x << oq.y << oq.z << oq.a << oq.c
                << std::endl;
      std::cout << v::dot(op.x, resc.p) << " " << v::dot(op.y, resc.p) << " "
                << v::dot(op.z, resc.p) << std::endl;
      std::cout << v::dot(oq.x, resc.p) << " " << v::dot(oq.y, resc.p) << " "
                << v::dot(oq.z, resc.p) << std::endl;
      Pose pp1 = pp;
      Pose qp1 = qp;
      apply2Pose(pp1, res1 * pv, res1 * omega1, pc);
      apply2Pose(qp1, res1 * qv, res1 * omega2, qc);
      Pose pp2 = pp;
      Pose qp2 = qp;
      apply2Pose(pp2, res2 * pv, res2 * omega1, pc);
      apply2Pose(qp2, res2 * qv, res2 * omega2, qc);
      std::cout << pp1.p << pp1.q << std::endl;
      std::cout << qp1.p << qp1.q << std::endl;
      std::cout << pp2.p << pp2.q << std::endl;
      std::cout << qp2.p << qp2.q << std::endl;
      std::cout << resc.p << resc.d << " " << resc.n << std::endl;

      std::cerr << "FAIL!" << std::endl;
      std::cerr << "Iteration #" << spam << std::endl;
      std::cerr << "obbRot result: " << res1 << std::endl;
      std::cerr << "stupid result: " << res2 << std::endl;
      std::cerr << "pp.p, pp.q: " << pp.p << pp.q << std::endl;
      std::cerr << "qp.p, qp.q: " << qp.p << qp.q << std::endl;
      std::cerr << "ps, qs: " << p.s << q.s << std::endl;
      std::cerr << "pv, pc, omega1: " << pv << pc << omega1 << std::endl;
      std::cerr << "qv, qc, omega2: " << qv << qc << omega2 << std::endl;
      return 1;
    }
  }
  return 0;
}
int testIntervalInt() {
  struct BasicTest {
    double operator()(double a, double b, double v1, double c, double d,
                      double v2) const {
      for (double t = 0; t <= 1; t += 0.0001) {
        double am = a + t * v1;
        double bm = b + t * v1;
        double cm = c + t * v2;
        double dm = d + t * v2;
        // std::cout << am << " " << bm << " " << cm << " " << dm << std::endl;
        if (am <= dm && bm >= cm) {
          return t;
        }
      }
      return 1.0001;
    }
  } basicTest;
  std::mt19937 mtrand(1);
  std::uniform_real_distribution<> distro(-100, 100);
  std::uniform_real_distribution<> pdistro(0, 50);
  for (int spam = 0; spam < 100000; spam++) {
    if (spam % 1000 == 0) {
      std::cout << "Iteration #" << spam << std::endl;
    }
    double a = distro(mtrand);
    double b = a + pdistro(mtrand);
    double v1 = distro(mtrand);
    double c = distro(mtrand);
    double d = c + pdistro(mtrand);
    double v2 = distro(mtrand);
    v::DVec<3> res3 = ccd::intervalInt(v2 - v1, a - d, b - c);
    double res1 = res3[0];
    if (res1 > 1) {
      res1 = 1.0001;
    }
    double res2 = basicTest(a, b, v1, c, d, v2);
    double df = res2 > res1 ? res2 - res1 : res1 - res2;
    if (df >= 0.001) {
      std::cerr << "FAILURE!!!" << std::endl;
      std::cerr << "res1, res2: " << res1 << " " << res2 << std::endl;
      std::cerr << "a,b,v1,c,d,v2: " << a << " " << b << " " << v1 << " " << c
                << " " << d << " " << v2 << std::endl;
      return 1;
    }
  }
  return 0;
}
int testLinear(int numIters, int hpthresh, int stupidthresh) {
  std::uniform_real_distribution<> distro(-100, 100);
  std::uniform_real_distribution<> smol(-0.05, 0.05);
  std::uniform_real_distribution<> pdistro(0, 100);
  for (int spam = 0; spam < numIters; spam++) {
    std::mt19937 mtrand(spam + 234567);
    if (spam % 10000 == 0) {
      std::cout << "Iteration #" << spam << std::endl;
    }
    Pose pp;
    pp.p = {distro(mtrand), distro(mtrand), distro(mtrand)};
    pp.q = normalizeQuaternion(
        {distro(mtrand), distro(mtrand), distro(mtrand), distro(mtrand)});
    Pose qp;
    qp.p = {distro(mtrand), distro(mtrand), distro(mtrand)};
    qp.q = normalizeQuaternion(
        {distro(mtrand), distro(mtrand), distro(mtrand), distro(mtrand)});
    v::DVec<3> pv{distro(mtrand), distro(mtrand), distro(mtrand)};
    v::DVec<3> qv{distro(mtrand), distro(mtrand), distro(mtrand)};
    OBB p =
        ccd::poseToOBB(pp, {pdistro(mtrand), pdistro(mtrand), pdistro(mtrand)});
    OBB q =
        ccd::poseToOBB(qp, {pdistro(mtrand), pdistro(mtrand), pdistro(mtrand)});
    ccd::Contact resc = ccd::obbLinear(p, pv, q, qv);
    // std::cout << "resc.p: " << resc.p << std::endl;
    double res1 = resc.t;
    if (res1 > 1) {
      res1 = 1.0001;
    }
    v::DVec<3> zero{0, 0, 0};
    double res2 = res1;
    double dt = 0.05;
    double dttol = 0.1;
    if (spam < hpthresh) {
      // because I want to at least do some higher-precision/slow tests
      dt = 0.001;
      dttol = 0.002;
    }
    if (spam < stupidthresh) {
      res2 = getTime(p, pv, zero, zero, q, qv, zero, zero, dt);
    }

    OBB op = p;
    apply(op, res1 * pv, zero, zero);
    OBB oq = q;
    apply(oq, res1 * qv, zero, zero);
    op.a -= 1e-4;
    op.c += 1e-4;
    oq.a -= 1e-4;
    oq.c += 1e-4;
    double df = res1 > res2 ? res1 - res2 : res2 - res1;
    bool kekis =
        res1 != 0 && res1 <= 1 && (!op.isIn(resc.p) || !oq.isIn(resc.p));
    if ((df >= dttol && res2 <= 1) || kekis) {
      std::cout.precision(17);
      std::cout << op.b << op.s << op.x << op.y << op.z << op.a << op.c
                << std::endl;
      std::cout << oq.b << oq.s << oq.x << oq.y << oq.z << oq.a << oq.c
                << std::endl;
      std::cout << v::dot(op.x, resc.p) << " " << v::dot(op.y, resc.p) << " "
                << v::dot(op.z, resc.p) << std::endl;
      std::cout << v::dot(oq.x, resc.p) << " " << v::dot(oq.y, resc.p) << " "
                << v::dot(oq.z, resc.p) << std::endl;
      std::cout << resc.p << " " << resc.d << " " << resc.n << std::endl;

      std::cerr.precision(17);
      std::cerr << "FAIL!" << std::endl;
      std::cerr << "Iteration #" << spam << std::endl;
      std::cerr << "obbLinear result: " << res1 << std::endl;
      std::cerr << "stupid result: " << res2 << std::endl;
      std::cerr << "pp.p, pp.q: " << pp.p << pp.q << std::endl;
      std::cerr << "qp.p, qp.q: " << qp.p << qp.q << std::endl;
      std::cerr << "p.b,p.s,p.x,p.y,p.z,p.a,p.c: " << p.b << p.s << p.x << p.y
                << p.z << p.a << p.c << std::endl;
      std::cerr << "q.b,q.s,q.x,q.y,q.z,q.a,q.c: " << q.b << q.s << q.x << q.y
                << q.z << q.a << q.c << std::endl;
      std::cerr << "ps, qs: " << p.s << q.s << std::endl;
      std::cerr << "pv: " << pv << std::endl;
      std::cerr << "qv: " << qv << std::endl;
      return 1;
    }
  }
  return 0;
}
int main() {
  // return testIntervalInt();
  // return testLinear(1e6, 1e3, 1e5);
  return testRot(1e5, 1e4, 1e5);
}
