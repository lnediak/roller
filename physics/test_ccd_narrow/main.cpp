// DISCLAIMER: This is by no means a thorough unit test.

#include <iostream>
#include <random>

#include "ccd_narrow.hpp"

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
  Plane planes[] = {{a.x, a.a[0]}, {a.x, a.c[0]}, {a.y, a.a[1]}, {a.y, a.c[1]},
                    {a.z, a.a[2]}, {a.z, a.a[2]}, {b.x, b.a[0]}, {b.x, b.c[0]},
                    {b.y, b.a[1]}, {b.y, b.c[1]}, {b.z, b.a[2]}, {b.z, b.a[2]}};
  OBB fa = a, fb = b;
  fa.a -= 1e-15;
  fa.c += 1e-15;
  fb.a -= 1e-15;
  fb.c += 1e-15;
  for (int i = 0; i < 10; i++) {
    for (int j = i + 1; j < 11; j++) {
      for (int k = j + 1; k < 12; k++) {
        v::DVec<3> pp = planes[i].get(planes[j], planes[k]);
        if (fa.isIn(pp) && fb.isIn(pp)) {
          return true;
        }
      }
    }
  }
  return false;
}

void apply(OBB &p, v::DVec<3> trans, v::DVec<3> omega, v::DVec<3> c) {
  v::DVec<4> orot = getRotQuaternion(omega);
  // TODO: REMOVE THIS
  if (orot[0] != 1 || orot[1] || orot[2] || orot[3]) {
    std::cerr << "FAILURE    >LLOL" << std::endl;
    std::terminate();
  }
  if (v::norm2(p.b - c - applyQuaternion(p.b - c, orot)) >= 1e-20) {
    std::cerr << "FAILURE    >LLOL" << std::endl;
    std::terminate();
  }

  p.b = trans + c + applyQuaternion(p.b - c, orot);
  p.x = applyQuaternion(p.x, orot);
  p.y = applyQuaternion(p.y, orot);
  p.z = applyQuaternion(p.z, orot);
  p.properify();
}

double getTime(const OBB &p, v::DVec<3> pv, v::DVec<3> omega1, v::DVec<3> pc,
               const OBB &q, v::DVec<3> qv, v::DVec<3> omega2, v::DVec<3> qc) {
  for (double t = 0; t < 1; t += 1e-4) {
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

int testRot() {
  std::mt19937 mtrand(1);
  std::uniform_real_distribution<> distro(-100, 100);
  std::uniform_real_distribution<> smol(-0.05, 0.05);
  std::uniform_real_distribution<> pdistro(0, 100);
  for (int spam = 0; spam < 100000; spam++) {
    if (spam % 1000 == 0) {
      std::cout << "Iteration #" << spam << std::endl;
    }
    Pose pp;
    pp.q = normalizeQuaternion(
        {distro(mtrand), distro(mtrand), distro(mtrand), distro(mtrand)});
    Pose qp;
    qp.q = normalizeQuaternion(
        {distro(mtrand), distro(mtrand), distro(mtrand), distro(mtrand)});
    v::DVec<3> omega1{smol(mtrand), smol(mtrand), smol(mtrand)};
    v::DVec<3> omega2{smol(mtrand), smol(mtrand), smol(mtrand)};
    v::DVec<3> pc{distro(mtrand), distro(mtrand), distro(mtrand)};
    v::DVec<3> qc{distro(mtrand), distro(mtrand), distro(mtrand)};

    CCDRotOBBIntersector inn(pp, {1, 1, 1}, omega1, pc, qp, {1, 1, 1}, omega2,
                             qc, 1e-5);
    for (int spam2 = 0; spam2 < 100; spam2++) {
      pp.p = {distro(mtrand), distro(mtrand), distro(mtrand)};
      qp.p = {distro(mtrand), distro(mtrand), distro(mtrand)};
      v::DVec<3> pv{distro(mtrand), distro(mtrand), distro(mtrand)};
      v::DVec<3> qv{distro(mtrand), distro(mtrand), distro(mtrand)};
      OBB p = CCDRotOBBIntersector::poseToOBB(
          pp, {pdistro(mtrand), pdistro(mtrand), pdistro(mtrand)});
      OBB q = CCDRotOBBIntersector::poseToOBB(
          qp, {pdistro(mtrand), pdistro(mtrand), pdistro(mtrand)});
      inn.updateOBBs(pp.p, p.s, qp.p, q.s);
      Contact resc = inn.getInt(pv, qv);
      double res1 = resc.t;
      if (res1 > 1) {
        res1 = 1.0001;
      }
      double res2 = getTime(p, pv, omega1, pc, q, qv, omega2, qc);
      double df = res1 > res2 ? res1 - res2 : res2 - res1;
      if (df >= 0.05) {

        // tmp code:
        OBB op = p;
        apply(op, res1 * pv, res1 * omega1, pc);
        OBB oq = q;
        apply(oq, res1 * qv, res1 * omega2, qc);
        std::cout << op.b << op.s << op.x << op.y << op.z << op.a << op.c
                  << std::endl;
        std::cout << oq.b << oq.s << oq.x << oq.y << oq.z << oq.a << oq.c
                  << std::endl;
        std::cout << v::dot(op.x, resc.p) << " " << v::dot(op.y, resc.p) << " "
                  << v::dot(op.z, resc.p) << std::endl;
        std::cout << v::dot(oq.x, resc.p) << " " << v::dot(oq.y, resc.p) << " "
                  << v::dot(oq.z, resc.p) << std::endl;
        std::cout << resc.p << std::endl;

        std::cerr << "FAIL!" << std::endl;
        std::cerr << "Iteration #" << spam << std::endl;
        std::cerr << "Sub-Iteration #" << spam2 << std::endl;
        std::cerr << "CCDRotOBBIntersector result: " << res1 << std::endl;
        std::cerr << "stupid result: " << res2 << std::endl;
        std::cerr << "pp.p, pp.q: " << pp.p << pp.q << std::endl;
        std::cerr << "qp.p, qp.q: " << qp.p << qp.q << std::endl;
        std::cerr << "ps, qs: " << p.s << q.s << std::endl;
        std::cerr << "pv, pc, omega1: " << pv << pc << omega1 << std::endl;
        std::cerr << "qv, qc, omega2: " << qv << qc << omega2 << std::endl;
        return 1;
      }
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
    v::DVec<3> res3 = CCDOBBIntersector::intervalInt(v2 - v1, a - d, b - c);
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
int testLinear() {
  std::mt19937 mtrand(1);
  std::uniform_real_distribution<> distro(-100, 100);
  std::uniform_real_distribution<> smol(-0.05, 0.05);
  std::uniform_real_distribution<> pdistro(0, 100);
  for (int spam = 0; spam < 100000; spam++) {
    if (spam % 1000 == 0) {
      std::cout << "Iteration #" << spam << std::endl;
    }
    Pose pp;
    pp.q = normalizeQuaternion(
        {distro(mtrand), distro(mtrand), distro(mtrand), distro(mtrand)});
    Pose qp;
    qp.q = normalizeQuaternion(
        {distro(mtrand), distro(mtrand), distro(mtrand), distro(mtrand)});
    OBB p = CCDRotOBBIntersector::poseToOBB(pp, {1, 1, 1});
    OBB q = CCDRotOBBIntersector::poseToOBB(qp, {1, 1, 1});

    CCDOBBIntersector inn(p, q);
    for (int spam2 = 0; spam2 < 100; spam2++) {
      pp.p = {distro(mtrand), distro(mtrand), distro(mtrand)};
      qp.p = {distro(mtrand), distro(mtrand), distro(mtrand)};
      v::DVec<3> pv{distro(mtrand), distro(mtrand), distro(mtrand)};
      v::DVec<3> qv{distro(mtrand), distro(mtrand), distro(mtrand)};
      OBB p = CCDRotOBBIntersector::poseToOBB(
          pp, {pdistro(mtrand), pdistro(mtrand), pdistro(mtrand)});
      OBB q = CCDRotOBBIntersector::poseToOBB(
          qp, {pdistro(mtrand), pdistro(mtrand), pdistro(mtrand)});
      inn.setOBBs(p, q);
      Contact resc = inn.getInt(pv, qv);
      double res1 = resc.t;
      if (res1 > 1) {
        res1 = 1.0001;
      }
      v::DVec<3> zero{0, 0, 0};
      double res2 = getTime(p, pv, zero, zero, q, qv, zero, zero);
      double df = res1 > res2 ? res1 - res2 : res2 - res1;
      if (df >= 0.05) {

        // tmp code:
        OBB op = p;
        apply(op, res1 * pv, zero, zero);
        OBB oq = q;
        apply(oq, res1 * qv, zero, zero);
        std::cout << op.b << op.s << op.x << op.y << op.z << op.a << op.c
                  << std::endl;
        std::cout << oq.b << oq.s << oq.x << oq.y << oq.z << oq.a << oq.c
                  << std::endl;
        std::cout << v::dot(op.x, resc.p) << " " << v::dot(op.y, resc.p) << " "
                  << v::dot(op.z, resc.p) << std::endl;
        std::cout << v::dot(oq.x, resc.p) << " " << v::dot(oq.y, resc.p) << " "
                  << v::dot(oq.z, resc.p) << std::endl;
        std::cout << resc.p << " " << resc.d << " " << resc.n << std::endl;
        // end tmp code

        std::cerr << "FAIL!" << std::endl;
        std::cerr << "Iteration #" << spam << std::endl;
        std::cerr << "Sub-Iteration #" << spam2 << std::endl;
        std::cerr << "CCDRotOBBIntersector result: " << res1 << std::endl;
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
  }
  return 0;
}
int main() {
  // return testIntervalInt();
  return testLinear();
  // return testRot();
}

