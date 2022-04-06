#ifndef ROLLER_PHYSICS_PRISM_HPP_
#define ROLLER_PHYSICS_PRISM_HPP_

#include "ccd_narrow.hpp"
#include "obb.hpp"
#include "phys_info.hpp"

namespace roller {

template <class Tag> struct Prism {

  v::DVec<3> s2; /// sidelengths divided by 2

  PhysInfo pi;

  bool doRender = true;
  Tag tag;
  bool meshed = false;

  Prism() : s2{0, 0, 0} {}

  /// s2 > 0, massi = 1 / mass
  Prism(const v::DVec<3> &s2, const v::DVec<3> &pos, double massi) : s2(s2) {
    pi.pose.p = pos;
    pi.massi = massi;
    if (pi.massi > 1e-12) {
      double xx = s2[1] * s2[1] + s2[2] * s2[2];
      double yy = s2[0] * s2[0] + s2[2] * s2[2];
      double zz = s2[0] * s2[0] + s2[1] * s2[1];
      double xxi = 1 / xx;
      double yyi = 1 / yy;
      double zzi = 1 / zz;
      double tmp = pi.mass = 1 / massi;
      tmp /= 3;
      pi.iner = {{xx * tmp, 0, 0}, {0, yy * tmp, 0}, {0, 0, zz * tmp}};
      tmp = massi * 3;
      pi.ineri = {{xxi * tmp, 0, 0}, {0, yyi * tmp, 0}, {0, 0, zzi * tmp}};
    }
  }

  Prism(const Prism &o) : s2(o.s2), pi(o.pi), doRender(o.doRender) {}
  Prism &operator=(const Prism &o) {
    s2 = o.s2;
    pi = o.pi;
    doRender = o.doRender;
    meshed = false;
    return *this;
  }

  void setPose(const Pose &pose) { pi.pose = pose; }
  PhysInfo getPhysInfo() const { return pi; }
  v::DVec<3> getSurfaceDetail() const { return {0.5, 0, 0}; }

  OBB getOBB() const {
    DMat3x3 rott = pi.pose.toRotationMatrix().transpose();
    v::DVec<3> mp = pi.pose.toWorldCoords(-s2);
    /*std::cout << "getOBB: " << mp
              << cross3(pi.getAuxInfo().omega, mp - pi.pose.p) << std::endl;
    std::cout << "riner: " << pi.getAuxInfo().riner.a << pi.getAuxInfo().riner.b
              << pi.getAuxInfo().riner.c << std::endl;*/
    return {mp, 2 * s2, rott.a, rott.b, rott.c};
  }

  // Prim interface functions
  Pose getPose() const { return pi.pose; }
  AABB getAABB(const ScrewM &sm) const { return getOBB().getAABB(sm); }
  ccd::Contact doCCD(double t1, double t2, const ScrewM &sm1,
                     const Prism &prim2, const ScrewM &sm2) const {
    Pose pose1 = getPose();
    pose1.p = pose1.toWorldCoords(-s2);
    Pose pose2 = prim2.getPose();
    pose2.p = pose2.toWorldCoords(-prim2.s2);
    if (t1 >= 1e-5) {
      apply2Pose(pose1, mult(t1, sm1));
      apply2Pose(pose2, mult(t2, sm2));
    }
    // XXX: Add tolerance as an option
    ccd::Contact toret = ccd::obbRot(pose1, 2 * s2, mult(t2 - t1, sm1), pose2,
                                     2 * prim2.s2, mult(t2 - t1, sm2), 1e-2);
    toret.t *= (t2 - t1);
    toret.t += t1;
    return toret;
  }

  template <class Fun> void exportAllTriangles(const SliceDirs &sd, Fun &&fun) {
    v::DVec<16> projMat = pi.pose.toProjMat(sd);
    fun(&projMat[0]);
    if (meshed) {
      fun(tag);
      return;
    }
    const Color c[] = {{0x00, 0x00, 0x00, 0xFF}, {0x22, 0x22, 0x22, 0xFF},
                       {0x44, 0x44, 0x44, 0xFF}, {0x66, 0x66, 0x66, 0xFF},
                       {0x88, 0x88, 0x88, 0xFF}, {0xAA, 0xAA, 0xAA, 0xFF},
                       {0xCC, 0xCC, 0xCC, 0xFF}, {0xFF, 0xFF, 0xFF, 0xFF}};

    float s2x = s2[0];
    float s2y = s2[1];
    float s2z = s2[2];
    float triangles[] = {
        // towards positive x-axis
        -s2x, +s2y, -s2z, c[2].f, // 1
        -s2x, -s2y, -s2z, c[0].f, // 2
        -s2x, -s2y, +s2z, c[1].f, // 3
        -s2x, -s2y, +s2z, c[1].f, // 1
        -s2x, +s2y, +s2z, c[3].f, // 2
        -s2x, +s2y, -s2z, c[2].f, // 3
        // towards negative x-axis
        +s2x, +s2y, -s2z, c[6].f, // 1
        +s2x, -s2y, +s2z, c[5].f, // 3
        +s2x, -s2y, -s2z, c[4].f, // 2
        +s2x, -s2y, +s2z, c[5].f, // 1
        +s2x, +s2y, -s2z, c[6].f, // 3
        +s2x, +s2y, +s2z, c[7].f, // 2
        // towards positive y-axis
        -s2x, -s2y, +s2z, c[1].f, // 1
        -s2x, -s2y, -s2z, c[0].f, // 2
        +s2x, -s2y, -s2z, c[4].f, // 3
        +s2x, -s2y, -s2z, c[4].f, // 1
        +s2x, -s2y, +s2z, c[5].f, // 2
        -s2x, -s2y, +s2z, c[1].f, // 3
        // towards negative y-axis
        -s2x, +s2y, +s2z, c[3].f, // 1
        +s2x, +s2y, -s2z, c[6].f, // 3
        -s2x, +s2y, -s2z, c[2].f, // 2
        +s2x, +s2y, -s2z, c[6].f, // 1
        -s2x, +s2y, +s2z, c[3].f, // 3
        +s2x, +s2y, +s2z, c[7].f, // 2
        // towards positive z-axis
        +s2x, -s2y, -s2z, c[4].f, // 1
        -s2x, -s2y, -s2z, c[0].f, // 2
        -s2x, +s2y, -s2z, c[2].f, // 3
        -s2x, +s2y, -s2z, c[2].f, // 1
        +s2x, +s2y, -s2z, c[6].f, // 2
        +s2x, -s2y, -s2z, c[4].f, // 3
        // towards negative z-axis
        +s2x, -s2y, +s2z, c[5].f, // 1
        -s2x, +s2y, +s2z, c[3].f, // 3
        -s2x, -s2y, +s2z, c[1].f, // 2
        -s2x, +s2y, +s2z, c[3].f, // 1
        +s2x, -s2y, +s2z, c[5].f, // 3
        +s2x, +s2y, +s2z, c[7].f, // 2
    };
    if (doRender) {
      fun(tag, triangles, sizeof(triangles) / sizeof(float));
    } else {
      fun(tag, 0, 0);
    }
    meshed = true;
    fun(tag);
  }
};

} // namespace roller

#endif // ROLLER_PHYSICS_PRISM_HPP_
