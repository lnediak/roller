#ifndef ROLLER_PHYSICS_PRISM_HPP_
#define ROLLER_PHYSICS_PRISM_HPP_

#include "phys_info.hpp"

namespace roller {

template <class Tag> struct Prism {

  v::DVec<3> s2; /// sidelengths divided by 2

  PhysInfo pi;

  Tag tag;
  bool meshed = false;

  /// ds2 > 0, massi = 1 / mass
  Prism(const v::DVec<3> &ds2, double massi) : s2(ds2) {
    pi.massi = massi;
    double xx = 1 / (s2[1] * s2[1] + s2[2] * s2[2]);
    double yy = 1 / (s2[0] * s2[0] + s2[2] * s2[2]);
    double zz = 1 / (s2[0] * s2[0] + s2[1] * s2[1]);
    massi *= 3; // because s2 is already divided by 2
    pi.ineri = {{xx * massi, 0, 0}, {0, yy * massi, 0}, {0, 0, zz * massi}};
  }

  OBB getOBB() const {
    DMat3x3 rott = pi.pose.toRotationMatrix().transpose();
    v::DVec<3> mp = pi.toWorldCoords(-s2);
    return {mp, 2 * s2[0] * rott.a, 2 * s2[1] * rott.b, 2 * s2[2] * rott.c};
  }
  AABB getAABB() const { return getOBB().getAABB(); }
  PhysInfo getPhysInfo() const { return pi; }
  void setPhysInfo(const PhysInfo &pi) { this->pi = pi; }

  Contact getContacts(const Prism &q) const {
    OBBIntersector inter(getOBB(), q.getOBB());
    return inter.getInts();
  }

  template <class Fun> void exportAllTriangles(const SliceDirs &sd, Fun &&fun) {
    v::DVec<16> projMat = pi.pose.toProjMat(sd);
    fun(&projMat[0]);
    if (meshed) {
      fun(tag);
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
    fun(tag, triangles, sizeof(triangles) / sizeof(float));
    meshed = true;
  }
};

} // namespace roller

#endif // ROLLER_PHYSICS_PRISM_HPP_

