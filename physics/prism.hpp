#ifndef ROLLER_PHYSICS_PRISM_HPP_
#define ROLLER_PHYSICS_PRISM_HPP_

#include "phys_info.hpp"
#include "vector.hpp"

namespace roller {

template <class Tag> struct Prism {

  v::FVec<3> s2; /// sidelengths divided by 2

  PhysInfo pi;

  Tag tag;
  bool meshed = false;

  Prism(const v::DVec<3> &ds2) {
    s2[0] = ds2[0];
    s2[1] = ds2[1];
    s2[2] = ds2[2];
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

    float triangles[] = {
        // towards positive x-axis
        -s2[0], +s2[1], -s2[2], c[2].f, // 1
        -s2[0], -s2[1], -s2[2], c[0].f, // 2
        -s2[0], -s2[1], +s2[2], c[1].f, // 3
        -s2[0], -s2[1], +s2[2], c[1].f, // 1
        -s2[0], +s2[1], +s2[2], c[3].f, // 2
        -s2[0], +s2[1], -s2[2], c[2].f, // 3
        // towards negative x-axis
        +s2[0], +s2[1], -s2[2], c[6].f, // 1
        +s2[0], -s2[1], +s2[2], c[5].f, // 3
        +s2[0], -s2[1], -s2[2], c[4].f, // 2
        +s2[0], -s2[1], +s2[2], c[5].f, // 1
        +s2[0], +s2[1], -s2[2], c[6].f, // 3
        +s2[0], +s2[1], +s2[2], c[7].f, // 2
        // towards positive y-axis
        -s2[0], -s2[1], +s2[2], c[1].f, // 1
        -s2[0], -s2[1], -s2[2], c[0].f, // 2
        +s2[0], -s2[1], -s2[2], c[4].f, // 3
        +s2[0], -s2[1], -s2[2], c[4].f, // 1
        +s2[0], -s2[1], +s2[2], c[5].f, // 2
        -s2[0], -s2[1], +s2[2], c[1].f, // 3
        // towards negative y-axis
        -s2[0], +s2[1], +s2[2], c[3].f, // 1
        +s2[0], +s2[1], -s2[2], c[6].f, // 3
        -s2[0], +s2[1], -s2[2], c[2].f, // 2
        +s2[0], +s2[1], -s2[2], c[6].f, // 1
        -s2[0], +s2[1], +s2[2], c[3].f, // 3
        +s2[0], +s2[1], +s2[2], c[7].f, // 2
        // towards positive z-axis
        +s2[0], -s2[1], -s2[2], c[4].f, // 1
        -s2[0], -s2[1], -s2[2], c[0].f, // 2
        -s2[0], +s2[1], -s2[2], c[2].f, // 3
        -s2[0], +s2[1], -s2[2], c[2].f, // 1
        +s2[0], +s2[1], -s2[2], c[6].f, // 2
        +s2[0], -s2[1], -s2[2], c[4].f, // 3
        // towards negative z-axis
        +s2[0], -s2[1], +s2[2], c[5].f, // 1
        -s2[0], +s2[1], +s2[2], c[3].f, // 3
        -s2[0], -s2[1], +s2[2], c[1].f, // 2
        -s2[0], +s2[1], +s2[2], c[3].f, // 1
        +s2[0], -s2[1], +s2[2], c[5].f, // 3
        +s2[0], +s2[1], +s2[2], c[7].f, // 2
    };
    fun(tag, triangles, sizeof(triangles) / sizeof(float));
    meshed = true;
  }
};

} // namespace roller

#endif // ROLLER_PHYSICS_PRISM_HPP_

