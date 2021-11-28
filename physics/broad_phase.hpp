#ifndef ROLLER_PHYSICS_BROAD_PHASE_HPP_
#define ROLLER_PHYSICS_BROAD_PHASE_HPP_

#include <priority_queue>

#include "aabb.hpp"
#include "pool.hpp"
#include "util.hpp"

namespace roller {

/// inspired by Bullet's b3DynamicBvh and Box2D's b2_dynamic_tree
class BroadPhaseAABBTree {
  struct Node {
    AABB aabb;
    int height;
    int parent;
    int child[2];
  };
  /// root is 0
  Pool<Node> p;

  /// this heuristic is from b3DynamicBvh
  double proximity(const AABB &a, const AABB &b) const {
    return v::sum(v::vabs((a.m[0] + a.m[1]) - (b.m[0] + b.m[1])));
  }

  static int imax(int a, int b) { return a > b ? a : b; }
  void writeAABB(int node) {
    p[node].aabb = p[p[node].left].aabb.combine(p[p[node].right].aabb);
  }
  /// does not fix the AABBs in the path to highest subtree; helper for insert
  void balance(int node) {
    Node &n = p[node];
    int nl = n.left;
    int nr = n.right;
    int hdf = p[nl].height - p[nr].height;
    // "Dynamic AABB Trees - GDC 2019"
    if (hdf >= 2) {
      int nll = p[nl].left;
      int nlr = p[nl].right;
      if (p[nll].height > p[nlr].height) {
        // swap nll with nr
        n.right = nll;
        p[nll].parent = node;
        p[nl].left = nr;
        p[nr].parent = nl;
        p[nl].height--;
        writeAABB(nl);
      } else {
        // swap nlr with nr
        n.right = nlr;
        p[nlr].parent = node;
        p[nl].right = nr;
        p[nr].parent = nl;
        p[nl].height--;
        writeAABB(nl);
      }
    } else if (hdf <= 2) {
      int nrl = p[nr].left;
      int nrr = p[nr].right;
      if (p[nrl].height > p[nrr].height) {
        // swap nrl with nl
        n.left = nrl;
        p[nrl].parent = node;
        p[nr].left = nl;
        p[nl].parent = nr;
        p[nr].height--;
        writeAABB(nr);
      } else {
        // swap nrr with nl
        n.left = nrr;
        p[nrr].parent = node;
        p[nr].right = nl;
        p[nl].parent = nr;
        p[nr].height--;
        writeAABB(nr);
      }
    }
    n.height = imax(p[n.left].height, p[n.right].height) + 1;
  }

  /// helper for insert
  bool properifyAABB(int node, const AABB &aabb) {
    if (p[node].aabb.contains(aabb)) {
      return false;
    }
    p[node].aabb = p[node].aabb.combine(aabb);
    return true;
  }

public:
  BroadPhaseAABBTree() {
    p.mkNew(); // should return 0
    p[0].height = -1;
    p[0].parent = 0;
    p[0].child[0] = 0;
    p[0].child[1] = 0;
  }

  void insert(const AABB &aabb) {
    // empty tree
    if (p[0].height < 0) {
      p[0] = {aabb, 0, 0, {0, 0}};
      return;
    }
    int sibling = 0;
    int lastStep = 0;
    while (p[sibling].height > 0) {
      sibling =
          p[sibling]
              .child[lastStep = (proximity(aabb, p[sibling].child[1].aabb) >
                                 proximity(aabb, p[sibling].child[0].aabb))];
    }
    int gpar = p[sibling].parent;
    int parent = p.mkNew();
    if (sibling) {
      p[gpar].child[lastStep] = parent;
    } else {
      sibling = parent;
      parent = 0;
      p[sibling] = p[0];
    }
    int leaf = p.mkNew();
    p[leaf] = {aabb, 0, parent, {0, 0}};
    p[sibling].parent = parent;
    p[parent].aabb = aabb.combine(p[sibling].aabb);
    p[parent].height = 1;
    p[parent] = gpar;
    p[parent].child[0] = sibling;
    p[parent].child[1] = leaf;
    int node = parent;
    while (node) {
      node = p[node].parent;
      balance(node);
    }
    node = parent;
    while (node) {
      node = p[node].parent;
      if (!properifyAABB(node, aabb)) {
        break;
      }
    }
  }
};

} // namespace roller

#endif // ROLLER_PHYSICS_BROAD_PHASE_HPP_
