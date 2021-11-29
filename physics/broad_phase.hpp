#ifndef ROLLER_PHYSICS_BROAD_PHASE_HPP_
#define ROLLER_PHYSICS_BROAD_PHASE_HPP_

#include <queue>
#include <unordered_set>
#include <vector>

#include "aabb.hpp"
#include "pool.hpp"
#include "util.hpp"

namespace roller {

/// inspired by Bullet's b3DynamicBvh and Box2D's b2_dynamic_tree
class BroadPhaseAABBTree {

  struct Node {
    AABB aabb;
    int primi; /// index of the primitive this is AABB around
    int height;
    int parent;
    int child[2];
  };
  Pool<Node> p;
  int root = -1;

  /// this heuristic is from b3DynamicBvh
  double proximity(const AABB &a, const AABB &b) const {
    return v::sum(v::vabs((a.m[0] + a.m[1]) - (b.m[0] + b.m[1])));
  }
  void writeAABB(int node) {
    p[node].aabb = p[p[node].child[0]].aabb.combine(p[p[node].child[1]].aabb);
  }
  static int imax(int a, int b) { return a > b ? a : b; }
  void writeHeight(int node) {
    p[node].height =
        imax(p[p[node].child[0]].height, p[p[node].child[1]].height) + 1;
  }
  /// does not fix the AABBs in the path to highest subtree; helper for insert
  void balance(int node) {
    Node &n = p[node];
    int nl = n.child[0];
    int nr = n.child[1];
    int hdf = p[nl].height - p[nr].height;
    // "Dynamic AABB Trees - GDC 2019"
    if (hdf >= 2) {
      int nll = p[nl].child[0];
      int nlr = p[nl].child[1];
      if (p[nll].height > p[nlr].height) {
        // swap nll with nr
        n.child[1] = nll;
        p[nll].parent = node;
        p[nl].child[0] = nr;
        p[nr].parent = nl;
        writeAABB(nl);
        writeHeight(nl);
      } else {
        // swap nlr with nr
        n.child[1] = nlr;
        p[nlr].parent = node;
        p[nl].child[1] = nr;
        p[nr].parent = nl;
        writeAABB(nl);
        writeHeight(nl);
      }
    } else if (hdf <= -2) {
      int nrl = p[nr].child[0];
      int nrr = p[nr].child[1];
      if (p[nrl].height > p[nrr].height) {
        // swap nrl with nl
        n.child[0] = nrl;
        p[nrl].parent = node;
        p[nr].child[0] = nl;
        p[nl].parent = nr;
        writeAABB(nr);
        writeHeight(nr);
      } else {
        // swap nrr with nl
        n.child[0] = nrr;
        p[nrr].parent = node;
        p[nr].child[1] = nl;
        p[nl].parent = nr;
        writeAABB(nr);
        writeHeight(nr);
      }
    }
    writeHeight(node);
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
  /// returns node index
  int insert(int primi, const AABB &aabb) {
    // empty tree
    if (root < 0) {
      root = p.mkNew();
      p[root] = {aabb, primi, 0, root, {0, 0}};
      return root;
    }
    int sibling = root;
    int lastStep = 0;
    while (p[sibling].height > 0) {
      sibling =
          p[sibling]
              .child[lastStep = (proximity(aabb, p[p[sibling].child[1]].aabb) >
                                 proximity(aabb, p[p[sibling].child[0]].aabb))];
    }
    int gpar = p[sibling].parent;
    int parent = p.mkNew();
    if (sibling != root) {
      p[gpar].child[lastStep] = parent;
    } else {
      root = parent;
    }
    int leaf = p.mkNew();
    p[leaf] = {aabb, primi, 0, parent, {0, 0}};
    p[sibling].parent = parent;
    p[parent].aabb = aabb.combine(p[sibling].aabb);
    p[parent].primi = -1;
    p[parent].height = 1;
    p[parent].parent = gpar;
    p[parent].child[0] = sibling;
    p[parent].child[1] = leaf;
    int node = parent;
    while (node != root) {
      node = p[node].parent;
      balance(node);
    }
    node = parent;
    while (node != root) {
      node = p[node].parent;
      if (!properifyAABB(node, aabb)) {
        break;
      }
    }
    return leaf;
  }

  void remove(int leaf) {
    int parent = p[leaf].parent;
    p.rem(leaf);
    if (leaf == root) {
      root = -1;
      return;
    }
    int sibling;
    if (leaf == p[parent].child[0]) {
      sibling = p[parent].child[1];
    } else {
      sibling = p[parent].child[0];
    }
    int gpar = p[parent].parent;
    p.rem(parent);
    if (parent != root) {
      if (parent == p[gpar].child[0]) {
        p[gpar].child[0] = sibling;
      } else {
        p[gpar].child[1] = sibling;
      }
      p[sibling].parent = gpar;
      int node = sibling;
      while (node != root) {
        node = p[node].parent;
        balance(node);
      }
      node = sibling;
      while (node != root) {
        node = p[node].parent;
        writeAABB(node);
      }
    } else {
      root = sibling;
    }
  }

  template <class Fun> int update(int leaf, const AABB &aabb, Fun &&fun) {
    if (p[leaf].aabb.contains(aabb)) {
      return leaf;
    }
    int primi = p[leaf].primi;
    remove(leaf);
    return insert(primi, fun(aabb));
  }

  /// calls fun(a, b) where a,b are primi of intersecting nodes
  template <class Fun> void exportInts(Fun &&fun) {
    if (root < 0) {
      return;
    }
    struct Entry {
      int a, b;
    };
    std::vector<Entry> stack;
    stack.push_back({root, root});
    while (stack.size()) {
      Entry e = stack.back();
      stack.pop_back();
      if (e.a == e.b) {
        if (p[e.a].primi >= 0) {
          continue;
        }
        stack.push_back({p[e.a].child[0], p[e.a].child[0]});
        stack.push_back({p[e.a].child[0], p[e.a].child[1]});
        stack.push_back({p[e.a].child[1], p[e.a].child[1]});
      } else {
        if (!p[e.a].aabb.intersects(p[e.b].aabb)) {
          continue;
        }
        if (p[e.a].primi >= 0) {
          if (p[e.b].primi >= 0) {
            fun(p[e.a].primi, p[e.b].primi);
          } else {
            stack.push_back({e.a, p[e.b].child[0]});
            stack.push_back({e.a, p[e.b].child[1]});
          }
        } else {
          if (p[e.b].primi >= 0) {
            stack.push_back({p[e.a].child[0], e.b});
            stack.push_back({p[e.a].child[1], e.b});
          } else {
            stack.push_back({p[e.a].child[0], p[e.b].child[0]});
            stack.push_back({p[e.a].child[0], p[e.b].child[1]});
            stack.push_back({p[e.a].child[1], p[e.b].child[0]});
            stack.push_back({p[e.a].child[1], p[e.b].child[1]});
          }
        }
      }
    }
  }

  /// returns false if it finds anything wrong with the structure
  bool testStructure(int printTree = 0) const {
    if (root < 0) {
      return true;
    }
    std::unordered_set<int> inds;
    struct Entry {
      int node, indent;
    };
    std::vector<Entry> stack;
    stack.push_back({root, 0});
    bool fail = true;
    while (stack.size()) {
      Entry e = stack.back();
      stack.pop_back();
      int node = e.node;
      if (printTree) {
        std::cout << std::string(e.indent, ' ') << node << "," << p[node].primi
                  << "," << p[node].height << std::endl;
        if (printTree > 1) {
          p[node].aabb.print();
        }
      }
      if (inds.count(node)) {
        std::cout << node << std::endl;
        std::cout << "cyclic" << std::endl;
        fail = false;
      }
      if (p[node].primi >= 0) {
        if (p[node].height) {
          std::cout << node << std::endl;
          std::cout << "leaf has height" << std::endl;
          fail = false;
        }
        continue;
      }
      int nl = p[node].child[0];
      int nr = p[node].child[1];
      if (p[node].height != imax(p[nl].height, p[nr].height) + 1) {
        std::cout << node << std::endl;
        std::cout << "height incorrect" << std::endl;
        fail = false;
      }
      if (p[nl].parent != node || p[nr].parent != node) {
        std::cout << node << std::endl;
        std::cout << "parent incorrect" << std::endl;
        fail = false;
      }
      if (!p[node].aabb.contains(p[nl].aabb) ||
          !p[node].aabb.contains(p[nr].aabb)) {
        std::cout << node << std::endl;
        std::cout << "aabb not bounding" << std::endl;
        fail = false;
      }
      stack.push_back({nl, e.indent + 8});
      stack.push_back({nr, e.indent + 8});
    }
    return fail;
  }
};

} // namespace roller

#endif // ROLLER_PHYSICS_BROAD_PHASE_HPP_
