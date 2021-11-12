#ifndef ROLLER_UNION_FIND_HPP_
#define ROLLER_UNION_FIND_HPP_

#include <vector>

class UnionFind {

  struct Node {
    int p;
    int sz;
  };
  std::vector<Node> p;

public:
  UnionFind() {}
  explicit UnionFind(std::size_t sz) : p(sz) {
    for (std::size_t i = 0; i < sz; i++) {
      p[i] = {i, 1};
    }
  }

  /// find representative index of the set i is in
  int find(int i) {
    // path halving
    while (p[i].p != i) {
      i = p[i].p = p[p[i].p].p;
    }
    return i;
  }
  /// returns true if a union was performed
  bool merge(int i, int j) {
    int pi = find(i);
    int pj = find(j);
    if (pi == pj) {
      return false;
    }
    if (p[pi].sz < p[pj].sz) {
      int tmp = pi;
      pi = pj;
      pj = tmp;
    }
    p[pj].p = pi;
    p[pi].sz += p[pj].sz;
    return true;
  }

  /// to identify singletons, of course
  int getWholeSz(int i) const { return p[find(i)].sz; }
  /// if you already have the representative index for some reason
  int getIndSz(int i) const { return p[i].sz; }
};

#endif // ROLLER_UNION_FIND_HPP_

