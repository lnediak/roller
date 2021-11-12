#ifndef ROLLER_PHYSICS_BROAD_PHASE_HPP_
#define ROLLER_PHYSICS_BROAD_PHASE_HPP_

#include <algorithm>
#include <unordered_set>
#include <vector>

#include "aabb.hpp"
#include "util.hpp"

namespace roller {

/**
  W template: Needs the follows
  typedef ? iterator;  // just needs to support ++ and !=
  typedef ? iter_hash; // a hash for iterator
  typedef ? iter_eq;   // a equality functor for iterator
  iterator begin() const ;
  iterator end() const ;
  AABB getAABB(iterator) const ;
  Note: getAABB should be fast

  IgnFun template (used in class body): Needs the follows
  bool operator()(W::iterator a, W::iterator b) const;
  which returns true if collision between a and b should be ignored
*/
template <class W> class BroadPhaseAABB {

  W w; /// access point for AABBs
  typedef rem_cvr<W> clnW;
  typedef typename clnW::iterator Witer;
  typedef typename clnW::iter_hash WiterHash;
  typedef typename clnW::iter_eq WiterEq;

  // AABB stuff
  struct Entry {
    Witer i;
    bool b; // is this the max of the AABB along the axis in question?
  };
  std::vector<Entry> xsort, ysort, zsort;
  /// pairs of intersecting AABBs
  std::unordered_set<UnorderedPair<Witer>, UnorderedPairHash<Witer, WiterHash>,
                     UnorderedPairEq<Witer, WiterEq>>
      aabbInts;

  template <class IgnFun> void initAABBStuff(IgnFun &&ignFun) {
    aabbInts.clear();
    xsort.clear();
    ysort.clear();
    zsort.clear();
    int nObjs = 0;
    for (Witer it = w.begin(), ite = w.end(); it != ite; ++it) {
      xsort.push_back({it, false});
      ysort.push_back({it, false});
      zsort.push_back({it, false});
      xsort.push_back({it, true});
      ysort.push_back({it, true});
      zsort.push_back({it, true});
      nObjs++;
    }
    if (!nObjs) {
      return;
    }
    int xyz;
    auto cmp = [this, &xyz](const Entry &a, const Entry &b) -> bool {
      return w.getAABB(a.i).m[a.b][xyz] < w.getAABB(b.i).m[b.b][xyz];
    };
    // sort...
    xyz = 0;
    std::sort(xsort.begin(), xsort.end(), cmp);
    xyz = 1;
    std::sort(ysort.begin(), ysort.end(), cmp);
    xyz = 2;
    std::sort(zsort.begin(), zsort.end(), cmp);
    // and sweep! (not most efficient lol)
    std::size_t bloat = 3 * nObjs / 2;
    decltype(aabbInts) tmpx(bloat), tmpy(bloat);
    aabbInts.reserve(bloat);
    std::unordered_set<Witer, WiterHash, WiterEq> active(bloat);
    for (const Entry &e : xsort) {
      if (e.b) {
        active.erase(e.i); // note: assuming no degenerate AABBs
      } else {
        for (const Witer &fi : active) {
          tmpx.insert({e.i, fi});
        }
        active.insert(e.i);
      }
    }
    // active should be empty right now
    for (const Entry &e : ysort) {
      if (e.b) {
        active.erase(e.i);
      } else {
        for (const Witer &fi : active) {
          if (tmpx.count({e.i, fi})) {
            tmpy.insert({e.i, fi});
          }
        }
        active.insert(e.i);
      }
    }
    for (const Entry &e : zsort) {
      if (e.b) {
        active.erase(e.i);
      } else {
        for (const Witer &fi : active) {
          if (!ignFun(e.i, fi) && tmpy.count({e.i, fi})) {
            aabbInts.insert({e.i, fi});
          }
        }
        active.insert(e.i);
      }
    }
  }

public:
  struct IgnNoOp {
    bool operator()(Witer, Witer) const { return false; }
  };
  BroadPhaseAABB(W &&w) : w(w) { initAABBStuff(IgnNoOp()); }
  template <class IgnFun = IgnNoOp &&>
  BroadPhaseAABB(W &&w, IgnFun &&ignFun = IgnNoOp()) : w(w) {
    initAABBStuff(std::forward<IgnFun>(ignFun));
  }

private:
  template <class IgnFun, class Fun>
  void updateAABBStuff0(std::vector<Entry> &dsort, std::size_t d,
                        IgnFun &&ignFun, Fun &&callback) {
    if (dsort.empty()) {
      return;
    }
    Entry rpe = dsort[0]; // real previous entry
    double rp = w.getAABB(rpe.i).m[rpe.b][d];
    // insertion sort
    for (std::size_t ind = 1, sz = dsort.size(); ind < sz; ind++) {
      std::size_t cind = ind;
      Entry tce = dsort[cind]; // temporary current entry
      Entry tpe = rpe;         // temporary previous entry
      double tc = w.getAABB(tce.i).m[tce.b][d];
      double tp = rp;
      if (tp <= tc) {
        rpe = tce;
        rp = tc;
        continue;
      }
      rpe = tpe;
      rp = tp;

      do {
        dsort[cind] = tpe;
        dsort[--cind] = tce;
        if (tpe.b != tce.b && !ignFun(tpe.i, tce.i)) {
          if (w.getAABB(tpe.i).intersects(w.getAABB(tce.i))) {
            if (aabbInts.insert({tpe.i, tce.i}).second) {
              callback(tpe.i, tce.i, true);
            }
          } else {
            if (aabbInts.erase({tpe.i, tce.i})) {
              callback(tpe.i, tce.i, false);
            }
          }
        }
        if (!cind) {
          break;
        }
        tpe = dsort[cind - 1];
        tp = w.getAABB(tpe.i).m[tpe.b][d];
      } while (tp > tc);
    }
  }

public:
  struct NoOp {
    void operator()(Witer, Witer, bool) const {}
  };
  /// TODO: callback doesn't work if changedW
  /*
    changedW: refers to whether w.begin() through w.end() have changed
    callback: calls callback(a, b, isAdd), where isAdd is true if {a, b} is a
    new collision, and false if {a, b} is now not a collision
  */
  template <class IgnFun = IgnNoOp &&, class Fun = NoOp &&>
  void updateAABBStuff(IgnFun &&ignFun = IgnNoOp(), bool changedW = true,
                       Fun &&callback = NoOp()) {
    if (changedW) {
      // TODO: do a faster incremental update
      initAABBStuff(std::forward<IgnFun>(ignFun));
      return;
    }
    updateAABBStuff0(xsort, 0, std::forward<IgnFun>(ignFun),
                     std::forward<Fun>(callback));
    updateAABBStuff0(ysort, 1, std::forward<IgnFun>(ignFun),
                     std::forward<Fun>(callback));
    updateAABBStuff0(zsort, 2, std::forward<IgnFun>(ignFun),
                     std::forward<Fun>(callback));
  }

  /// calls fun(i, j) for all i,j indices of intersecting object AABBs
  template <class Fun> void exportInts(Fun &&fun) {
    for (auto &pair : aabbInts) {
      fun(pair.a, pair.b);
    }
  }
};

} // namespace roller

#endif // ROLLER_PHYSICS_BROAD_PHASE_HPP_
