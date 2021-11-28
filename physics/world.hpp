#ifndef ROLLER_PHYSICS_WORLD_HPP_
#define ROLLER_PHYSICS_WORLD_HPP_

#include <iterator>
#include <memory>
#include <type_traits>
#include <types>
#include <utility>

#include "aabb.hpp"
#include "phys_info.hpp"
#include "pool.hpp"

namespace roller {

/// inserts ni after i
template <class Indexable> void dllAddAfter(Indexable &p, int i, int ni) {
  p[ni].previ = i;
  p[ni].nexti = p[i].nexti;
  p[p[i].nexti].previ = ni;
  p[i].nexti = ni;
}
template <class Indexable> void dllRemove(Indexable &p, int i) {
  p[p[i].previ].nexti = p[i].nexti;
  p[p[i].nexti].previ = p[i].previ;
}
template <class Indexable, class Fun> struct DllIterator {
  typedef std::ptrdiff_t difference_type;
  typedef typename std::decay<Fun>::type::return_type value_type;
  typedef value_type *pointer;
  typedef value_type &reference;
  typedef std::bidirectional_iterator_tag iterator_category;

  Indexable *p;
  int i;
  Fun fun;

  // say hi to DefaultConstructible
  DllIterator() : p(nullptr), i(-1) {}
  DllIterator(Indexable *p, int i, Fun &&fun) : p(p), i(i), fun(fun) {}

  reference operator*() { return fun((*p)[i]); }
  pointer operator->() { return &fun((*p)[i]); }
  DllIterator<Indexable> &operator++() {
    i = p[i].nexti;
    return *this;
  }
  DllIterator<Indexable> operator++(int) {
    const DllIterator<Indexable> &toret = *this;
    operator++();
    return toret;
  }
  DllIterator<Indexable> &operator--() {
    i = p[i].previ;
    return *this;
  }
  DllIterator<Indexable> operator--(int) {
    const DllIterator<Indexable> &toret = *this;
    operator--();
    return toret;
  }
  /// don't use this between different indexables please...
  bool operator==(DllIterator<Indexable> o) const { return i == o.i; }
  bool operator!=(DllIterator<Indexable> o) const { return i != o.i; }
};

/*
  Prim template:
  Pose getPose() const;
  AABB getAABB(v::DVec<3> velo, v::DVec<3> omega, v::DVec<3> center) const;

  Obj template:
  void setCPhysInfo(const CPhysInfo &, primIter);
  which should modify the Prims in this Obj so that getPose retrieves them,
  optionally using primIter, which is also the end iterator (the linked lists
  are circular, after all)
*/
template <class Prim, class Obj> struct World {
  struct PrimEntry {
    Prim prim;
    AABB aabb;
    int obji; /// the index of the object this primitive is in
    /// doubly linked list values for primitives in same object
    int nexti;
    int previ;
  };
  struct GetPrim {
    typedef Prim return_type;
    return_type &operator()(PrimEntry &e) const { return e.prim; }
  };
  struct ObjEntry {
    Obj obj;
    CPhysInfo cpi;
    Pose relP; /// relative pose to that of parent coalesced object
    int primi; /// index of one of the primitives this object contains
    int cobji; /// the index of the coalesced object this is in, or else 0
    /// dllist of either independent objs or all objs in same cobj
    int nexti;
    int previ;
  };
  struct CObjEntry {
    CPhysInfo cpi;
    int obji; /// the index of one of the objects this contains
    /// doubly linked list for all coalesced objects
    int nexti;
    int previ;
  };
  Pool<PrimEntry> prims;
  // element 0 of this is dummy for independent objs dllist, 1 is ground
  Pool<ObjEntry> objs;
  // element 0 of this is dummy for its dllist
  Pool<CObjEntry> cobjs;

  template <class PrimIter>
  std::enable_if<
      std::is_same<typename std::iterator_traits<PrimIter>::value_type>::value>
  addObj(const CPhysInfo &cpi, PrimIter iter, const PrimIter &iter_end) {
    int obji = objs.mkNew();
    objs[obji].cpi = cpi;
    objs[obji].cobji = 0;
    dllAddAfter(objs, 0, obji);
    int primi = prims.mkNew();
    objs[obji].primi = primi;
    prims[primi].nexti = primi;
    prims[primi].previ = primi;
    while (true) {
      prims[primi].prim = *iter;
      prims[primi].obji = obji;
      // TODO: BROAD PHASE ADD, AABB CALCULATE
      if (++iter == iter_end) {
        break;
      }
      int tprimi = prims.mkNew();
      dllAddAfter(prims, primi, tprimi);
    }
  }
  /// assumes obji is independent
  void remObj(int obji) {
    int primi = objs[obji].primi;
    int tprimi = primi;
    do {
      int tmp = prims[tprimi].nexti;
      dllRemove(prims, tprimi);
      tprimi = tmp;
    } while (tprimi != primi);
    dllRemove(objs, obji);
    objs.rem(obji);
  }

private:
  void updateObj0(int obji) {
    objs[obji].obj.setCPhysInfo(cpi, DllIterator<Pool<PrimEntry>, GetPrim>(
                                         &prims, objs[obji].primi, GetPrim()));
  }

public:
  /// assumes obji is independent
  void updateObj(const CPhysInfo &cpi, int obji) {
    objs[obji].cpi = cpi;
    updateObj0(obji);
  }
  /// assumes obji is independent
  void updateObj(const Pose &pose, int obji) {
    objs[obji].cpi.pi.pose = pose;
    updateObj0(obji);
  }

  void updateCObj(const CPhysInfo &cpi, int cobji) {
    cobjs[cobji].cpi = cpi;
    int obji = cobjs[cobji].obji;
    int tobji = obji;
    do {
      updateObj(combinePose(cpi.pi.pose, objs[tobji].relP), tobji);
      tobji = objs[tobji].nexti;
    } while (tobji != obji);
  }
};

} // namespace roller

#endif // ROLLER_PHYSICS_WORLD_HPP_
