#ifndef ROLLER_PHYSICS_WORLD_HPP_
#define ROLLER_PHYSICS_WORLD_HPP_

#include <iterator>
#include <limits>
#include <memory>
#include <type_traits>
#include <types>
#include <unordered_map>
#include <unordered_set>
#include <utility>

#include "broad_phase.hpp"
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
/// inserts the dll specified by ni after i
template <class Indexable> void dllAddAllAfter(Indexable &p, int i, int ni) {
  int nie = p[ni].previ; // end of dll specified by ni
  p[ni].previ = i;
  p[nie].nexti = p[i].nexti;
  p[p[i].nexti].previ = nie;
  p[i].nexti = ni;
}
template <class Indexable> void dllRemove(Indexable &p, int i) {
  p[p[i].previ].nexti = p[i].nexti;
  p[p[i].nexti].previ = p[i].previ;
}
template <class T> struct RetSame {
  T operator()(const T &t) const noexcept { return t; }
};
/// iterator for a (circular) doubly linked list; Fun cannot be pure reference
template <class Indexable, class Fun = RetSame,
          class T = decltype(std::declval<Fun>()(std::declval<Indexable>()[0]))>
struct DllIterator {
  typedef std::ptrdiff_t difference_type;
  typedef T value_type;
  typedef value_type *pointer;
  typedef value_type &reference;
  typedef std::bidirectional_iterator_tag iterator_category;

  Indexable *p;
  int i;
  Fun fun;

  // say hi to DefaultConstructible
  DllIterator() : p(nullptr), i(-1) {}
  DllIterator(Indexable *p, int i) : p(p), i(i) {}
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
  bool operator==(const DllIterator<Indexable> &o) const { return i == o.i; }
  bool operator!=(const DllIterator<Indexable> &o) const { return i != o.i; }
};
/// wraps a LegacyInputIterator
template <
    class Iterator, class Fun = RetSame,
    class T = decltype(std::declval<Fun>()(
        std::declval<typename std::iterator_traits<Iterator>::value_type>()))>
class IteratorWrapper {
  typedef std::ptrdiff_t difference_type;
  typedef T value_type;
  typedef value_type *pointer;
  typedef value_type &reference;
  typedef std::input_iterator_tag iterator_category;

  Iterator i;
  Fun fun;

  IteratorWrapper() {}
  IteratorWrapper(const Iterator &i) : i(i) {}
  IteratorWrapper(const Iterator &i, Fun &&fun) : i(i), fun(fun) {}

  reference operator*() { return *i; }
  pointer operator->() { return &*i; }
  decltype(*this) &operator++() {
    ++i;
    return *this;
  }
  decltype(*this) &operator++(int) {
    const auto &toret = *this;
    ++i;
    return toret;
  }
  bool operator==(const decltype(*this) &o) const { return i == o.i; }
  bool operator!=(const decltype(*this) &o) const { return i != o.i; }
};

/// df is the new center relative to the center of mass
DMat3x3 adjustCenter(const DMat3x3 &iner, double mass, const v::DVec<3> &df) {
  v::DVec<3> mm = mass * df;
  double xx = mm[0] * df[0];
  double xy = mm[0] * df[1];
  double xz = mm[0] * df[2];
  double yy = mm[1] * df[1];
  double yz = mm[1] * df[2];
  double zz = mm[2] * df[2];
  DMat3x3 dblCrossMat = {
      {yy + zz, -xy, -xz}, {-xy, xx + zz, -yz}, {-xz, -yz, xx + yy}};
  return iner + dblCrossMat;
}

template <class Iterator>
std::enable_if<std::is_same<typename std::iterator_traits<Iterator>::value_type,
                            CPhysInfo>::value,
               CPhysInfo>
combineCPI(const Iterator &it, const Iterator &ite) {
  PhysInfo pi;
  double tm = 0;             // total mass
  v::DVec<3> cm = {0, 0, 0}; // center of mass
  for (Iterator o = it; o != ite; ++o) {
    double m = o->pi.mass;
    v::DVec<3> c = o->pi.pose.p;
    if (!m) {
      tm = 0;
      cm = c;
      break;
    }
    tm += m;
    cm += m * c;
  }
  if (!tm) {
    pi.pose.p = cm;
    return pi; // the rest are all 0 by default so we're fine
  }
  pi.mass = tm;
  double tmi = 1 / tm;
  cm *= tmi;
  pi.massi = tmi;
  pi.pose.p = cm;
  v::DVec<3> tlm = {0, 0, 0};
  DMat3x3 iner = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
  v::DVec<3> tam = {0, 0, 0};
  for (Iterator o = it; o != ite; ++o) {
    v::DVec<3> relp = cm - o->pi.pose.p;
    if (!fat) {
      iner += adjustCenter(o->aux.riner, o->pi.mass, relp);
      tlm += o->pi.lm;
      tam += o->pi.am - cross3(relp, o->pi.lm);
    }
  }
  pi.lm = tlm;
  pi.am = tam;
  pi.iner = iner;
  pi.ineri = iner.inverse();
  return pi;
}

/*
  Prim template:
  Pose getPose() const;
  AABB getAABB(const ScrewM &) const;
  // returned contact time in world time plz
  Contact doCCD(double t1, double t2, ScrewM, const Prim &, ScrewM) const;

  Obj template:
  void setPose(const Pose &, primIter);
  which should modify the Prims in this Obj so that getPose retrieves them,
  optionally using primIter, which is also the end iterator (the linked lists
  are circular, after all)
*/
template <class Prim, class Obj> class World {
  struct PrimEntry {
    Prim prim;
    int leafi; /// index of leaf in broad phase
    int obji;  /// the index of the object this primitive is in
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
    int cobji; /// the index of the cobj this is in, or else -(this's index)
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

  BroadPhaseAABBTree broad;

  struct Collision {
    Contact c;        /// contact time in world time
    int occobji;      /// other ccobji
    int obji1, obji2; /// obji1 is for key, obji2 is for occobji
  };
  /// map from ccobjis to valid earliest-time collisions
  std::unordered_map<int, Collision> collisions;
  /// objects coalesced with the ground; includes the ground (1)
  std::unordered_set<int> groundCoal;

  double g;   // acceleration of gravity in negative z
  double pad; // aabb padding

public:
  // TODO: CONSTRUCTOR

  template <class PrimIter>
  std::enable_if<
      std::is_same<typename std::iterator_traits<PrimIter>::value_type,
                   Prim>::value,
      int>
  addObj(const CPhysInfo &cpi, PrimIter iter, const PrimIter &iter_end) {
    int obji = objs.mkNew();
    objs[obji].cpi = cpi;
    objs[obji].cobji = -obji;
    dllAddAfter(objs, 0, obji);
    int primi = prims.mkNew();
    objs[obji].primi = primi;
    prims[primi].nexti = primi;
    prims[primi].previ = primi;
    while (true) {
      prims[primi].prim = *iter;
      prims[primi].obji = obji;
      AABB aabb = prims[primi].prim.getAABB({});
      aabb.m[0] -= pad;
      aabb.m[1] += pad;
      prims[primi].leafi = broad.insert(primi, aabb);
      if (++iter == iter_end) {
        break;
      }
      int tprimi = prims.mkNew();
      dllAddAfter(prims, primi, tprimi);
      primi = tprimi;
    }
    return obji;
  }
  /// assumes obji is independent
  void remObj(int obji) {
    int primi = objs[obji].primi;
    int tprimi = primi;
    do {
      int tmp = prims[tprimi].nexti;
      broad.remove(prims[tprimi].leafi);
      dllRemove(prims, tprimi);
      prims.rem(tprimi);
      tprimi = tmp;
    } while (tprimi != primi);
    dllRemove(objs, obji);
    objs.rem(obji);
  }

private:
  void updateObj0(const Pose &pose, int obji) {
    objs[obji].obj.setPose(
        pose, DllIterator<Pool<PrimEntry>, GetPrim>(&prims, objs[obji].primi));
  }
  void updateCObj0(const Pose &pose, int obji) {
    int tobji = obji;
    do {
      updateObj(pose.toWorldCoords(objs[tobji].relP), tobji);
      tobji = objs[tobji].nexti;
    } while (tobji != obji);
  }

public:
  /// assumes obji is independent
  void updateObj(const CPhysInfo &cpi, int obji) {
    objs[obji].cpi = cpi;
    updateObj0(cpi.pi.pose, obji);
  }
  /// assumes obji is independent
  void updateObj(const Pose &pose, int obji) {
    objs[obji].cpi.pi.pose = pose;
    updateObj0(pose, obji);
  }

  void updateCObj(const CPhysInfo &cpi, int cobji) {
    cobjs[cobji].cpi = cpi;
    updateCObj0(cpi.pi.pose, cobjs[cobji].obji);
  }
  void updateCObj(const Pose &pose, int cobji) {
    cobjs[cobji].cpi.pi.pose = pose;
    updateCObj0(pose, cobjs[cobji].obji);
  }

  /// takes iterators for object indices
  template <class Iterator>
  std::enable_if<
      std::is_same<typename std::iterator_traits<Iterator>::value_type,
                   int>::value,
      int>
  addCObj(const Iterator &it, const Iterator &ite) {
    int cobji = cobjs.mkNew();
    dllAddAfter(cobjs, 0, cobji);
    struct GetCPI {
      Pool<ObjEntry> *objs;
      CPhysInfo operator()(int i) const { return (*objs)[i].cpi; }
    } fun{&objs};
    cobjs[cobji].cpi =
        combinedCPI(IteratorWrapper<Iterator, GetCPI>(it, fun), ite);
    for (Iterator o = it; o != ite; ++o) {
      objs[*o].cobji = cobji;
      objs[*o].relP = objs[*o].cpi.pi.pose;
      objs[*o].relP.p -= cobjs[cobji].cpi.pi.pose.p; // a.k.a. center of mass
    }
    int obji = *it;
    objs[obji].nexti = obji;
    objs[obji].previ = obji;
    Iterator o = it;
    ++o;
    for (; o != ite; ++o) {
      // btw this ends up putting the objects in reverse order lol
      dllRemove(*o);
      dllAddAfter(objs, obji, *o);
    }
    cobjs[cobji].obji = obji;
    return cobji;
  }
  /// for exactly two objects
  int addCObj(int obji1, int obji2) {
    int objis[2] = {obji1, obji2};
    return addCObj(objis, objis + 2);
  }
  void addObjToCObj(int cobji, int obji) {
    CPhysInfo ncpi = combinedCPI(cobjs[cobji].cpi, objs[obji].cpi);
    ncpi.pi.pose.q = cobjs[cobji].cpi.pi.pose.q;
    v::DVec<3> transl = ncpi.pi.pose.p - cobjs[cobji].cpi.pi.pose.p;
    cobjs[cobji].cpi = ncpi;
    int oi = cobjs[cobji].obji;
    int toi = oi;
    do {
      objs[toi].relP -= transl;
      toi = objs[toi].nexti;
    } while (toi != oi);
    dllRemove(objs, obji);
    dllAddAfter(objs, cobjs[cobji].obji, obji);
    objs[obji].cobji = cobji;
    objs[obji].relP = ncpi.pi.pose.fromWorldCoords(objs[obji].cpi.pi.pose);
  }
  int mergeCObjs(int cobji1, int cobji2) {
    struct GetRelP {
      typedef Pose return_type;
      return_type &operator()(ObjEntry &e) const { return e.relP; }
    };
    CPhysInfo ncpi = combinedCPI(cobjs[cobji1].cpi, cobjs[cobji2].cpi);
    ncpi.pi.pose.q = cobjs[cobji1].cpi.pi.pose.q;
    v::DVec<3> transl = ncpi.pi.pose.p - cobjs[cobji1].cpi.pi.pose.p;
    cobjs[cobji1].cpi = ncpi;
    int oi = cobjs[cobji1].obji;
    int toi = oi;
    do {
      objs[toi].relP -= transl;
      toi = objs[toi].nexti;
    } while (toi != oi);
    Pose shiftP = ncpi.pi.pose.fromWorldCoords(cobjs[cobji2].cpi.pi.pose);
    oi = cobjs[cobji2].obji;
    toi = oi;
    do {
      objs[toi].cobji = cobji1;
      // relying on associativity of Pose::toWorldCoords
      objs[toi].relP = shiftP.toWorldCoords(objs[toi].relP);
      toi = objs[toi].nexti;
    } while (toi != oi);
    dllAddAllAfter(objs, cobjs[cobji1].obji, cobjs[cobji2].obji);
    dllRemove(cobji2);
    cobjs.rem(cobji2);
    return cobji1;
  }
  int mergeCCObjs(int ccobji1, int ccobji2) {
    if (ccobji1 < 0) {
      if (ccobji2 < 0) {
        return addCObj(-ccobji1, -ccobji2);
      }
      addObjToCObj(ccobji2, -ccobji1);
      return ccobji2;
    }
    if (ccobji2 < 0) {
      addObjToCObj(ccobji1, -ccobji2);
      return ccobji1;
    }
    return mergeCObjs(ccobji1, ccobji2);
  }

private:
  CPhysInfo getCCObjCPI(int ccobji) const {
    return ccobji < 0 ? objs[-ccobji].cpi : cobjs[ccobji].cpi;
  }

  static bool isSmol(const ScrewM &sm) {
    return l1norm(sm.velo) < 1e-2 && l1norm(sm.omega) < 1e-3;
  }
  static ScrewM getScrewM(const CPhysInfo &cpi) {
    return {cpi.aux.velo, cpi.aux.omega, cpi.pi.pose.p};
  }
  ScrewM getScrewM(int ccobji) const { return getScrewM(getCCObjCPI(ccobji)); }

  double getColTime(int ccobji) const {
    auto res = collisions.find(ccobji);
    if (res == collisions.end()) {
      return std::numeric_limits<double>::max();
    }
    return res->second.c.t;
  }
  static double dmin(double a, double b) { return a < b ? a : b; }
  /// a helper for step
  template <bool doC>
  void checkCollision(double time, double dt, int primi1, int primi2) {
    int obji1 = prims[primi1].obji;
    int obji2 = prims[primi2].obji;
    int ccobji1 = objs[obji1].cobji;
    int ccobji2 = objs[obji2].cobji;
    if (ccobji1 == ccobji2) {
      return;
    }
    auto coli1 = collisions.find(ccobji1);
    auto coli2 = collisions.find(ccobji2);
    double colt1 = std::numeric_limits<double>::max();
    double colt2 = colt1;
    if (coli1 != collisions.end()) {
      colt1 = coli1->second.c.t;
    }
    if (coli2 != collisions.end()) {
      colt2 = coli2->second.c.t;
    }
    double ccdIntE = dmin(dt, dmin(colt1, colt2));
    CPhysInfo tcpi1 = getCCObjCPI(ccobji1);
    ScrewM sm1 = getScrewM(tcpi1);
    CPhysInfo tcpi2 = getCCObjCPI(ccobji2);
    ScrewM sm2 = getScrewM(tcpi2);
    Contact con = prims[primi1].prim.doCCD(
        time, ccdIntE, getScrewM(tcpi1), prims[primi2].prim, getScrewM(tcpi2));
    if (con.t >= ccdIntE) {
      return;
    }
    int gccobji = getCCObj(1);
    if (doC && ((ccobji1 == gccobji && isSmol(sm2)) ||
                (ccobji2 == gccobji && isSmol(sm1)))) {
      int ngccobji = ccobji1 == gccobji ? ccobji2 : ccobji1;
      if (ngccobji < 0) {
        groundCoal.insert(-ngccobji);
      } else {
        int oi = cobjs[ngccobji].obji;
        int toi = oi;
        do {
          groundCoal.insert(toi);
          toi = objs[toi].nexti;
        } while (toi != oi);
      }
      int ncobji = mergeCCObjs(ccobji1, ccobji2);
      if (coli1 != collisions.end()) {
        int occobji = coli1->second.occobji;
        if (occobji == ccobji2) {
          // a rather interesting situation...
          collisions.erase(coli1);
          collisions.erase(coli2);
          return;
        }
        if (ncobji != ccobji1) {
          // like if ccobji1 was not a cobj
          collisions[ncobji] = coli1->second;
          collisions.erase(coli1);
          collisions[occobji].occobji = ncobji;
        }
      }
      if (coli2 != collisions.end()) {
        int occobji = coli2->second.occobji;
        collisions[ncobji] = coli2->second;
        collisions.erase(coli2);
        collisions[occobji].occobji = ncobji;
      }
      return;
    }
    // note that these do not affect the aux, so updateAux is not necessary
    tcpi1.pi.pose.p += ccdIntE * sm1.velo;
    tcpi2.pi.pose.p += ccdIntE * sm2.velo;
    v::DVec<3> urel = tcpi1.getVelocity(con.p) - tcpi2.getVelocity(con.p);
    if (v::dot(urel, con.n) > 0) {
      return;
    }
    if (ccobji1 == gccobji || ccobji2 == gccobji) {
      if (prims[primi1].obji != 1 && prims[primi2].obji != 1) {
        groundCoal.clear();
      }
    }
    if (coli1 != collisions.end()) {
      int occobji = coli1->second.occobji;
      if (occobji != ccobji2) {
        collisions.erase(occobji);
      }
    }
    if (coli2 != collisions.end()) {
      int occobji = coli2->second.occobji;
      if (occobji != ccobji1) {
        collisions.erase(occobji);
      }
    }
    coli1->second.c = coli2->second.c = con;
    coli1->second.obji1 = coli2->second.obji2 = obji1;
    coli2->second.obji1 = coli1->second.obji2 = obji2;
    coli1->second.occobji = ccobji2;
    coli2->second.occobji = ccobji1;
  }
};

} // namespace roller

#endif // ROLLER_PHYSICS_WORLD_HPP_
