#ifndef ROLLER_PHYSICS_WORLD_HPP_
#define ROLLER_PHYSICS_WORLD_HPP_

#include <iterator>
#include <limits>
#include <memory>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
#include <utility>

#include "broad_phase.hpp"
#include "contact.hpp"
#include "phys_info.hpp"
#include "pool.hpp"
#include "screw.hpp"

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
  p[i].previ = p[i].nexti = i;
}
template <class T> struct RetSame {
  T operator()(T t) const noexcept { return t; }
};
/*
  iterator for a (circular) doubly linked list; Fun cannot be pure reference
  Fun takes reference as input, produces reference
*/
template <class Indexable,
          class Fun = RetSame<decltype(std::declval<Indexable>()[0]) &>,
          class T = decltype(std::declval<Fun>()(std::declval<Indexable>()[0]))>
struct DllIterator {
  typedef std::ptrdiff_t difference_type;
  typedef typename std::remove_reference<T>::type value_type;
  typedef value_type *pointer;
  typedef value_type &reference;
  typedef std::bidirectional_iterator_tag iterator_category;

  typedef DllIterator<Indexable, Fun, T> ThisType;

private:
  Indexable *p;
  int i;
  Fun fun;

public:
  // say hi to DefaultConstructible
  DllIterator() : p(nullptr), i(-1) {}
  DllIterator(Indexable *p, int i) : p(p), i(i) {}
  DllIterator(Indexable *p, int i, Fun &&fun) : p(p), i(i), fun(fun) {}

  reference operator*() { return fun((*p)[i]); }
  pointer operator->() { return &fun((*p)[i]); }
  ThisType &operator++() {
    i = (*p)[i].nexti;
    return *this;
  }
  ThisType operator++(int) {
    ThisType toret = *this;
    operator++();
    return toret;
  }
  ThisType &operator--() {
    i = (*p)[i].previ;
    return *this;
  }
  ThisType operator--(int) {
    ThisType toret = *this;
    operator--();
    return toret;
  }
  bool operator==(const ThisType &o) const { return i == o.i; }
  bool operator!=(const ThisType &o) const { return i != o.i; }
};
/// wraps a LegacyInputIterator with a function returning a reference
template <
    class Iterator,
    class Fun = RetSame<typename std::iterator_traits<Iterator>::value_type &>,
    class T = decltype(std::declval<Fun>()(
        std::declval<typename std::iterator_traits<Iterator>::value_type &>()))>
struct IteratorWrapper {
  typedef std::ptrdiff_t difference_type;
  typedef typename std::remove_reference<T>::type value_type;
  typedef value_type *pointer;
  typedef value_type &reference;
  typedef std::input_iterator_tag iterator_category;

private:
  Iterator i;
  Fun fun;

  typedef IteratorWrapper<Iterator, Fun, T> ThisType;

public:
  IteratorWrapper() {}
  IteratorWrapper(const Iterator &i) : i(i) {}
  IteratorWrapper(const Iterator &i, Fun &&fun) : i(i), fun(fun) {}

  reference operator*() { return fun(*i); }
  pointer operator->() { return &fun(*i); }
  ThisType &operator++() {
    ++i;
    return *this;
  }
  ThisType &operator++(int) {
    const auto &toret = *this;
    ++i;
    return toret;
  }
  bool operator==(const ThisType &o) const { return i == o.i; }
  bool operator!=(const ThisType &o) const { return i != o.i; }
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
typename std::enable_if<
    std::is_same<typename std::iterator_traits<Iterator>::value_type,
                 CPhysInfo>::value,
    CPhysInfo>::type
combinedCPI(const Iterator &it, const Iterator &ite) {
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
    iner += adjustCenter(o->aux.rot * o->pi.iner * o->aux.rot.transpose(),
                         o->pi.mass, relp);
    tlm += o->pi.lm;
    tam += o->pi.am - cross3(relp, o->pi.lm);
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
  // return contact time on same scale as t1 and t2 plz
  ccd::Contact doCCD(double t1, double t2, ScrewM, const Prim &, ScrewM) const;

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
    Pose relP; /// cobji>0?cobjs[cobji].cpi.pi.pose.fromWorldCoords(cpi.pi.pose)
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
  // element 0 of this is dummy for its dllist (only meaningful during step)
  Pool<CObjEntry> cobjs;

  BroadPhaseAABBTree broad;

  struct Collision {
    ccd::Contact c;   /// contact time in world time
    int occobji;      /// other ccobji (cobji or -obji)
    int prim1, prim2; /// prim1 is for key, prim2 is for occobji
  };
  /// map from ccobjis to valid earliest-time collisions (nonempty only in step)
  std::unordered_map<int, Collision> collisions;
  /// objects coalesced with the ground; includes the ground (1)
  std::unordered_set<int> groundCoal;

  double g;   // acceleration of gravity in negative z
  double pad; // aabb padding

public:
  World(double g, double pad) : g(g), pad(pad) {
    objs.mkNew(); // should be 0
    objs[0].cobji = objs[0].nexti = objs[0].previ = 0;
    cobjs.mkNew(); // should be 0
    cobjs[0].obji = cobjs[0].nexti = cobjs[0].previ = 0;
  }

  /// all the public and private object addition or removal methods require that
  /// their respective ccobjis are not in collisions
  template <class PrimIter>
  typename std::enable_if<
      std::is_same<typename std::iterator_traits<PrimIter>::value_type,
                   Prim>::value,
      int>::type
  addObj(const CPhysInfo &cpi, Obj &&obj, PrimIter iter,
         const PrimIter &iter_end) {
    int obji = objs.mkNew();
    objs[obji].obj = obj;
    objs[obji].cpi = cpi;
    objs[obji].cobji = -obji;
    objs[obji].previ = objs[obji].nexti = obji;
    dllAddAfter(objs, 0, obji);
    int primi = prims.mkNew();
    objs[obji].primi = primi;
    prims[primi].nexti = primi;
    prims[primi].previ = primi;
    while (true) {
      prims[primi].prim = *iter;
      prims[primi].obji = obji;
      AABB aabb = prims[primi].prim.getAABB(getScrewM(cpi));
      prims[primi].leafi = broad.insert(primi, aabb);
      if (++iter == iter_end) {
        break;
      }
      int tprimi = prims.mkNew();
      prims[tprimi].previ = prims[tprimi].nexti = tprimi;
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
  void breakCObj(int cobji) {
    int oi = cobjs[cobji].obji;
    int toi = oi;
    do {
      int tmp = objs[toi].nexti;
      dllRemove(objs, toi);
      dllAddAfter(objs, 0, toi);
      objs[toi].cobji = -toi;
      toi = tmp;
    } while (toi != oi);
    dllRemove(cobjs, cobji);
    cobjs.rem(cobji);
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

private:
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
  typename std::enable_if<
      std::is_same<typename std::iterator_traits<Iterator>::value_type,
                   int>::value,
      int>::type
  addCObj(const Iterator &it, const Iterator &ite) {
    int cobji = cobjs.mkNew();
    cobjs[cobji].previ = cobjs[cobji].nexti = cobji;
    dllAddAfter(cobjs, 0, cobji);
    struct GetCPI {
      Pool<ObjEntry> *objs;
      CPhysInfo &operator()(int i) const { return (*objs)[i].cpi; }
    } fun{&objs};
    cobjs[cobji].cpi =
        combinedCPI(IteratorWrapper<Iterator, GetCPI &>(it, fun),
                    IteratorWrapper<Iterator, GetCPI &>(ite, fun));
    for (Iterator o = it; o != ite; ++o) {
      objs[*o].cobji = cobji;
      objs[*o].relP = objs[*o].cpi.pi.pose;
      objs[*o].relP.p -= cobjs[cobji].cpi.pi.pose.p; // a.k.a. center of mass
    }
    int obji = *it;
    dllRemove(objs, obji); // obji becomes a length-1 dll from this
    Iterator o = it;
    ++o;
    for (; o != ite; ++o) {
      // btw this ends up putting the objects in reverse order lol
      dllRemove(objs, *o);
      dllAddAfter(objs, obji, *o);
    }
    cobjs[cobji].obji = obji;
    return cobji;
  }
  /// for exactly two objects
  int addCObj(int obji1, int obji2) {
    int objis[2] = {obji1, obji2};
    return addCObj(objis + 0, objis + 2);
  }
  void addObjToCObj(int cobji, int obji) {
    CPhysInfo cpis[2] = {cobjs[cobji].cpi, objs[obji].cpi};
    CPhysInfo ncpi = combinedCPI(cpis + 0, cpis + 2);
    ncpi.pi.pose.q = cobjs[cobji].cpi.pi.pose.q;
    v::DVec<3> transl = ncpi.pi.pose.p - cobjs[cobji].cpi.pi.pose.p;
    cobjs[cobji].cpi = ncpi;
    int oi = cobjs[cobji].obji;
    int toi = oi;
    do {
      objs[toi].relP.p -= transl;
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
    CPhysInfo cpis[2] = {cobjs[cobji1].cpi, cobjs[cobji2].cpi};
    CPhysInfo ncpi = combinedCPI(cpis + 0, cpis + 2);
    ncpi.pi.pose.q = cobjs[cobji1].cpi.pi.pose.q;
    v::DVec<3> transl = ncpi.pi.pose.p - cobjs[cobji1].cpi.pi.pose.p;
    cobjs[cobji1].cpi = ncpi;
    int oi = cobjs[cobji1].obji;
    int toi = oi;
    do {
      objs[toi].relP.p -= transl;
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
    dllRemove(cobjs, cobji2);
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
    ccd::Contact con = prims[primi1].prim.doCCD(
        time, ccdIntE, getScrewM(tcpi1), prims[primi2].prim, getScrewM(tcpi2));
    if (con.t >= ccdIntE) {
      return;
    }
    int gccobji = objs[1].cobji;
    // XXX: I'm pretty sure this doesn't work as intended
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
    tcpi1.pi.pose.p += (con.t - time) * sm1.velo;
    tcpi2.pi.pose.p += (con.t - time) * sm2.velo;
    // XXX: angular velocity here is approx
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
    } else {
      coli1 = collisions.emplace(ccobji1, Collision()).first;
    }
    if (coli2 != collisions.end()) {
      int occobji = coli2->second.occobji;
      if (occobji != ccobji1) {
        collisions.erase(occobji);
      }
    } else {
      coli2 = collisions.emplace(ccobji2, Collision()).first;
    }
    coli1->second.c = coli2->second.c = con;
    coli1->second.prim1 = coli2->second.prim2 = primi1;
    coli2->second.prim1 = coli1->second.prim2 = primi2;
    coli1->second.occobji = ccobji2;
    coli2->second.occobji = ccobji1;
  }

  template <bool doPad> void updateAllAABBs(int obji, const ScrewM &sm) {
    double pd = pad;
    int primi = objs[obji].primi;
    int tprimi = primi;
    do {
      if (doPad) {
        broad.update(prims[tprimi].leafi, prims[tprimi].prim.getAABB(sm),
                     [pd](AABB aabb) -> AABB {
                       aabb.m[0] -= pd;
                       aabb.m[1] += pd;
                       return aabb;
                     });
      } else {
        broad.update(prims[tprimi].leafi, prims[tprimi].prim.getAABB(sm),
                     [](const AABB &aabb) -> AABB { return aabb; });
      }
      tprimi = prims[tprimi].nexti;
    } while (tprimi != primi);
  }

public:
  template <class Fun>
  typename std::enable_if<std::is_same<
      v::DVec<3>, decltype(std::declval<Fun>()(
                      ccd::Contact(), std::declval<Prim>(), CPhysInfo(),
                      std::declval<Prim>(), CPhysInfo()))>::value>::type
  step(double dt, Fun &&calculateImpulse) {
    for (int toi = objs[0].nexti; toi; toi = objs[toi].nexti) {
      CPhysInfo tmp = objs[toi].cpi;
      tmp.stepVelo<true>(dt, g);
      updateObj(tmp, toi);
    }
    if (groundCoal.size() > 1) {
      int gcobji = addCObj(groundCoal.begin(), groundCoal.end());
      int oi = cobjs[gcobji].obji;
      int toi = oi;
      do {
        updateAllAABBs<false>(toi, {});
        toi = objs[toi].nexti;
      } while (toi != oi);
      for (toi = objs[0].nexti; toi; toi = objs[toi].nexti) {
        updateAllAABBs<true>(toi, getScrewM(objs[toi].cpi));
      }
    } else {
      updateAllAABBs<false>(1, {});
      for (int toi = objs[0].nexti; toi; toi = objs[toi].nexti) {
        if (toi == 1) {
          continue;
        }
        updateAllAABBs<true>(toi, getScrewM(objs[toi].cpi));
      }
    }
    broad.exportInts([this, dt](int primi1, int primi2) -> void {
      // XXX: when I feel like it, change false into true
      checkCollision<false>(0, dt, primi1, primi2);
    });
    double time = 0;
    while (collisions.size()) {
      auto citer = collisions.end();
      double ctime = std::numeric_limits<double>::max();
      int ccobji = 0;
      for (auto it = collisions.begin(), ite = citer; it != ite; ++it) {
        if (it->second.c.t < ctime) {
          citer = it;
          ctime = it->second.c.t;
          ccobji = it->first;
        }
      }
      Collision col = citer->second;
      int occobji = col.occobji;
      collisions.erase(citer);
      collisions.erase(occobji);

      double adt = ctime - time;
      for (int tobji = objs[0].nexti; tobji; tobji = objs[tobji].nexti) {
        CPhysInfo tmp = objs[tobji].cpi;
        tmp.stepTime<true>(adt);
        updateObj(tmp, tobji);
      }
      for (int tcobji = cobjs[0].nexti; tcobji; tcobji = cobjs[tcobji].nexti) {
        CPhysInfo tmp = cobjs[tcobji].cpi;
        tmp.stepTime<true>(adt);
        updateCObj(tmp, tcobji);
      }
      v::DVec<3> imp =
          calculateImpulse(col.c, prims[col.prim1].prim, getCCObjCPI(ccobji),
                           prims[col.prim2].prim, getCCObjCPI(occobji));
      int obji1 = prims[col.prim1].obji;
      int obji2 = prims[col.prim2].obji;
      objs[obji1].cpi.pi.lm += imp;
      objs[obji2].cpi.pi.lm -= imp;
      objs[obji1].cpi.pi.am += cross3(col.c.p - objs[obji1].cpi.pi.pose.p, imp);
      objs[obji2].cpi.pi.am -= cross3(col.c.p - objs[obji2].cpi.pi.pose.p, imp);
      objs[obji1].cpi.updateAux();
      objs[obji2].cpi.updateAux();

      int ncobji = mergeCCObjs(ccobji, occobji);
      time = ctime;
      if (ncobji != objs[1].cobji) {
        ScrewM sm = getScrewM(cobjs[ncobji].cpi);
        int oi = cobjs[ncobji].obji;
        int toi = oi;
        do {
          updateAllAABBs<true>(toi, sm);
          toi = objs[toi].nexti;
        } while (toi != oi);
      }
      broad.exportInts(
          [this, time, dt, ncobji](int primi1, int primi2) -> void {
            if (objs[prims[primi1].obji].cobji == ncobji ||
                objs[prims[primi2].obji].cobji == ncobji) {
              checkCollision<false>(time, dt, primi1, primi2);
            }
          });
    }
    double adt = dt - time;
    if (adt > 0) {
      for (int tobji = objs[0].nexti; tobji; tobji = objs[tobji].nexti) {
        CPhysInfo tmp = objs[tobji].cpi;
        tmp.stepTime<true>(adt);
        updateObj(tmp, tobji);
      }
      for (int tcobji = cobjs[0].nexti; tcobji; tcobji = cobjs[tcobji].nexti) {
        CPhysInfo tmp = cobjs[tcobji].cpi;
        tmp.stepTime<true>(adt);
        updateCObj(tmp, tcobji);
      }
    }
    for (int tcobji = cobjs[0].nexti; tcobji;) {
      int tmp = cobjs[tcobji].nexti;
      breakCObj(tcobji);
      tcobji = tmp;
    }
  }

  template <class Fun> void applyToPrims(Fun &&fun) {
    for (int obji = objs[0].nexti; obji; obji = objs[obji].nexti) {
      int primi = objs[obji].primi;
      do {
        fun(prims[primi].prim);
        primi = prims[primi].nexti;
      } while (primi != objs[obji].primi);
    }
    for (int cobji = cobjs[0].nexti; cobji; cobji = cobjs[cobji].nexti) {
      int obji = cobjs[cobji].obji;
      do {
        int primi = objs[obji].primi;
        do {
          fun(prims[primi].prim);
          primi = prims[primi].nexti;
        } while (primi != objs[obji].primi);
        obji = objs[obji].nexti;
      } while (obji != cobjs[cobji].obji);
    }
  }
};

} // namespace roller

#endif // ROLLER_PHYSICS_WORLD_HPP_
