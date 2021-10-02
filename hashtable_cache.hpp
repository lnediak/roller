#ifndef ROLLER_HASHTABLE_CACHE_HPP_
#define ROLLER_HASHTABLE_CACHE_HPP_

#include <functional>
#include <memory>
#include <unordered_map>
#include <utility>

#include "pool_linked_list.hpp"

namespace roller {

template <class T> struct NoOpFunctor {
  void operator()(T) const noexcept {}
};

template <class T, class Fun> struct Cacher {

  PoolLinkedList<T> list;
  typedef typename decltype(list)::Node lnode;
  std::size_t sz, maxSz;

  Fun fun;

  Cacher(std::size_t maxSize, Fun &&fun)
      : list(maxSize), sz(0), maxSz(maxSize), fun(fun) {}

  lnode *push(const T &t) {
    if (sz >= maxSz) {
      fun(list.back()->obj);
      list.removeFromEnd();
    } else {
      sz++;
    }
    return list.addToBeg(t);
  }
  void refresh(lnode *p) {
    T obj = p->obj;
    list.remove(p);
    list.addToBeg(obj);
  }
};

template <class U, class V, class H = std::hash<U>, class E = std::equal_to<U>>
class HashtableCache {

  struct Entry {
    V v;
    void *lp; /// pointer to node in list
  };
  typedef std::unordered_map<U, Entry, H, E> cmap;
  typedef typename cmap::iterator cmap_iter;
  typedef typename cmap::const_iterator cmap_const_iter;
  cmap map;

  struct DelCallback {
    cmap &map;
    void operator()(const cmap_const_iter &i) noexcept { map.erase(i); }
  };
  Cacher<cmap_const_iter, DelCallback> list;

public:
  HashtableCache(std::size_t maxSz) : map(3 * maxSz / 2), list(maxSz, {map}) {}

  template <class Fun> V &applyIfNew(const U &u, Fun &&fun) {
    auto res = map.emplace(std::piecewise_construct_t{}, std::make_tuple(u),
                           std::make_tuple());
    if (res.second) {
      res.first->second.lp = list.push(res.first);
      fun(res.first->second.v);
    } else {
      list.refresh((typename decltype(list)::lnode *)res.first->second.lp);
    }
    return res.first->second.v;
  }

  V &operator[](const U &u) { return applyIfNew(u, NoOpFunctor<const V &>{}); }

  V *get(const U &u) {
    auto iter = map.find(u);
    return iter == map.end() ? nullptr : &iter->second.v;
  }
};

} // namespace roller

#endif // ROLLER_HASHTABLE_CACHE_HPP_

