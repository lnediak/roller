#ifndef ROLLER_HASHTABLE_CACHE_HPP_
#define ROLLER_HASHTABLE_CACHE_HPP_

#include <functional>
#include <memory>
#include <unordered_map>
#include <utility>

namespace roller {

template <class T> struct NoOpFunctor {
  void operator()(T *) const noexcept {}
};

/**
  A weird queue (that doesn't have a pop operation) that deletes an element for
  every push after a certain size, and calls a callback each time that is done.
  `maxSzMask` needs to be a power of two minus 1.
*/
template <class T, class F = NoOpFunctor<T>> class CacheQueue {

  std::unique_ptr<T[]> data;
  /**
    lower `index` means newer, and the idem at `index` is stored at
    `data[(i + index) & maxSzMask]`. sz is current size.
  */
  std::size_t i, sz, maxSzMask;

  F fun;

public:
  CacheQueue(std::size_t maxSzMask, const F &fun)
      : data(new T[maxSzMask + 1]), i(0), sz(0), maxSzMask(maxSzMask),
        fun(fun) {}
  explicit CacheQueue(std::size_t maxSzMask)
      : data(new T[maxSzMask + 1]), i(0), sz(0), maxSzMask(maxSzMask), fun() {}
  CacheQueue() : data(new T[8]), i(0), sz(0), maxSzMask(7), fun() {}

  void clear() { i = sz = 0; }

  T *push(const T &t) {
    i = (i - 1) & maxSzMask;
    if (sz > maxSzMask) {
      fun(data[i]);
    } else {
      sz++;
    }
    data[i] = t;
    return &data[i];
  }
};

template <class U, class V, class H = std::hash<U>, class E = std::equal_to<U>>
class HashtableCache {

  typedef std::unordered_map<U, V, H, E> cmap;
  typedef typename cmap::iterator cmap_iter;
  typedef typename cmap::const_iterator cmap_const_iter;
  cmap map;

  struct DelCallback {
    cmap &map;
    void operator()(const cmap_const_iter &i) noexcept { map.erase(i); }
  };
  CacheQueue<cmap_const_iter, DelCallback> list;

public:
  HashtableCache(std::size_t maxSzMask)
      : map(3 * maxSzMask / 2), list(maxSzMask, {map}) {}

  V &operator[](const U &u) {
    auto res = map.emplace(std::piecewise_construct_t{}, std::make_tuple(u),
                           std::make_tuple());
    if (res.second) {
      list.push(res.first);
    }
    return res.first->second;
  }

  template <class Fun> V &applyIfNew(const U &u, Fun &&fun) {
    auto res = map.emplace(std::piecewise_construct_t{}, std::make_tuple(u),
                           std::make_tuple());
    if (res.second) {
      list.push(res.first);
      fun(res.first->second);
    }
    return res.first->second;
  }

  V *get(const U &u) {
    auto iter = map.find(u);
    return iter == map.end() ? nullptr : &iter->second;
  }
};

} // namespace roller

#endif // ROLLER_HASHTABLE_CACHE_HPP_

