#ifndef ROLLER_PHYSICS_POOL_HPP_
#define ROLLER_PHYSICS_POOL_HPP_

#include <cstring>
#include <memory>

namespace roller {

template <class Entry> class Pool {
  std::unique_ptr<Entry[]> a;
  std::size_t sz;
  std::unique_ptr<int[]> st; /// stack for free entries
  std::size_t tz;            /// stack size

  /// note: only called when tz == 0
  void resize() {
    std::size_t osz = sz;
    sz *= 2;
    st.reset(new int[sz]);
    for (std::size_t i = 0; i < osz; i++) {
      st[i] = sz - 1 - i;
    }
    tz = osz;
    std::unique_ptr<Entry[]> aNew(new Entry[sz]);
    for (std::size_t i = 0; i < osz; i++) {
      aNew[i] = std::move(a[i]);
    }
    a = std::move(aNew);
  }

  void init() {
    for (std::size_t i = 0; i < sz; i++) {
      st[i] = sz - 1 - i;
    }
  }

public:
  typedef Entry value_type;

  Pool() : a(new Entry[256]), sz(256), st(new int[256]), tz(256) { init(); }
  explicit Pool(std::size_t sz)
      : a(new Entry[sz]), sz(sz), st(new int[sz]), tz(sz) {
    init();
  }

  Entry &operator[](int i) { return a[i]; }
  const Entry &operator[](int i) const { return a[i]; }
  std::size_t size() const { return sz; }

  int mkNew() {
    if (!tz) {
      resize();
    }
    return st[--tz];
  }
  void rem(int i) { st[tz++] = i; }
};

} // namespace roller

#endif // ROLLER_PHYSICS_POOL_HPP_
