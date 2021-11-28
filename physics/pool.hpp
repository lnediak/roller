#ifndef ROLLER_PHYSICS_POOL_HPP_
#define ROLLER_PHYSICS_POOL_HPP_

#include <memory>

namespace roller {

template <class Entry> class Pool {
  std::unique_ptr<Entry> a;
  std::size_t sz;
  std::unique_ptr<int> st; /// stack for free entries

  void resize() {
    std::size_t osz = sz;
    sz *= 2;
    st.reset(new int[sz]);
    for (std::size_t i = osz; i < sz; i++) {
      st[i] = sz - 1 - i;
    }
    std::unique_ptr<Entry> aNew(new Entry[sz]);
    std::memcpy(aNew.get(), a.get(), osz * sizeof(Entry));
    a = std::move(aNew);
  }

  void init() {
    for (std::size_t i = 0; i < sz; i++) {
      st[i] = sz - 1 - i;
    }
  }

public:
  typedef Entry value_type;

  Pool() : a(new Entry[256]), sz(256), st(new int[256]) { init(); }
  explicit Pool(std::size_t sz) : a(new Entry[sz]), sz(sz), st(new int[sz]) {
    init();
  }

  Entry &operator[](int i) { return a[i]; }
  const Entry &operator[](int i) { return a[i]; }
  std::size_t size() const { return sz; }

  int mkNew() {
    if (st.empty()) {
      resize();
    }
    int toret = st.back();
    st.pop_back();
  }
  void rem(int i) { st.push_back(i); }
};

} // namespace roller

#endif // ROLLER_PHYSICS_POOL_HPP_
