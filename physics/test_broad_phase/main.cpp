#include <functional>
#include <iostream>
#include <random>

#include "broad_phase.hpp"

// I'm lazy lol
using namespace roller;

struct EasyW {
  struct Entry {
    AABB aabb;
    int leafi;
    int next;
  };
  /// 0 is head, aabb meaningless
  Pool<Entry> d;

  EasyW() : d(256) {
    d.mkNew();
    d[0].next = 0;
  }

  struct LazyIter {
    typedef Entry value_type;
    Pool<Entry> *d;
    int i;
    value_type &operator*() { return (*d)[i]; }
    const value_type &operator*() const { return (*d)[i]; }
    value_type *operator->() { return &(*d)[i]; }
    const value_type *operator->() const { return &(*d)[i]; }
    bool operator==(const LazyIter &o) const noexcept { return i == o.i; }
    bool operator!=(const LazyIter &o) const noexcept { return !operator==(o); }
    LazyIter &operator++() {
      i = (*d)[i].next;
      return *this;
    }
    LazyIter operator++(int) {
      LazyIter toret = *this;
      operator++();
      return toret;
    }
  };
  struct LazyIterHash {
    std::size_t operator()(const LazyIter &o) const noexcept { return o.i; }
  };
  typedef LazyIter iterator;
  typedef LazyIterHash iter_hash;
  typedef std::equal_to<iterator> iter_eq;

  iterator atInd(int i) { return {&d, i}; }
  iterator begin() { return atInd(d[0].next); }
  iterator end() { return atInd(0); }
  iterator before_begin() { return end(); }
  iterator insert_after(iterator iter, const AABB &value, int leafi) {
    int newn = d.mkNew();
    d[newn].aabb = value;
    d[newn].leafi = leafi;
    d[newn].next = iter->next;
    iter->next = newn;
    return {&d, newn};
  }
  iterator push_front(const AABB &value, int leafi) {
    return insert_after(before_begin(), value, leafi);
  }
  iterator erase_after(iterator iter) {
    iter->next = d[iter->next].next;
    return ++iter;
  }
};

/// brute force AABB broad phase
template <class Fun>
void primitiveBroadPhase(EasyW::iterator iter, EasyW::iterator iter_end,
                         Fun &&callback) {
  for (auto ii = iter; ii != iter_end; ++ii) {
    auto jj = ii;
    for (++jj; jj != iter_end; ++jj) {
      if (ii->aabb.intersects(jj->aabb)) {
        callback(ii, jj);
      }
    }
  }
}

AABB randAABB(std::mt19937 &mtrand, std::uniform_real_distribution<> &rdistro,
              std::uniform_real_distribution<> &pdistro) {
  AABB aabb;
  for (int i = 0; i < 3; i++) {
    aabb.m[0][i] = rdistro(mtrand);
    aabb.m[1][i] = aabb.m[0][i] + pdistro(mtrand);
  }
  return aabb;
}

template <class Witer> void printUPairWiter(UnorderedPair<Witer> up) {
  std::cout << "unordered pair: " << up.a.i << " " << up.b.i << std::endl;
  up.a->aabb.print();
  up.b->aabb.print();
}

/// massively slow equality comparison between two unordered_sets that doesn't
/// rely on Key's operator==
template <class Key, class KeyHash, class KeyEq>
bool unorderedSetEq(const std::unordered_set<Key, KeyHash, KeyEq> &a,
                    const std::unordered_set<Key, KeyHash, KeyEq> &b) {
  // this order makes debugging easier (vs just doing the size comparison first)
  for (auto &e : a) {
    if (a.count(e) != b.count(e)) {
      std::cout << "found mismatch! in a" << std::endl;
      printUPairWiter(e);
      return false;
    }
  }
  for (auto &e : b) {
    if (a.count(e) != b.count(e)) {
      std::cout << "found mismatch! in b" << std::endl;
      printUPairWiter(e);
      return false;
    }
  }
  if (a.size() != b.size()) {
    std::cout << "size mismatch: " << a.size() << " " << b.size() << std::endl;
    return false;
  }
  return true;
}
/// determines if each element in a is also in b
template <class Key, class KeyHash, class KeyEq>
bool uSetSubsetOf(const std::unordered_set<Key, KeyHash, KeyEq> &a,
                  const std::unordered_set<Key, KeyHash, KeyEq> &b) {
  for (auto &e : a) {
    if (a.count(e) != b.count(e)) {
      std::cout << "found element in a, not in b" << std::endl;
      printUPairWiter(e);
      return false;
    }
  }
  return true;
}

int reluI(int a) { return a < 0 ? 0 : a; }

int testBroadPhaseAABB(int numIters, int numWats, int numSubIters) {
  std::mt19937 mtrand(1); // seed
  std::uniform_real_distribution<> rdistro(-100, 100);
  std::uniform_real_distribution<> pdistro(0.1, 100);
  std::uniform_int_distribution<> ldistro(-50, 200);
  // std::uniform_int_distribution<> ldistro(-5, 20);
  std::uniform_int_distribution<> goofy(0, 10);

  for (int spam = 0; spam < numIters; spam++) {
    if (spam % 1 == 0) {
      std::cout << "Iteration #" << spam << std::endl;
    }
    EasyW w;
    std::size_t nObjs = 0;
    BroadPhaseAABBTree broadPhase;
    typedef EasyW::iterator Witer;
    for (int wat = 0; wat < numWats; wat++) {
      for (auto iter = w.before_begin(), iter_end = w.end();
           iter != iter_end;) {
        if (goofy(mtrand)) {
          ++iter;
        } else {
          auto tmp = iter++;
          if (iter == iter_end) {
            break;
          }
          broadPhase.remove(iter->leafi);
          iter = w.erase_after(tmp);
          nObjs--;
        }
      }
      for (int ii = 0, sz = reluI(ldistro(mtrand)); ii < sz; ii++) {
        AABB aabb = randAABB(mtrand, rdistro, pdistro);
        auto iter = w.push_front(aabb, 0);
        iter->leafi = broadPhase.insert(iter.i, aabb);
        nObjs++;
      }
      if (!broadPhase.testStructure()) {
        std::cout << "broadPhase tree invalid" << std::endl;
        return 1;
      }
      std::unordered_set<UnorderedPair<Witer>,
                         UnorderedPairHash<Witer, EasyW::iter_hash>,
                         UnorderedPairEq<Witer, EasyW::iter_eq>>
          reported, bruteForced;
      broadPhase.exportInts([&w, &reported](int a, int b) -> void {
        reported.insert({w.atInd(a), w.atInd(b)});
      });
      primitiveBroadPhase(w.begin(), w.end(),
                          [&bruteForced](Witer a, Witer b) -> void {
                            bruteForced.insert({a, b});
                          });
      if (!wat && !unorderedSetEq(reported, bruteForced)) {
        std::cout << "initial report did not match bruteForced" << std::endl;
        std::cout << "This was iteration #" << spam << ", wat #" << wat
                  << std::endl;
        return 1;
      }
      if (!uSetSubsetOf(bruteForced, reported)) {
        std::cout << "bruteForced not subset of reported" << std::endl;
        std::cout << "This was iteration #" << spam << ", wat #" << wat
                  << std::endl;
        return 1;
      }
      for (int spam2 = 0; spam2 < numSubIters; spam2++) {
        if (!broadPhase.testStructure()) {
          std::cout << "broadPhase tree invalid" << std::endl;
        }
        for (auto iter = w.begin(), iter_end = w.end(); iter != iter_end;
             ++iter) {
          if (!goofy(mtrand)) {
            iter->aabb = randAABB(mtrand, rdistro, pdistro);
            iter->leafi = broadPhase.update(
                iter->leafi, iter->aabb,
                [](const AABB &aabb) -> AABB { return aabb; });
          }
        }
        bool fail = false;
        if (!broadPhase.testStructure()) {
          std::cout << "broadPhase tree invalid" << std::endl;
          fail = true;
        }
        reported.clear();
        broadPhase.exportInts([&w, &reported](int a, int b) -> void {
          reported.insert({w.atInd(a), w.atInd(b)});
        });
        bruteForced.clear();
        primitiveBroadPhase(w.begin(), w.end(),
                            [&bruteForced](Witer a, Witer b) -> void {
                              bruteForced.insert({a, b});
                            });
        if (!uSetSubsetOf(bruteForced, reported)) {
          std::cout << "bruteForced not subset of reported" << std::endl;
          fail = true;
        }
        if (fail) {
          std::cout << "This was iteration #" << spam << ", wat #" << wat
                    << ", sub-iteration #" << spam2 << std::endl;
          return 1;
        }
      }
    }
  }
  return 0;
}

int main() { testBroadPhaseAABB(30, 10, 10); }
