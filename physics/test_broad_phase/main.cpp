#include <forward_list>
#include <iostream>
#include <random>

#include "broad_phase.hpp"

// I'm lazy lol
using namespace roller;

template <class Iter> struct EzIterHash {
  std::hash<const typename std::iterator_traits<Iter>::value_type *> base;
  std::size_t operator()(Iter iter) const noexcept { return base(&*iter); }
};
struct EasyW {
  std::forward_list<AABB> d;
  typedef decltype(d.cbegin()) iterator;
  typedef EzIterHash<iterator> iter_hash;
  typedef std::equal_to<iterator> iter_eq;

  iterator begin() const { return d.begin(); }
  iterator end() const { return d.end(); }
  AABB getAABB(iterator iter) const { return *iter; }
};

/// brute force AABB broad phase
template <class Fun, class IgnFun>
void primitiveBroadPhase(EasyW::iterator iter, EasyW::iterator iter_end,
                         IgnFun &&ignFun, Fun &&callback) {
  for (auto ii = iter; ii != iter_end; ++ii) {
    auto jj = ii;
    for (++jj; jj != iter_end; ++jj) {
      if (ii->intersects(*jj) && !ignFun(ii, jj)) {
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

void printAABB(AABB a) {
  std::cout << "min: " << a.m[0] << std::endl;
  std::cout << "max: " << a.m[1] << std::endl;
}
template <class Witer> void printUPairWiter(UnorderedPair<Witer> up) {
  std::cout << "unordered pair:" << std::endl;
  printAABB(*up.a);
  printAABB(*up.b);
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

int reluI(int a) { return a < 0 ? 0 : a; }

int testBroadPhaseAABB(int numIters, int numWats, int numSubIters) {
  std::mt19937 mtrand(1); // seed
  std::uniform_real_distribution<> rdistro(-100, 100);
  std::uniform_real_distribution<> pdistro(0.1, 100);
  std::uniform_int_distribution<> ldistro(-50, 200);
  std::uniform_int_distribution<> goofy(0, 10);

  for (int spam = 0; spam < numIters; spam++) {
    if (spam % 1 == 0) {
      std::cout << "Iteration #" << spam << std::endl;
    }
    EasyW w;
    std::size_t nObjs = 0;
    BroadPhaseAABB<EasyW &> broadPhase(w);
    typedef EasyW::iterator Witer;
    std::unordered_set<UnorderedPair<Witer>,
                       UnorderedPairHash<Witer, EasyW::iter_hash>,
                       UnorderedPairEq<Witer, EasyW::iter_eq>>
        ignoreList;
    auto ignFun = [&ignoreList](Witer a, Witer b) -> bool {
      return ignoreList.count({a, b});
    };
    for (int wat = 0; wat < numWats; wat++) {
      for (auto iter = w.d.before_begin(), iter_end = w.d.end();
           iter != iter_end;) {
        if (goofy(mtrand)) {
          ++iter;
        } else {
          auto tmp = iter++;
          if (iter == iter_end) {
            break;
          }
          for (auto it = ignoreList.begin(), it_e = ignoreList.end();
               it != it_e; ++it) {
            if (it->a == iter || it->b == iter) {
              ignoreList.erase(it);
              break;
            }
          }
          iter = w.d.erase_after(tmp);
          nObjs--;
        }
      }
      for (int ii = 0, sz = reluI(ldistro(mtrand)); ii < sz; ii++) {
        w.d.push_front(randAABB(mtrand, rdistro, pdistro));
        if (nObjs && !goofy(mtrand)) {
          std::uniform_int_distribution<> veryGoofy(1, nObjs);
          std::size_t otherI = veryGoofy(mtrand);
          Witer otherIt = w.begin();
          for (std::size_t jj = 0; jj < otherI; jj++) {
            ++otherIt;
          }
          ignoreList.insert({w.begin(), otherIt});
        }
        nObjs++;
      }
      // TODO: update when callback works with changedW
      broadPhase.updateAABBStuff(ignFun);
      decltype(ignoreList) reported, adjusted, bruteForced;
      broadPhase.exportInts([&adjusted](Witer a, Witer b) -> void {
        adjusted.insert({a, b});
      });
      primitiveBroadPhase(w.begin(), w.end(), ignFun,
                          [&bruteForced](Witer a, Witer b) -> void {
                            bruteForced.insert({a, b});
                          });
      if (!unorderedSetEq(adjusted, bruteForced)) {
        std::cout << "initial report did not match bruteForced" << std::endl;
        std::cout << "This was iteration #" << spam << ", wat #" << wat
                  << std::endl;
        return 1;
      }
      for (int spam2 = 0; spam2 < numSubIters; spam2++) {
        for (auto iter = w.d.begin(), iter_end = w.d.end(); iter != iter_end;
             ++iter) {
          if (!goofy(mtrand)) {
            // normally the changes are going to be much more coherent than this
            *iter = randAABB(mtrand, rdistro, pdistro);
          }
        }
        bool fail = false;
        broadPhase.updateAABBStuff(
            ignFun, false,
            [&adjusted, &fail](Witer a, Witer b, bool isAdd) -> void {
              if (isAdd) {
                if (!adjusted.insert({a, b}).second) {
                  std::cout << "adjusted already contained {a, b}" << std::endl;
                  fail = true;
                }
              } else {
                if (!adjusted.erase({a, b})) {
                  std::cout << "adjusted did not contain {a, b}" << std::endl;
                  fail = true;
                }
              }
            });
        reported.clear();
        broadPhase.exportInts([&reported](Witer a, Witer b) -> void {
          reported.insert({a, b});
        });
        bruteForced.clear();
        primitiveBroadPhase(w.begin(), w.end(), ignFun,
                            [&bruteForced](Witer a, Witer b) -> void {
                              bruteForced.insert({a, b});
                            });
        if (!unorderedSetEq(reported, bruteForced)) {
          std::cout << "reported did not match bruteForced" << std::endl;
          fail = true;
        }
        if (!unorderedSetEq(reported, adjusted)) {
          std::cout << "reported did not match adjusted" << std::endl;
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

int main() { testBroadPhaseAABB(10, 10, 10); }
