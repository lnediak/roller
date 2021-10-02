#include <iostream>

#include "hashtable_cache.hpp"

struct A {
  int i = 0;
  A() {}
  A(const A &) = delete;
  A &operator=(const A &) = delete;
  const A &operator=(const A &) const = delete;
};

std::ostream &operator<<(std::ostream &os, const A &a) {
  os << a.i;
  return os;
}

int main() {
  roller::HashtableCache<int, A> cacher(4);
  cacher[-1].i = 1;
  cacher[-2].i = 2;
  cacher[-3].i = 3;
  cacher[-4].i = 4;
  std::cout << cacher[-1] << " " << cacher[-2] << " " << cacher[-3] << " "
            << cacher[-4] << std::endl;
  cacher[-1].i = 2;
  std::cout << cacher[-1] << " " << cacher[-2] << " " << cacher[-3] << " "
            << cacher[-4] << std::endl;
  cacher[1].i = 0;
  std::cout << !!cacher.get(-1) << " " << !!cacher.get(-2) << " "
            << !!cacher.get(-3) << " " << !!cacher.get(-4) << std::endl;
  cacher[1].i = 2;
  std::cout << !!cacher.get(-1) << " " << !!cacher.get(-2) << " "
            << !!cacher.get(-3) << " " << !!cacher.get(-4) << std::endl;
  cacher[2].i = 3;
  std::cout << !!cacher.get(-1) << " " << !!cacher.get(-2) << " "
            << !!cacher.get(-3) << " " << !!cacher.get(-4) << std::endl;
  cacher[3].i = 4;
  std::cout << !!cacher.get(-1) << " " << !!cacher.get(-2) << " "
            << !!cacher.get(-3) << " " << !!cacher.get(-4) << std::endl;

  /*std::cout << cacher[1] << " " << cacher[2] << " " << cacher[3] << " "
            << cacher[4] << std::endl;*/
  std::cout << *cacher.get(1) << " " << *cacher.get(2) << " " << *cacher.get(3)
            << " " << std::endl;
}

