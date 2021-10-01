#include <iostream>

#include "hashtable_cache.hpp"

int main() {
  roller::HashtableCache<int, int> cacher(3);
  cacher[-1] = 1;
  cacher[-2] = 2;
  cacher[-3] = 3;
  cacher[-4] = 4;
  std::cout << cacher[-1] << " " << cacher[-2] << " " << cacher[-3] << " "
            << cacher[-4] << std::endl;
  cacher[-1] = 2;
  std::cout << cacher[-1] << " " << cacher[-2] << " " << cacher[-3] << " "
            << cacher[-4] << std::endl;
  cacher[1] = 0;
  std::cout << !!cacher.get(-1) << " " << !!cacher.get(-2) << " "
            << !!cacher.get(-3) << " " << !!cacher.get(-4) << std::endl;
  cacher[1] = 2;
  std::cout << !!cacher.get(-1) << " " << !!cacher.get(-2) << " "
            << !!cacher.get(-3) << " " << !!cacher.get(-4) << std::endl;
  cacher[2] = 3;
  std::cout << !!cacher.get(-1) << " " << !!cacher.get(-2) << " "
            << !!cacher.get(-3) << " " << !!cacher.get(-4) << std::endl;
  cacher[3] = 4;
  std::cout << !!cacher.get(-1) << " " << !!cacher.get(-2) << " "
            << !!cacher.get(-3) << " " << !!cacher.get(-4) << std::endl;
  cacher[4] = 5;
  std::cout << !!cacher.get(-1) << " " << !!cacher.get(-2) << " "
            << !!cacher.get(-3) << " " << !!cacher.get(-4) << std::endl;

  std::cout << cacher[1] << " " << cacher[2] << " " << cacher[3] << " "
            << cacher[4] << std::endl;
  std::cout << *cacher.get(1) << " " << *cacher.get(2) << " " << *cacher.get(3)
            << " " << *cacher.get(4) << std::endl;
}

