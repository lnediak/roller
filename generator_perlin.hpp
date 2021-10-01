#ifndef GENERATOR_PERLIN_HPP_
#define GENERATOR_PERLIN_HPP_

#include <cmath>
#include <memory>
#include <random>
#include <type_traits>

#include "util.hpp"

namespace roller {

/// numGradVecs needs to be a power of 2
std::unique_ptr<double[]> getGradVecs(std::size_t numGradVecs,
                                      std::size_t numDims, unsigned seed) {
  std::unique_ptr<double[]> gradVecs(new double[numGradVecs * numDims]);
  std::mt19937 mtrand(seed);
  double pi2 = 2 * 3.141592653589792653589793238462643383;
  // assuming numGradVecs is even (which it is)
  for (std::size_t i = 0; i < numGradVecs * numDims; i += 2) {
    double u1 = (1 - mtrand() / 4294967296.0);
    double u2 = (1 - mtrand() / 4294967296.0);
    double r = std::sqrt(-2 * std::log(u1));
    gradVecs[i] = r * std::cos(pi2 * u2);
    gradVecs[i + 1] = r * std::sin(pi2 * u2);
  }
  for (std::size_t i = 0; i < numGradVecs; i++) {
    double norm = 0;
    for (std::size_t j = 0; j < numDims; j++) {
      norm += gradVecs[numDims * i + j] * gradVecs[numDims * i + j];
    }
    norm = std::sqrt(norm);
    // the likelihood of norm being less than 1e-12 in 3 dimensions
    // (chi-squared random variable) is insanely small
    if (norm < 1e-12) {
      norm = 1e-12;
    }
    for (std::size_t j = 0; j < numDims; j++) {
      gradVecs[numDims * i + j] /= norm;
    }
  }
  return gradVecs;
}

template <std::size_t N> class GeneratorPerlin {

  struct GetResult {
    typedef void thisisavvec;
    typedef double value_type;
    static const std::size_t size = 1 << N;

    double *gradVecs;
    std::size_t numGradVecsMask;
    std::size_t octave;
    const v::IVec<N> &posf;
    const double *vecs[2];

    value_type operator[](std::size_t i) const {
      const v::IVec<N> bits = BitsVec<N>{i};
      /*std::cout << "YAHALLO FOLKS!!!! GETRESULT IN THE ROOM!!!" << std::endl;
      std::cout << "i: " << i << std::endl;
      std::cout << "posf + bits: " << (posf + bits) << std::endl;*/
      const SplitVec<N> svecs{bits, {vecs[0], vecs[1]}};
      const v::DVec<N> gvec(
          gradVecs +
          N * ((v::IVecHash<N>{}(posf + bits) + octave) & numGradVecsMask));
      return v::sum(svecs * gvec);
    }
  };

public:
  struct Options {
    v::DVec<N> scale;
    std::unique_ptr<double[]> gradVecs;
    std::size_t numGradVecsMask;
    std::size_t numOctaves;
    double persistence;
  } options;

  double operator()(const v::IVec<N> &coord) const {
    double total = 0;
    double amplitude = 1;
    v::DVec<N> pos = (vint2double(coord) + 0.50001) / options.scale;
    for (std::size_t i = options.numOctaves; i--;) {
      pos += 0.5;
      v::IVec<N> posf = vfastFloorD(pos);
      v::DVec<N> vec0 = pos - vint2double(posf);
      v::DVec<N> vec1 = vec0 - 1.;
      v::DVec<N> lerp = vec0 * vec0 * (3. - 2. * vec0);
      total += vlerpAll(GetResult{options.gradVecs.get(),
                                  options.numGradVecsMask,
                                  i,
                                  posf,
                                  {vec0.data, vec1.data}},
                        lerp) *
               amplitude;
      amplitude *= options.persistence;
      pos *= 2.;
    }
    return total - 0.1;
  }
};

} // namespace hypervoxel

#endif // GENERATOR_PERLIN_HPP_

