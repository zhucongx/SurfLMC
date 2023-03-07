#ifndef SLMC_SLMC_CFG_INCLUDE_CONSTANTS_HPP_
#define SLMC_SLMC_CFG_INCLUDE_CONSTANTS_HPP_

#include <cmath>
#include <array>

constexpr double kEpsilon = 1e-8;



namespace constants {
constexpr double kLatticeConstant = 4.046;
constexpr double kFirstNearestNeighborsCutoff = 3.5;
constexpr double kSecondNearestNeighborsCutoff = 4.8;
constexpr double kThirdNearestNeighborsCutoff = 5.3;
// length between fourth nearest neighbors is double of length of first nearest neighbors
constexpr double kNearNeighborsCutoff = kThirdNearestNeighborsCutoff;

constexpr size_t kNumThirdNearestSetSizeOfPair = 60;
constexpr size_t kNumThirdNearestSetSizeOfSite = 1 + 12 + 6 + 24;

constexpr size_t kNumFirstNearestNeighbors = 12;
constexpr size_t kNumSecondNearestNeighbors = 6;
constexpr size_t kNumThirdNearestNeighbors = 24;





constexpr double kBoltzmann = 8.617333262145e-5;
constexpr double kPrefactor = 1e13;
} // constants
#endif // SLMC_SLMC_CFG_INCLUDE_CONSTANTS_HPP_
