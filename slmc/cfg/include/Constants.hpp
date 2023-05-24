#ifndef SLMC_SLMC_CFG_INCLUDE_CONSTANTS_HPP_
#define SLMC_SLMC_CFG_INCLUDE_CONSTANTS_HPP_

#include <cmath>
#include <array>

constexpr double kEpsilon = 1e-8;



namespace constants {
constexpr double kHexagonalLatticeConstant = 3.21; // Angstrom
constexpr double kLatticeLayerDistance = 2.62; // Angstrom

constexpr size_t kNumFirstNearestNeighbors = 4;
constexpr size_t kNumSecondNearestNeighbors = 12;


constexpr double kFirstNearestNeighborsCutoff = 2; // 1.95 Ga_N bond length
constexpr double kSecondNearestNeighborsCutoff = 3.2; // 3.18 Ga_Ga bond length


constexpr double kBoltzmann = 8.617333262145e-5; // eV/K
constexpr double kPrefactor = 1e13;
} // constants
#endif // SLMC_SLMC_CFG_INCLUDE_CONSTANTS_HPP_
