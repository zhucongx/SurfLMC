#ifndef SLMC_SLMC_CFG_INCLUDE_LATTICE_HPP_
#define SLMC_SLMC_CFG_INCLUDE_LATTICE_HPP_
#include "Constants.hpp"
#include <Eigen/Dense>
#include <utility>
using Eigen::Matrix3d;
using Eigen::Vector3d;
using Eigen::Vector3i;
namespace cfg {

class Lattice {
  public:
    /// Constructor
    Lattice() = default;
    Lattice(size_t id, const Vector3d &position)
        : id_(id),
          cartesian_position_(position),
          relative_position_(position) {}
    Lattice(size_t id,
            Vector3d cartesian_position,
            Vector3d relative_position)
        : id_(id),
          cartesian_position_(std::move(cartesian_position)),
          relative_position_(std::move(relative_position)) {}
    Lattice(size_t id, double x, double y, double z)
        : id_(id),
          cartesian_position_{x, y, z},
          relative_position_{x, y, z} {}
    /// Getter
    [[nodiscard]] size_t GetId() const {
      return id_;
    }
    [[nodiscard]] const Vector3d &GetCartesianPosition() const {
      return cartesian_position_;
    }
    [[nodiscard]] const Vector3d &GetRelativePosition() const {
      return relative_position_;
    }
    /// Setter
    void SetId(size_t id) {
      id_ = id;
    }
    void SetCartesianPosition(const Vector3d &cartesian_position) {
      cartesian_position_ = cartesian_position;
    }
    void SetRelativePosition(const Vector3d &relative_position) {
      relative_position_ = relative_position;
    }
  private:
    // lattice id
    size_t id_{};
    // absolute position
    Vector3d cartesian_position_{};
    // relative position in the box
    Vector3d relative_position_{};
};

inline Vector3d GetRelativeDistanceVectorLattice(const Lattice &first, const Lattice &second) {
  Vector3d relative_distance_vector = second.GetRelativePosition() - first.GetRelativePosition();
  // periodic boundary conditions
  for (const auto kDim: {0, 1, 2}) {
    while (relative_distance_vector[kDim] >= 0.5) { relative_distance_vector[kDim] -= 1; }
    while (relative_distance_vector[kDim] < -0.5) { relative_distance_vector[kDim] += 1; }
  }
  return relative_distance_vector;
}
} // cfg
#endif // SLMC_SLMC_CFG_INCLUDE_LATTICE_HPP_
