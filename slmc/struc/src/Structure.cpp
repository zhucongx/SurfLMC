#include "Structure.h"
namespace cfg {
std::pair<std::vector<Lattice>, std::vector<Atom> > CreateOneLayerA(size_t n_x, size_t n_y, size_t layer_index) {
  Matrix3d basis{{static_cast<double>(n_x) * constants::kHexagonalLatticeConstant, 0, 0},
                 {-static_cast<double>(n_y) * constants::kHexagonalLatticeConstant / 2,
                  static_cast<double>(n_y) * constants::kHexagonalLatticeConstant * std::sqrt(3) / 2, 0},
                 {0, 0, constants::kLatticeLayerDistance}};
  std::vector<Lattice> lattice_vector;
  lattice_vector.reserve(n_x * n_y * 2);
  std::vector<Atom> atom_vector;
  atom_vector.reserve(n_x * n_y * 2);
  for (size_t j = 0; j < n_y; ++j) {
    for (size_t i = 0; i < n_x; ++i) {
      Vector3d relative_position1{(static_cast<double>(i) - 1. / 3) / static_cast<double>(n_x),
                                  (static_cast<double>(j) + 1. / 3) / static_cast<double>(n_y),
                                  static_cast<double>(layer_index)}; // Ga
      lattice_vector.emplace_back(0, basis * relative_position1, relative_position1);
      atom_vector.emplace_back(0, "Ga");
      Vector3d relative_position2{(static_cast<double>(i) - 1. / 3) / static_cast<double>(n_x),
                                  (static_cast<double>(j) + 1. / 3) / static_cast<double>(n_y),
                                  static_cast<double>(layer_index) + 0.75}; // N
      lattice_vector.emplace_back(0, basis * relative_position2, relative_position2);
      atom_vector.emplace_back(0, "N");
    }
  }
  return {lattice_vector, atom_vector};
}
std::pair<std::vector<Lattice>, std::vector<Atom> > CreateOneLayerB(size_t n_x, size_t n_y, size_t layer_index) {
  Matrix3d basis{{static_cast<double>(n_x) * constants::kHexagonalLatticeConstant, 0, 0},
                 {-static_cast<double>(n_y) * constants::kHexagonalLatticeConstant / 2,
                  static_cast<double>(n_y) * constants::kHexagonalLatticeConstant * std::sqrt(3) / 2, 0},
                 {0, 0, constants::kLatticeLayerDistance}};
  std::vector<Lattice> lattice_vector;
  lattice_vector.reserve(n_x * n_y * 2);
  std::vector<Atom> atom_vector;
  atom_vector.reserve(n_x * n_y * 2);
  for (size_t j = 0; j < n_y; ++j) {
    for (size_t i = 0; i < n_x; ++i) {
      Vector3d relative_position1{(static_cast<double>(i) + 1. / 3) / static_cast<double>(n_x),
                                  (static_cast<double>(j) - 1. / 3) / static_cast<double>(n_y),
                                  static_cast<double>(layer_index)}; // Ga
      lattice_vector.emplace_back(0, basis * relative_position1, relative_position1);
      atom_vector.emplace_back(0, "Ga");
      Vector3d relative_position2{(static_cast<double>(i) + 1. / 3) / static_cast<double>(n_x),
                                  (static_cast<double>(j) - 1. / 3) / static_cast<double>(n_y),
                                  static_cast<double>(layer_index) + 0.75}; // N
      lattice_vector.emplace_back(0, basis * relative_position2, relative_position2);
      atom_vector.emplace_back(0, "N");
    }
  }
  return {lattice_vector, atom_vector};
}
std::pair<std::vector<Lattice>, std::vector<Atom> > CreateOneLayerC(size_t n_x, size_t n_y, size_t layer_index) {
  Matrix3d basis{{static_cast<double>(n_x) * constants::kHexagonalLatticeConstant, 0, 0},
                 {-static_cast<double>(n_y) * constants::kHexagonalLatticeConstant / 2,
                  static_cast<double>(n_y) * constants::kHexagonalLatticeConstant * std::sqrt(3) / 2, 0},
                 {0, 0, constants::kLatticeLayerDistance}};
  std::vector<Lattice> lattice_vector;
  lattice_vector.reserve(n_x * n_y * 2);
  std::vector<Atom> atom_vector;
  atom_vector.reserve(n_x * n_y * 2);
  for (size_t j = 0; j < n_y; ++j) {
    for (size_t i = 0; i < n_x; ++i) {
      Vector3d relative_position1{(static_cast<double>(i)) / static_cast<double>(n_x),
                                  (static_cast<double>(j)) / static_cast<double>(n_y),
                                  static_cast<double>(layer_index)}; // Ga
      lattice_vector.emplace_back(0, basis * relative_position1, relative_position1);
      atom_vector.emplace_back(0, "Ga");
      Vector3d relative_position2{(static_cast<double>(i)) / static_cast<double>(n_x),
                                  (static_cast<double>(j)) / static_cast<double>(n_y),
                                  static_cast<double>(layer_index) + 0.75}; // N
      lattice_vector.emplace_back(0, basis * relative_position2, relative_position2);
      atom_vector.emplace_back(0, "N");
    }
  }
  return {lattice_vector, atom_vector};
}
cfg::Config CreateLayers(size_t n_x, size_t ny, const std::string &layers_type) {
  size_t n_layers = layers_type.size();
  Matrix3d basis{{static_cast<double>(n_x) * constants::kHexagonalLatticeConstant, 0, 0},
                 {-static_cast<double>(ny) * constants::kHexagonalLatticeConstant / 2,
                  static_cast<double>(ny) * constants::kHexagonalLatticeConstant * std::sqrt(3) / 2, 0},
                 {0, 0, static_cast<double>(n_layers) * constants::kLatticeLayerDistance}};
  auto inverse_basis = basis.inverse();
  std::vector<Lattice> lattice_vector;
  lattice_vector.reserve(n_x * ny * 2 * n_layers);
  std::vector<Atom> atom_vector;
  atom_vector.reserve(n_x * ny * 2 * n_layers);
  for (size_t layer_index = 0; layer_index < n_layers; ++layer_index) {
    switch (layers_type[layer_index]) {
      case 'A': {
        auto [lattice, atom] = CreateOneLayerA(n_x, ny, layer_index);
        lattice_vector.insert(lattice_vector.end(), lattice.begin(), lattice.end());
        atom_vector.insert(atom_vector.end(), atom.begin(), atom.end());
        break;
      }
      case 'B': {
        auto [lattice, atom] = CreateOneLayerB(n_x, ny, layer_index);
        lattice_vector.insert(lattice_vector.end(), lattice.begin(), lattice.end());
        atom_vector.insert(atom_vector.end(), atom.begin(), atom.end());
        break;
      }
      case 'C': {
        auto [lattice, atom] = CreateOneLayerC(n_x, ny, layer_index);
        lattice_vector.insert(lattice_vector.end(), lattice.begin(), lattice.end());
        atom_vector.insert(atom_vector.end(), atom.begin(), atom.end());
        break;
      }
      default: {
        throw std::runtime_error("Unknown layer type");
      }
    }
  }
  for (size_t i = 0; i < lattice_vector.size(); ++i) {
    lattice_vector[i].SetId(i);
    lattice_vector[i].SetRelativePosition(inverse_basis * lattice_vector[i].GetCartesianPosition()
                                              + Vector3d{0, 0, 1 / (static_cast<double>(n_layers))});
  }
  for (size_t i = 0; i < atom_vector.size(); ++i) {
    atom_vector[i].SetId(i);
  }
  return Config{basis, lattice_vector, atom_vector, false};
}
// cfg::Config CreateLayers(size_t n_x, size_t ny, const std::string &layers_type) {
//   size_t n_layers = layers_type.size();
//   Matrix3d basis{{static_cast<double>(n_x) * constants::kHexagonalLatticeConstant, 0, 0},
//                  {-static_cast<double>(ny) * constants::kHexagonalLatticeConstant / 2,
//                   static_cast<double>(ny) * constants::kHexagonalLatticeConstant *  std::sqrt(3) / 2, 0},
//                  {0, 0, static_cast<double>(n_layers + 2) * constants::kLatticeLayerDistance}};
//   auto inverse_basis = basis.inverse();
//   std::vector<Lattice> lattice_vector;
//   lattice_vector.reserve(n_x * ny * 2 * n_layers);
//   std::vector<Atom> atom_vector;
//   atom_vector.reserve(n_x * ny * 2 * n_layers);
//   for (size_t layer_index = 0; layer_index < n_layers; ++layer_index) {
//     switch (layers_type[layer_index]) {
//       case 'A': {
//         auto [lattice, atom] = CreateOneLayerA(n_x, ny, layer_index);
//         lattice_vector.insert(lattice_vector.end(), lattice.begin(), lattice.end());
//         atom_vector.insert(atom_vector.end(), atom.begin(), atom.end());
//         break;
//       }
//       case 'B': {
//         auto [lattice, atom] = CreateOneLayerB(n_x, ny, layer_index);
//         lattice_vector.insert(lattice_vector.end(), lattice.begin(), lattice.end());
//         atom_vector.insert(atom_vector.end(), atom.begin(), atom.end());
//         break;
//       }
//       case 'C': {
//         auto [lattice, atom] = CreateOneLayerC(n_x, ny, layer_index);
//         lattice_vector.insert(lattice_vector.end(), lattice.begin(), lattice.end());
//         atom_vector.insert(atom_vector.end(), atom.begin(), atom.end());
//         break;
//       }
//       default: {
//         throw std::runtime_error("Unknown layer type");
//       }
//     }
//   }
//   for (size_t i = 0; i < lattice_vector.size(); ++i) {
//     lattice_vector[i].SetId(i);
//     lattice_vector[i].SetRelativePosition(inverse_basis * lattice_vector[i].GetCartesianPosition()
//                                               + Vector3d{0, 0, 1 / (static_cast<double>(n_layers + 2))});
//   }
//   for (size_t i = 0; i < atom_vector.size(); ++i) {
//     atom_vector[i].SetId(i);
//   }
//   return Config{basis, lattice_vector, atom_vector, false};
// }
std::pair<std::vector<Lattice>, std::vector<Atom> > CreateOneOrthLayerA(size_t n_x, size_t n_y, size_t layer_index) {
  Matrix3d basis{{static_cast<double>(n_x) * constants::kHexagonalLatticeConstant * 2, 0, 0},
                 {0, static_cast<double>(n_y) * constants::kHexagonalLatticeConstant * std::sqrt(3), 0},
                 {0, 0, constants::kLatticeLayerDistance}};
  std::vector<Lattice> lattice_vector;
  lattice_vector.reserve(n_x * n_y * 2);
  std::vector<Atom> atom_vector;
  atom_vector.reserve(n_x * n_y * 2);
  for (size_t j = 0; j < n_y; ++j) {
    for (size_t i = 0; i < n_x; ++i) {
      const Vector3d relative_position{(static_cast<double>(i)) / static_cast<double>(n_x),
                                       (static_cast<double>(j)) / static_cast<double>(n_y),
                                       static_cast<double>(layer_index)};
      std::vector<Vector3d> position_list = {{0, 0, 0},
                                             {0.5, 0, 0},
                                             {0.25, 0.5, 0},
                                             {0.75, 0.5, 0}};
      for (const auto &position: position_list) {
        auto relative_position1 = relative_position + position;
        lattice_vector.emplace_back(0, basis * relative_position1, relative_position1);
        atom_vector.emplace_back(0, "Ga");
        auto relative_position2 = relative_position1 + Vector3d{0, 0, 0.75};
        lattice_vector.emplace_back(0, basis * relative_position2, relative_position2);
        atom_vector.emplace_back(0, "N");
      }
    }
  }
  return {lattice_vector, atom_vector};
}
std::pair<std::vector<Lattice>, std::vector<Atom> > CreateOneOrthLayerB(size_t n_x, size_t n_y, size_t layer_index) {
  Matrix3d basis{{static_cast<double>(n_x) * constants::kHexagonalLatticeConstant * 2, 0, 0},
                 {0, static_cast<double>(n_y) * constants::kHexagonalLatticeConstant * std::sqrt(3), 0},
                 {0, 0, constants::kLatticeLayerDistance}};
  std::vector<Lattice> lattice_vector;
  lattice_vector.reserve(n_x * n_y * 2);
  std::vector<Atom> atom_vector;
  atom_vector.reserve(n_x * n_y * 2);
  for (size_t j = 0; j < n_y; ++j) {
    for (size_t i = 0; i < n_x; ++i) {
      const Vector3d relative_position{(static_cast<double>(i)) / static_cast<double>(n_x),
                                       (static_cast<double>(j) + 1. / 3) / static_cast<double>(n_y),
                                       static_cast<double>(layer_index)};
      std::vector<Vector3d> position_list = {{0, 0, 0},
                                             {0.5, 0, 0},
                                             {0.25, 0.5, 0},
                                             {0.75, 0.5, 0}};
      for (const auto &position: position_list) {
        auto relative_position1 = relative_position + position;
        lattice_vector.emplace_back(0, basis * relative_position1, relative_position1);
        atom_vector.emplace_back(0, "Ga");
        auto relative_position2 = relative_position1 + Vector3d{0, 0, 0.75};
        lattice_vector.emplace_back(0, basis * relative_position2, relative_position2);
        atom_vector.emplace_back(0, "N");
      }
    }
  }
  return {lattice_vector, atom_vector};
}
std::pair<std::vector<Lattice>, std::vector<Atom> > CreateOneOrthLayerC(size_t n_x, size_t n_y, size_t layer_index) {
  Matrix3d basis{{static_cast<double>(n_x) * constants::kHexagonalLatticeConstant * 2, 0, 0},
                 {0, static_cast<double>(n_y) * constants::kHexagonalLatticeConstant * std::sqrt(3), 0},
                 {0, 0, constants::kLatticeLayerDistance}};
  std::vector<Lattice> lattice_vector;
  lattice_vector.reserve(n_x * n_y * 2);
  std::vector<Atom> atom_vector;
  atom_vector.reserve(n_x * n_y * 2);
  for (size_t j = 0; j < n_y; ++j) {
    for (size_t i = 0; i < n_x; ++i) {
      const Vector3d relative_position{(static_cast<double>(i) + 1. / 4) / static_cast<double>(n_x),
                                       (static_cast<double>(j) + 1. / 6) / static_cast<double>(n_y),
                                       static_cast<double>(layer_index)};
      std::vector<Vector3d> position_list = {{0, 0, 0},
                                             {0.5, 0, 0},
                                             {0.25, 0.5, 0},
                                             {0.75, 0.5, 0}};
      for (const auto &position: position_list) {
        auto relative_position1 = relative_position + position;
        lattice_vector.emplace_back(0, basis * relative_position1, relative_position1);
        atom_vector.emplace_back(0, "Ga");
        auto relative_position2 = relative_position1 + Vector3d{0, 0, 0.75};
        lattice_vector.emplace_back(0, basis * relative_position2, relative_position2);
        atom_vector.emplace_back(0, "N");
      }
    }
  }
  return {lattice_vector, atom_vector};
}

cfg::Config CreateOrthLayers(size_t n_x, size_t ny, const std::string &layers_type) {
  size_t n_layers = layers_type.size();
  Matrix3d basis{{static_cast<double>(n_x) * constants::kHexagonalLatticeConstant * 2, 0, 0},
                 {0, static_cast<double>(ny) * constants::kHexagonalLatticeConstant * std::sqrt(3), 0},
                 {0, 0, static_cast<double>(n_layers) * constants::kLatticeLayerDistance}};
  auto inverse_basis = basis.inverse();
  std::vector<Lattice> lattice_vector;
  lattice_vector.reserve(n_x * ny * 2 * n_layers);
  std::vector<Atom> atom_vector;
  atom_vector.reserve(n_x * ny * 2 * n_layers);
  for (size_t layer_index = 0; layer_index < n_layers; ++layer_index) {
    switch (layers_type[layer_index]) {
      case 'A': {
        auto [lattice, atom] = CreateOneOrthLayerA(n_x, ny, layer_index);
        lattice_vector.insert(lattice_vector.end(), lattice.begin(), lattice.end());
        atom_vector.insert(atom_vector.end(), atom.begin(), atom.end());
        break;
      }
      case 'B': {
        auto [lattice, atom] = CreateOneOrthLayerB(n_x, ny, layer_index);
        lattice_vector.insert(lattice_vector.end(), lattice.begin(), lattice.end());
        atom_vector.insert(atom_vector.end(), atom.begin(), atom.end());
        break;
      }
      case 'C': {
        auto [lattice, atom] = CreateOneOrthLayerC(n_x, ny, layer_index);
        lattice_vector.insert(lattice_vector.end(), lattice.begin(), lattice.end());
        atom_vector.insert(atom_vector.end(), atom.begin(), atom.end());
        break;
      }
      default: {
        throw std::runtime_error("Unknown layer type");
      }
    }
  }
  for (size_t i = 0; i < lattice_vector.size(); ++i) {
    lattice_vector[i].SetId(i);
    lattice_vector[i].SetRelativePosition(inverse_basis * lattice_vector[i].GetCartesianPosition());
  }
  for (size_t i = 0; i < atom_vector.size(); ++i) {
    atom_vector[i].SetId(i);
  }
  return Config{basis, lattice_vector, atom_vector, false};
}

} // cfg