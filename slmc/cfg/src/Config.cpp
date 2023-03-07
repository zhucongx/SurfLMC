#include "Config.h"

#include <random>
#include <chrono>
#include <utility>
#include <algorithm>
#include <boost/functional/hash.hpp>

#include <omp.h>

namespace cfg {
Config::Config() = default;
Config::Config(Matrix3d basis,
               std::vector<Lattice> lattice_vector,
               std::vector<Atom> atom_vector,
               bool update_neighbor)
    : basis_(std::move(basis)),
      atom_vector_(std::move(atom_vector)) {

  if (lattice_vector.size() != atom_vector_.size()) {
    throw std::runtime_error("Lattice vector and atom vector size do not match");
  }
  relative_position_matrix_ = Matrix3Xd(3, static_cast<int>(lattice_vector.size()));
  cartesian_position_matrix_ = Matrix3Xd(3, static_cast<int>(lattice_vector.size()));
  for (size_t i = 0; i < lattice_vector.size(); ++i) {
    auto lattice_id = lattice_vector.at(i).GetId();
    auto atom_id = atom_vector_.at(i).GetId();
    lattice_to_atom_hashmap_.emplace(lattice_id, atom_id);
    atom_to_lattice_hashmap_.emplace(atom_id, lattice_id);
    relative_position_matrix_.col(static_cast<int>(i)) = lattice_vector.at(i).GetRelativePosition();
    cartesian_position_matrix_.col(static_cast<int>(i)) = lattice_vector.at(i).GetCartesianPosition();
  }
  if (update_neighbor) {
    UpdateNeighbors();
  }
}
size_t Config::GetNumAtoms() const {
  return atom_vector_.size();
}
const Matrix3d &Config::GetBasis() const {
  return basis_;
}
const std::unordered_map<size_t, size_t> &Config::GetLatticeToAtomHashmap() const {
  return lattice_to_atom_hashmap_;
}
const std::unordered_map<size_t, size_t> &Config::GetAtomToLatticeHashmap() const {
  return atom_to_lattice_hashmap_;
}
const Matrix3Xd &Config::GetRelativePositionMatrix() const {
  return relative_position_matrix_;
}
const Matrix3Xd &Config::GetCartesianPositionMatrix() const {
  return cartesian_position_matrix_;
}
const std::vector<Atom> &Config::GetAtomVector() const {
  return atom_vector_;
}
const std::vector<std::vector<size_t> > &Config::GetFirstNeighborsAdjacencyList() const {
  return first_neighbors_adjacency_list_;
}
std::vector<size_t> Config::GetFirstNeighborsAtomIdVectorOfAtom(size_t atom_id) const {
  auto lattice_id = atom_to_lattice_hashmap_.at(atom_id);
  std::vector<size_t> first_neighbors_atom_id_vector;
  first_neighbors_atom_id_vector.reserve(constants::kNumFirstNearestNeighbors);
  for (auto neighbor_lattice_id: first_neighbors_adjacency_list_[lattice_id]) {
    first_neighbors_atom_id_vector.push_back(lattice_to_atom_hashmap_.at(neighbor_lattice_id));
  }
  return first_neighbors_atom_id_vector;
}
size_t Config::GetAtomIdFromLatticeId(size_t lattice_id) const {
  return lattice_to_atom_hashmap_.at(lattice_id);
}
size_t Config::GetLatticeIdFromAtomId(size_t atom_id) const {
  return atom_to_lattice_hashmap_.at(atom_id);
}
Element Config::GetElementAtAtomId(size_t atom_id) const {
  return atom_vector_[atom_id].GetElement();
}
Element Config::GetElementAtLatticeId(size_t lattice_id) const {
  auto atom_id = lattice_to_atom_hashmap_.at(lattice_id);
  return atom_vector_[atom_id].GetElement();
}
std::set<Element> Config::GetElementSetWithoutVacancy() const {
  std::set<Element> res;
  for (const auto &atom: atom_vector_) {
    if (atom.GetElement() == ElementName::X) { continue; }
    res.insert(atom.GetElement());
  }
  return res;
}
std::map<Element, std::vector<size_t> > Config::GetElementAtomIdVectorMap() const {
  std::map<Element, std::vector<size_t> > element_list_map;
  for (const auto &atom: atom_vector_) {
    element_list_map[atom.GetElement()].push_back(atom.GetId());
  }
  return element_list_map;
}
size_t Config::GetStateHash() const {
  size_t seed = 0;
  for (size_t i = 0; i < GetNumAtoms(); ++i) {
    boost::hash_combine(seed, lattice_to_atom_hashmap_.at(i));
  }
  return seed;
}

void Config::AtomJump(const std::pair<size_t, size_t> &atom_id_jump_pair) {
  const auto [atom_id_lhs, atom_id_rhs] = atom_id_jump_pair;
  const auto lattice_id_lhs = atom_to_lattice_hashmap_.at(atom_id_lhs);
  const auto lattice_id_rhs = atom_to_lattice_hashmap_.at(atom_id_rhs);

  atom_to_lattice_hashmap_.at(atom_id_lhs) = lattice_id_rhs;
  atom_to_lattice_hashmap_.at(atom_id_rhs) = lattice_id_lhs;
  lattice_to_atom_hashmap_.at(lattice_id_lhs) = atom_id_rhs;
  lattice_to_atom_hashmap_.at(lattice_id_rhs) = atom_id_lhs;
}
void Config::LatticeJump(const std::pair<size_t, size_t> &lattice_id_jump_pair) {
  const auto [lattice_id_lhs, lattice_id_rhs] = lattice_id_jump_pair;
  const auto atom_id_lhs = lattice_to_atom_hashmap_.at(lattice_id_lhs);
  const auto atom_id_rhs = lattice_to_atom_hashmap_.at(lattice_id_rhs);

  atom_to_lattice_hashmap_.at(atom_id_lhs) = lattice_id_rhs;
  atom_to_lattice_hashmap_.at(atom_id_rhs) = lattice_id_lhs;
  lattice_to_atom_hashmap_.at(lattice_id_lhs) = atom_id_rhs;
  lattice_to_atom_hashmap_.at(lattice_id_rhs) = atom_id_lhs;
}

void Config::ChangeAtomElementTypeAtAtom(size_t atom_id, Element element) {
  atom_vector_.at(atom_id).SetElement(element);
}
void Config::ChangeAtomElementTypeAtLattice(size_t lattice_id, Element element) {
  atom_vector_.at(lattice_to_atom_hashmap_.at(lattice_id)).SetElement(element);
}
Config Config::ReadConfig(const std::string &filename) {
  std::ifstream ifs(filename, std::ifstream::in);
  if (!ifs.is_open()) {
    throw std::runtime_error("Cannot open " + filename);
  }
  // "Number of particles = %i"
  ifs.ignore(std::numeric_limits<std::streamsize>::max(), '=');
  size_t num_atoms;
  ifs >> num_atoms;
  // A = 1.0 Angstrom (basic length-scale)
  ifs.ignore(std::numeric_limits<std::streamsize>::max(), '=');
  double basis_xx, basis_xy, basis_xz,
      basis_yx, basis_yy, basis_yz,
      basis_zx, basis_zy, basis_zz;
  // "H0(1,1) = %lf A"
  ifs.ignore(std::numeric_limits<std::streamsize>::max(), '=');
  ifs >> basis_xx;
  // "H0(1,2) = %lf A"
  ifs.ignore(std::numeric_limits<std::streamsize>::max(), '=');
  ifs >> basis_xy;
  // "H0(1,3) = %lf A"
  ifs.ignore(std::numeric_limits<std::streamsize>::max(), '=');
  ifs >> basis_xz;
  // "H0(2,1) = %lf A"
  ifs.ignore(std::numeric_limits<std::streamsize>::max(), '=');
  ifs >> basis_yx;
  // "H0(2,2) = %lf A"
  ifs.ignore(std::numeric_limits<std::streamsize>::max(), '=');
  ifs >> basis_yy;
  // "H0(2,3) = %lf A"
  ifs.ignore(std::numeric_limits<std::streamsize>::max(), '=');
  ifs >> basis_yz;
  // "H0(3,1) = %lf A"
  ifs.ignore(std::numeric_limits<std::streamsize>::max(), '=');
  ifs >> basis_zx;
  // "H0(3,2) = %lf A"
  ifs.ignore(std::numeric_limits<std::streamsize>::max(), '=');
  ifs >> basis_zy;
  // "H0(3,3) = %lf A"
  ifs.ignore(std::numeric_limits<std::streamsize>::max(), '=');
  ifs >> basis_zz;
  // finish this line
  ifs.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  // .NO_VELOCITY.
  ifs.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  // "entry_count = 3"
  ifs.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  auto basis = Matrix3d{{basis_xx, basis_xy, basis_xz},
                        {basis_yx, basis_yy, basis_yz},
                        {basis_zx, basis_zy, basis_zz}};
  std::vector<Atom> atom_vector;
  atom_vector.reserve(num_atoms);
  std::vector<Lattice> lattice_vector;
  lattice_vector.reserve(num_atoms);

  double mass;
  std::string type;
  Vector3d relative_position;

  std::vector<std::vector<size_t> > first_neighbors_adjacency_list,
      second_neighbors_adjacency_list, third_neighbors_adjacency_list;

  for (size_t lattice_id = 0; lattice_id < num_atoms; ++lattice_id) {
    ifs >> mass >> type >> relative_position(0) >> relative_position(1) >> relative_position(2);
    atom_vector.emplace_back(lattice_id, type);
    lattice_vector.emplace_back(lattice_id, basis * relative_position, relative_position);
  }
  return Config{basis, lattice_vector, atom_vector, true};
}
void Config::WriteConfig(const std::string &filename) const {
  Eigen::IOFormat fmt(Eigen::StreamPrecision, Eigen::DontAlignCols, " ", " ", "", "", "", "");
  std::ofstream ofs(filename, std::ofstream::out);
  ofs.precision(16);
  ofs << "Number of particles = " << GetNumAtoms() << '\n';
  ofs << "A = 1.0 Angstrom (basic length-scale)\n";
  ofs << "H0(1,1) = " << basis_(0, 0) << " A\n";
  ofs << "H0(1,2) = " << basis_(0, 1) << " A\n";
  ofs << "H0(1,3) = " << basis_(0, 2) << " A\n";
  ofs << "H0(2,1) = " << basis_(1, 0) << " A\n";
  ofs << "H0(2,2) = " << basis_(1, 1) << " A\n";
  ofs << "H0(2,3) = " << basis_(1, 2) << " A\n";
  ofs << "H0(3,1) = " << basis_(2, 0) << " A\n";
  ofs << "H0(3,2) = " << basis_(2, 1) << " A\n";
  ofs << "H0(3,3) = " << basis_(2, 2) << " A\n";
  ofs << ".NO_VELOCITY.\n";
  ofs << "entry_count = 3\n";
  for (const auto &atom: atom_vector_) {
    size_t lattice_id = atom_to_lattice_hashmap_.at(atom.GetId());
    const auto &relative_position = relative_position_matrix_.col(static_cast<int>(lattice_id));
    ofs << atom.GetMass() << '\n'
        << atom.GetElementString() << '\n'
        << relative_position.transpose().format(fmt) << std::endl;
  }
}
void Config::WriteExtendedConfig(
    const std::string &filename,
    const std::map<std::string, std::vector<double> > &auxiliary_lists) const {
  Eigen::IOFormat fmt(Eigen::StreamPrecision, Eigen::DontAlignCols, " ", " ", "", "", "", "");
  std::ofstream ofs(filename, std::ofstream::out);
  ofs.precision(16);
  ofs << "Number of particles = " << GetNumAtoms() << '\n';
  ofs << "A = 1.0 Angstrom (basic length-scale)\n";
  ofs << "H0(1,1) = " << basis_(0, 0) << " A\n";
  ofs << "H0(1,2) = " << basis_(0, 1) << " A\n";
  ofs << "H0(1,3) = " << basis_(0, 2) << " A\n";
  ofs << "H0(2,1) = " << basis_(1, 0) << " A\n";
  ofs << "H0(2,2) = " << basis_(1, 1) << " A\n";
  ofs << "H0(2,3) = " << basis_(1, 2) << " A\n";
  ofs << "H0(3,1) = " << basis_(2, 0) << " A\n";
  ofs << "H0(3,2) = " << basis_(2, 1) << " A\n";
  ofs << "H0(3,3) = " << basis_(2, 2) << " A\n";
  ofs << ".NO_VELOCITY.\n";
  ofs << "entry_count = " << 3 + auxiliary_lists.size() << "\n";
  size_t auxiliary_index = 0;
  for (const auto &auxiliary_list: auxiliary_lists) {
    ofs << "auxiliary[" << auxiliary_index << "] = " << auxiliary_list.first << " [reduced unit]\n";
    ++auxiliary_index;
  }
  for (size_t it = 0; it < atom_vector_.size(); ++it) {
    const auto &atom = atom_vector_[it];
    const auto &relative_position
        = relative_position_matrix_.col(static_cast<int>(atom_to_lattice_hashmap_.at(atom.GetId())));
    ofs << atom.GetMass() << '\n'
        << atom.GetElementString() << '\n'
        << relative_position.transpose().format(fmt);
    for (const auto &auxiliary_list: auxiliary_lists) {
      ofs << ' ' << auxiliary_list.second[it];
    }
    ofs << std::endl;
  }
}

void Config::InitializeNeighborsList(size_t num_atoms) {
  first_neighbors_adjacency_list_.resize(num_atoms);
  for (auto &neighbor_list: first_neighbors_adjacency_list_) {
    neighbor_list.clear();
    neighbor_list.reserve(constants::kNumFirstNearestNeighbors);
  }
}
void Config::UpdateNeighbors() {

}


// Config GenerateClusteredConfigFromExcitingPure(Config config,
//                                                const std::map<Element, size_t> &solute_atom_count) {
//   return Config();
// }
// Config GenerateClusteredConfig(const Factor_t &factors,
//                                Element solvent_element,
//                                const std::map<Element, size_t> &solute_atom_count) {
//   return Config();
// }


} // cfg
