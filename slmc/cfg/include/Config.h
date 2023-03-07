#ifndef SLMC_SLMC_CFG_INCLUDE_CONFIG_H_
#define SLMC_SLMC_CFG_INCLUDE_CONFIG_H_
#include <vector>
#include <set>
#include <map>
#include <unordered_set>
#include <unordered_map>
#include "Atom.hpp"
#include "Lattice.hpp"
// using Graph = boost::adjacency_list<boost::vecS, boost::vecS>;
using Eigen::Matrix3Xd;
namespace cfg {
class Config {
  public:
    /// Constructor
    Config();
    Config(Matrix3d basis,
           std::vector<Lattice> lattice_vector,
           std::vector<Atom> atom_vector,
           bool update_neighbor);
    /// Getter
    [[nodiscard]] size_t GetNumAtoms() const;
    [[nodiscard]] const Matrix3d &GetBasis() const;
    [[nodiscard]] const Matrix3Xd &GetRelativePositionMatrix() const;
    [[nodiscard]] const Matrix3Xd &GetCartesianPositionMatrix() const;
    [[nodiscard]] const std::vector<Atom> &GetAtomVector() const;
    [[nodiscard]] const std::vector<std::vector<size_t> > &GetFirstNeighborsAdjacencyList() const;
    [[nodiscard]] std::vector<size_t> GetFirstNeighborsAtomIdVectorOfAtom(size_t atom_id) const;
    [[nodiscard]] size_t GetAtomIdFromLatticeId(size_t lattice_id) const;
    [[nodiscard]] size_t GetLatticeIdFromAtomId(size_t atom_id) const;
    [[nodiscard]] Element GetElementAtAtomId(size_t atom_id) const;
    [[nodiscard]] Element GetElementAtLatticeId(size_t lattice_id) const;
    [[nodiscard]] std::set<Element> GetElementSetWithoutVacancy() const;
    [[nodiscard]] std::map<Element, std::vector<size_t> > GetElementAtomIdVectorMap() const;
    [[nodiscard]] size_t GetStateHash() const;
    /// Modify config
    void AtomJump(const std::pair<size_t, size_t> &atom_id_jump_pair);
    void LatticeJump(const std::pair<size_t, size_t> &lattice_id_jump_pair);
    void ChangeAtomElementTypeAtAtom(size_t atom_id, Element element);
    void ChangeAtomElementTypeAtLattice(size_t lattice_id, Element element);

    /// IO
    static Config ReadConfig(const std::string &filename);
    void WriteConfig(const std::string &filename) const;
    void WriteExtendedConfig(
        const std::string &filename,
        const std::map<std::string, std::vector<double> >& auxiliary_lists) const;
  private:
    /// Private Getter
    [[nodiscard]] const std::unordered_map<size_t, size_t> &GetLatticeToAtomHashmap() const;
    [[nodiscard]] const std::unordered_map<size_t, size_t> &GetAtomToLatticeHashmap() const;
    /// Modify config

    void InitializeNeighborsList(size_t num_atoms);
    void UpdateNeighbors();
    /// Properties
    Matrix3d basis_{};
    Matrix3Xd cartesian_position_matrix_{};
    Matrix3Xd relative_position_matrix_{};
    std::vector<Atom> atom_vector_{};
    std::unordered_map<size_t, size_t> lattice_to_atom_hashmap_{};
    std::unordered_map<size_t, size_t> atom_to_lattice_hashmap_{};
    // nearest neighbor lists
    std::vector<std::vector<size_t> > first_neighbors_adjacency_list_{};

};

} // cfg
#endif // SLMC_SLMC_CFG_INCLUDE_CONFIG_H_
