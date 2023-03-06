#ifndef LMC_LMC_CFG_INCLUDE_ELEMENTCLUSTER_HPP_
#define LMC_LMC_CFG_INCLUDE_ELEMENTCLUSTER_HPP_
#include <utility>
#include <vector>
#include <map>
#include <set>
#include <string>
#include <unordered_set>
#include <unordered_map>
#include <type_traits>
#include <boost/functional/hash.hpp>
#include "Element.hpp"
namespace cfg {

class ElementCluster {
  public:
    /// Constructor
    ElementCluster() = default;
    ElementCluster(int label, std::vector<Element> element_vector)
        : label_(label), element_vector_(std::move(element_vector)) {
      std::sort(element_vector_.begin(), element_vector_.end());
    }
    template<typename ... Ts>
    explicit ElementCluster(int label, Ts &&... ts)
        : label_(label), element_vector_{std::forward<Ts>(ts)...} {
      std::sort(element_vector_.begin(), element_vector_.end());
    }
    /// Getter
    [[nodiscard]] size_t GetSize() const {
      return element_vector_.size();
    }
    [[nodiscard]] int GetLabel() const {
      return label_;
    }
    [[nodiscard]] const std::vector<Element> &GetElementVector() const {
      return element_vector_;
    }
    /// Operators
    friend bool operator<(const ElementCluster &lhs, const ElementCluster &rhs) {
      if (lhs.element_vector_.size() < rhs.element_vector_.size()) { return true; }
      if (rhs.element_vector_.size() < lhs.element_vector_.size()) { return false; }
      if (lhs.label_ < rhs.label_) { return true; }
      if (rhs.label_ < lhs.label_) { return false; }
      for (size_t i = 0; i < lhs.element_vector_.size(); ++i) {
        if (lhs.element_vector_[i] < rhs.element_vector_[i]) { return true; }
        if (rhs.element_vector_[i] < lhs.element_vector_[i]) { return false; }
      }
      return false;
    }
    friend bool operator==(const ElementCluster &lhs,
                           const ElementCluster &rhs) {
      if (lhs.element_vector_.size() != rhs.element_vector_.size()) { return false; }
      if (lhs.label_ != rhs.label_) { return false; }
      for (size_t i = 0; i < lhs.element_vector_.size(); ++i) {
        if (lhs.element_vector_[i] != rhs.element_vector_[i])
          return false;
      }
      return true;
    }
    friend std::ostream &operator<<(std::ostream &os, const ElementCluster &element_cluster) {
      os << element_cluster.label_ << ' ';
      for (auto element: element_cluster.element_vector_) {
        os << element.GetString() << '-';
      }
      return os;
    }
    friend size_t hash_value(const ElementCluster &element_cluster) {
      size_t seed = 0;
      boost::hash_combine(seed, element_cluster.element_vector_.size());
      boost::hash_combine(seed, element_cluster.label_);
      for (auto element: element_cluster.element_vector_) {
        boost::hash_combine(seed, element);
      }
      return seed;
    }
  private:
    int label_{};
    std::vector<Element> element_vector_;
};

} // cfg
#endif //LMC_LMC_CFG_INCLUDE_ELEMENTCLUSTER_HPP_
