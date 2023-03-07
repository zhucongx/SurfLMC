#ifndefSLMC_SLMC_CFG_INCLUDE_ELEMENT_HPP_
#defineSLMC_SLMC_CFG_INCLUDE_ELEMENT_HPP_
#include <cmath>
#include <string>
#include <map>
#include <iostream>

enum class ElementName {
  X, B, Al, Ga, In, N, P, As, Sb
};

class Element {
  public:
    Element() = default;
    constexpr explicit Element(ElementName element_type) : element_name_(element_type) {}
    explicit Element(const std::string &element_string) {
      const std::map<std::string, ElementName> ElementStrings{
          {"X", ElementName::X},
          {"B", ElementName::B},
          {"Al", ElementName::Al},
          {"Ga", ElementName::Ga},
          {"In", ElementName::In},
          {"N", ElementName::N},
          {"P", ElementName::P},
          {"As", ElementName::As},
          {"Sb", ElementName::Sb}

      };
      auto it = ElementStrings.find(element_string);
      element_name_ = it == ElementStrings.end() ? ElementName::X : it->second;
    }

    constexpr explicit operator ElementName() const { return element_name_; }
    explicit operator bool() = delete;
    constexpr bool operator==(Element rhs) const {
      return element_name_ == rhs.element_name_;
    }
    constexpr bool operator!=(Element rhs) const {
      return element_name_ != rhs.element_name_;
    }
    constexpr bool operator==(ElementName rhs) const {
      return element_name_ == rhs;
    }
    constexpr bool operator!=(ElementName rhs) const {
      return element_name_ != rhs;
    }
    bool operator<(Element rhs) const {
      return GetString() < rhs.GetString();
    }
    bool operator<(const ElementName rhs) const {
      return GetString() < Element(rhs).GetString();
    }
    friend size_t hash_value(Element element) {
      return static_cast<std::size_t>(element.element_name_);
    }
    [[nodiscard]] std::string GetString() const {
      switch (element_name_) {
        case ElementName::X : return "X";
        case ElementName::B: return "B";
        case ElementName::Al: return "Al";
        case ElementName::Ga: return "Ga";
        case ElementName::In: return "In";
        case ElementName::N: return "N";
        case ElementName::P: return "P";
        case ElementName::As: return "As";
        case ElementName::Sb: return "Sb";
        default: throw std::invalid_argument("Unexpected pseudo element");
          // omit default case to trigger compiler warning for missing cases
      }
    }
    [[nodiscard]] double GetMass() const {
      switch (element_name_) {
        case ElementName::X : return 0.000;
        case ElementName::B: return 10.811;
        case ElementName::Al: return 26.982;
        case ElementName::Ga: return 69.723;
        case ElementName::In: return 114.818;
        case ElementName::N: return 14.007;
        case ElementName::P: return 30.974;
        case ElementName::As: return 74.922;
        case ElementName::Sb: return 121.760;
        default: throw std::invalid_argument("Unexpected pseudo element");
          // omit default case to trigger compiler warning for missing cases
      }
    }
  private:
    ElementName element_name_{};
};

#endif //LMC_LMC_CFG_INCLUDE_ELEMENT_HPP_
