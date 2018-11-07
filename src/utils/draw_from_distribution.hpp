#pragma once

#include <random>
#include <vector>
#include "tlx/sort/strings/string_set.hpp"

using namespace tlx;
using namespace sort_strings_detail;

namespace dss {
  template <typename StringSet>
    std::vector<double> get_weights(const StringSet& string_set)
    {
      std::vector<double> weights;
      weights.reserve(string_set.size());
      for(const auto& str : string_set)
        weights.push_back(str);
      return weights;
    }

  template <typename WeightsContainer>
    std::vector<int> draw_elements_from_distribution(size_t number_of_elements, const WeightsContainer& weights)
    {
      std::random_device rd;
      std::mt19937 gen(rd());
      std::discrete_distribution<> d(weights.begin(), weights.end());
      std::vector<int> chosenElements;
      chosenElements.reserve(number_of_elements);
      for(size_t i = 0; i < number_of_elements; ++i)
        chosenElements.push_back(d(gen));
      return chosenElements;
    }
}
