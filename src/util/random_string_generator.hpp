/*******************************************************************************
 * tests/util/random_string_generator.hpp
 *
 * Copyright (C) 2018 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include <cstdint>
#include <limits>
#include <random>

#include "util/indexed_string_set.hpp"
#include "util/string.hpp"
#include "util/string_set.hpp"

#include "strings/stringptr.hpp"
#include "strings/stringcontainer.hpp"

#include "mpi/environment.hpp"

namespace dss_schimek {
  using namespace dss_schimek;

  template <typename StringSet>
  class PrefixNumberStringLcpContainer : public StringLcpContainer<StringSet>
  {
    using Char = typename StringSet::Char;
    public:
    PrefixNumberStringLcpContainer(const size_t size, const Char prefix)
    {
      std::vector<Char> raw_string_data;
      for(size_t i = 1; i <= size; ++i)
      {
        raw_string_data.emplace_back(prefix);
        size_t number_to_print = i;
        while(number_to_print > 0)
        {
          Char last_digit = number_to_print % 10;
          raw_string_data.emplace_back(48 + last_digit);
          number_to_print /= 10;
        }
        raw_string_data.emplace_back(Char(0));
      }
      this->update(std::move(raw_string_data));
    }
  };

  template <typename StringSet>
  class RandomStringLcpContainer : public StringLcpContainer<StringSet>
  {
    using Char = typename StringSet::Char;
    public:
    RandomStringLcpContainer(const size_t size,
          const size_t min_length = 10,
          const size_t max_length = 20)
    {
      std::vector<Char> random_raw_string_data;
      std::random_device rand_seed;
      std::mt19937 rand_gen(rand_seed());
      std::uniform_int_distribution<Char> char_dis(65, 90);
      
      std::uniform_int_distribution<size_t> length_dis(min_length, max_length);
      random_raw_string_data.reserve(size + 1);
      for (size_t i = 0; i < size; ++i) {
        size_t length = length_dis(rand_gen);
        for (size_t j = 0; j < length; ++j)
          random_raw_string_data.emplace_back(char_dis(rand_gen));
        random_raw_string_data.emplace_back(Char(0));
      }
      this->update(std::move(random_raw_string_data));
    }
  };

  template <typename StringSet>
  class SkewedRandomStringLcpContainer : public StringLcpContainer<StringSet>
  {
    using Char = typename StringSet::Char;
    public:
    SkewedRandomStringLcpContainer(const size_t size,
          const size_t min_length = 10,
          const size_t max_length = 20)
    {
      std::vector<Char> random_raw_string_data;
      std::random_device rand_seed;
      std::mt19937 rand_gen(rand_seed());
      std::uniform_int_distribution<Char> small_char_dis(65, 70);
      std::uniform_int_distribution<Char> char_dis(65, 90);
      
      std::uniform_int_distribution<size_t> normal_length_dis(min_length, max_length);
      std::uniform_int_distribution<size_t> large_length_dis(min_length + 100, max_length + 100);

      random_raw_string_data.reserve(size + 1);
      for (size_t i = 0; i < size / 4; ++i) {
        size_t length = large_length_dis(rand_gen);
        for (size_t j = 0; j < length; ++j)
          random_raw_string_data.emplace_back(small_char_dis(rand_gen));
        random_raw_string_data.emplace_back(Char(0));
      }
      for (size_t i = size / 4; i < size; ++i) {
        size_t length = normal_length_dis(rand_gen);
        for (size_t j = 0; j < length; ++j)
          random_raw_string_data.emplace_back(char_dis(rand_gen));
        random_raw_string_data.emplace_back(Char(0));
      }
      this->update(std::move(random_raw_string_data));
    }
  };
}
namespace dsss {

class random_string_set : public string_set {

public:
  random_string_set(const size_t size,
    const size_t alphabet_size =
      std::numeric_limits<dsss::char_type>::max());

  random_string_set(const size_t size, const size_t min_length,
    const size_t max_length);

  random_string_set(random_string_set&&) = default;
  random_string_set& operator =(random_string_set&&) = default; 
  random_string_set(const random_string_set&) = delete;
  random_string_set& operator =(const random_string_set&) = delete;

}; // class random_string_set

template <typename IndexType>
class random_indexed_string_set : public indexed_string_set<IndexType> {

public:
  random_indexed_string_set(const size_t size,
    const size_t alphabet_size =
      std::numeric_limits<dsss::char_type>::max()) {

    std::random_device rand_seed;
    std::mt19937 rand_gen(rand_seed());
    std::uniform_int_distribution<dsss::char_type> char_dis(1,
      std::min<size_t>(
        std::numeric_limits<dsss::char_type>::max(), alphabet_size + 1));

    this->strings_raw_data_.reserve(size + 1);
    for (size_t i = 0; i < size; ++i) {
      this->strings_raw_data_.emplace_back(char_dis(rand_gen));
    }
    this->strings_raw_data_.emplace_back(dsss::char_type(0));
    this->idxd_strings_.emplace_back(indexed_string<IndexType> { 0,
      this->strings_raw_data_.data() });
  }

  random_indexed_string_set(const size_t size,
    const size_t min_length, const size_t max_length) {

    std::random_device rand_seed;
    std::mt19937 rand_gen(rand_seed());
    std::uniform_int_distribution<size_t> length_dis(
      min_length, max_length);
    std::uniform_int_distribution<dsss::char_type> char_dis(1,
      std::numeric_limits<dsss::char_type>::max());

    for (size_t i = 0; i < size; ++i) {
      const size_t string_length = length_dis(rand_gen);
      for (size_t j = 0; j < string_length; ++j) {
        this->strings_raw_data_.emplace_back(char_dis(rand_gen));
      }
      this->strings_raw_data_.emplace_back(dsss::char_type(0));
    }

    IndexType index = { 0 };
    this->idxd_strings_.emplace_back(
      indexed_string<IndexType> { index++, this->strings_raw_data_.data() });
    for (size_t i = 0; i < this->strings_raw_data_.size(); ++i) {
      while (this->strings_raw_data_[i++] != 0) { }
      this->idxd_strings_.emplace_back(indexed_string<IndexType> {
        index++, this->strings_raw_data_.data() + i });
    }
   this->idxd_strings_.pop_back();
  }

  random_indexed_string_set(random_indexed_string_set&&) = default;
  random_indexed_string_set& operator =(random_indexed_string_set&&) = default;

  random_indexed_string_set(const random_indexed_string_set&) = delete;
  random_indexed_string_set& operator =(
    const random_indexed_string_set&) = delete;

}; // class random_indexed_string_set

} // namespace dsss

/******************************************************************************/
