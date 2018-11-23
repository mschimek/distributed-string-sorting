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

#include "mpi/environment.hpp"

namespace dss_schimek {
  using namespace dss_schimek;

  template <typename CharType>
  class PrefixNumberStringLcpContainer : public StringLcpContainer<CharType>
  {
    public:
    PrefixNumberStringLcpContainer(const size_t size, const CharType prefix)
    {
      std::vector<CharType> raw_string_data;
      for(size_t i = 1; i <= size; ++i)
      {
        raw_string_data.emplace_back(prefix);
        size_t number_to_print = i;
        while(number_to_print > 0)
        {
          CharType last_digit = number_to_print % 10;
          raw_string_data.emplace_back(48 + last_digit);
          number_to_print /= 10;
        }
        raw_string_data.emplace_back(CharType(0));
      }
      this->update(std::move(raw_string_data));
    }
  };

  template <typename CharType>
  class RandomStringLcpContainer : public StringLcpContainer<CharType>
  {
    public:
    RandomStringLcpContainer(const size_t size,
          const size_t min_length = 10,
          const size_t max_length = 20)
    {
      std::vector<CharType> random_raw_string_data;
      std::random_device rand_seed;
      std::mt19937 rand_gen(dsss::mpi::environment().rank());
      std::uniform_int_distribution<CharType> char_dis(65, 90);
      
      std::uniform_int_distribution<size_t> length_dis(min_length, max_length);
      random_raw_string_data.reserve(size + 1);
      for (size_t i = 0; i < size; ++i) {
        size_t length = length_dis(rand_gen);
        for (size_t j = 0; j < length; ++j)
          random_raw_string_data.emplace_back(char_dis(rand_gen));
        random_raw_string_data.emplace_back(CharType(0));
      }
      std::cout << "hallo" << std::endl;
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
