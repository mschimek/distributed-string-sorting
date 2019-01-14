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
#include <cmath>
#include <algorithm>

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

    static std::string getName() {
      return "PrefixStringGenerator";
    }
  };
  
  template <typename StringSet>
  class DNRationGeneratorStringLcpContainer : public StringLcpContainer<StringSet>
  {
    using Char = typename StringSet::Char;
    public:
    std::vector<unsigned char> nextChar(const std::vector<unsigned char>& lastChar, const size_t min, const size_t max) {
      std::vector<unsigned char> nextChar(lastChar.size(), min);
      int64_t i = lastChar.size() - 1;
      for (; i >= 0; --i) {
        if (lastChar[i] < max) {
          nextChar[i] = lastChar[i] + 1;
          std::copy(lastChar.begin(), lastChar.begin() + i, nextChar.begin());
          break;
        }
      }
      return nextChar;
    }

    std::pair<std::vector<unsigned char>, size_t> rawStrings(size_t numStrings, size_t desiredStringLength, double dToN) {
      std::vector<unsigned char> rawStrings;
      const size_t minInternChar = 65;
      const size_t maxInternChar = 90;
      const size_t numberInternChars = maxInternChar - minInternChar + 1;
      const size_t charLength = std::ceil(0.5 * std::log(numStrings) / log(numberInternChars));
      std::cout << "charLength: " << charLength << std::endl;
      const size_t commonPrefixLength = std::max(static_cast<int64_t>(desiredStringLength * dToN - 2 * charLength), 0l);
      std::cout << "commonPrefixLength: " << commonPrefixLength << std::endl;
      const size_t paddingLength = std::max(static_cast<int64_t>(desiredStringLength - (commonPrefixLength + 2 * charLength)), 0l);
      std::cout << "paddingLength: " << paddingLength << std::endl;
      const size_t stringLength = commonPrefixLength + 2 * charLength + paddingLength;
      const size_t wrap = std::pow(static_cast<double>(maxInternChar - minInternChar + 1), charLength);
      std::cout << "wrap: " << wrap << std::endl;
      std::vector<unsigned char> curFirstChar(charLength, minInternChar);
      std::vector<unsigned char> curSecondChar(charLength, minInternChar);

      for (size_t i = 0; i < numStrings; ++i) {
        for (size_t j = 0; j < commonPrefixLength; ++j)
          rawStrings.emplace_back(maxInternChar);

        for (size_t j = 0; j < charLength; ++j) 
          rawStrings.emplace_back(curFirstChar[j]);
        for (size_t j = 0; j < charLength; ++j) 
          rawStrings.emplace_back(curSecondChar[j]);

        for (size_t j = 0; j < paddingLength; ++j)
          rawStrings.emplace_back(maxInternChar);
        rawStrings.emplace_back(0);

        if ((i + 1)% wrap == 0)
          curFirstChar = nextChar(curFirstChar, minInternChar, maxInternChar);
        curSecondChar = nextChar(curSecondChar, minInternChar, maxInternChar);
      }
      std::cout << "finished call" << std::endl;
      return make_pair(rawStrings, stringLength);
    }

    public:
    DNRationGeneratorStringLcpContainer(const size_t size,
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

    static std::string getName() {
      return "DNRationGeneratorStringLcpContainer";
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

    static std::string getName() {
      return "RandomStringGenerator";
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
      dsss::mpi::environment env;
      std::mt19937 rand_gen(env.rank());//rand_seed());
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

    static std::string getName() {
      return "SkewedStringGenerator";
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
