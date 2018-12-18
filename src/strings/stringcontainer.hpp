#pragma once

#include <cassert>
#include <stdint.h>
#include <memory>
#include <numa.h>
#include "strings/stringset.hpp"
#include <iostream>
#include <tlx/logger.hpp>
#include <algorithm>
#include "strings/stringptr.hpp"
#include "strings/stringtools.hpp"

namespace dss_schimek {

  template <typename StringSet>
    class InitPolicy
    {
      using Char = typename StringSet::Char;
      using String = typename StringSet::String;

      static constexpr size_t approx_string_length = 10;
      public: 
      std::vector<String> init_strings(std::vector<Char>& raw_strings) {
        std::vector<String> strings;
        size_t approx_string_size = raw_strings.size() / approx_string_length;
        strings.reserve(approx_string_size);
        for (size_t i = 0; i < raw_strings.size(); ++i)
        {
          strings.emplace_back(raw_strings.data() + i);
          while(raw_strings[i] != 0) ++i;
        }
        return strings;
      }
    };

  template <typename CharType>
    class InitPolicy<GenericCharLengthStringSet<CharType>> 
    {
      using StringSet = GenericCharLengthStringSet<CharType>;
      using Char = typename StringSet::Char;
      using String = typename StringSet::String;

      static constexpr size_t approx_string_length = 10;
      public: 
      std::vector<String> init_strings(std::vector<Char>& raw_strings) {
        std::vector<String> strings;
        size_t approx_string_size = raw_strings.size() / approx_string_length;
        strings.reserve(approx_string_size);
        for (size_t i = 0; i < raw_strings.size(); ++i)
        {
          strings.emplace_back(raw_strings.data() + i, i);
          while(raw_strings[i] != 0) ++i;
          strings.back().length = i - strings.back().length; 
        }
        return strings;
      }
    };

  template <typename StringSet_>
    class StringLcpContainer : private InitPolicy<StringSet_>
  {
    public:
      using StringSet = StringSet_;
      using Char = typename StringSet::Char; 
      using CharIterator = typename StringSet::CharIterator;
      using String = typename StringSet::String;

      StringLcpContainer() : raw_strings_(std::make_unique<std::vector<Char>>()),
      strings_(),
      lcps_() {}

      StringLcpContainer(std::vector<Char>&& raw_strings) : 
        raw_strings_(std::make_unique<std::vector<Char>>(std::move(raw_strings))) {
          update_strings(); 
          lcps_.resize(size(), 0);
        }

      explicit StringLcpContainer(std::vector<Char>&& raw_strings, std::vector<size_t>&& lcp) : 
        raw_strings_(std::make_unique<std::vector<Char>>(std::move(raw_strings))) {

          update_strings();
          lcps_ = std::move(lcp);
        } 

      String operator[] (size_t i) { return strings_[i]; }
      String front() { return strings_.front(); } 
      String back() { return strings_.back(); }
      String* strings() { return strings_.data();}
      size_t size() const { return strings_.size(); }
      size_t char_size() const { return raw_strings_->size(); }
      std::vector<size_t>& lcps() { return lcps_; }
      const std::vector<size_t>& lcps() const { return lcps_; }
      size_t* lcp_array() { return lcps_.data(); }
      std::vector<Char>& raw_strings() { return *raw_strings_; }
      const std::vector<Char>& raw_strings() const { return *raw_strings_; }
      
      

      StringSet make_string_set() {
        return StringSet(strings(), strings() + size());
      }

      dss_schimek::StringPtr<StringSet> make_string_ptr() {
        return dss_schimek::StringPtr(make_string_set());
      }

      dss_schimek::StringLcpPtr<StringSet> make_string_lcp_ptr() {
        return dss_schimek::StringLcpPtr(make_string_set(), lcp_array());
      }


      void set(std::vector<Char>&& raw_strings) { *raw_strings_ = std::move(raw_strings); }
      void set(std::vector<String>&& strings) { strings_ = std::move(strings);  }
      void set(std::vector<size_t>&& lcps) { lcps_ = std::move(lcps); }
      
      bool operator==(const StringLcpContainer<StringSet_>& other) {
        return (raw_strings() == other.raw_strings()) && (lcps() == other.lcps());
      }

      

      void update(std::vector<Char>&& raw_strings)
      {
        set(std::move(raw_strings));
        update_strings(); 
        if (lcps_.size() != size())
          lcps_.resize(size(), 0);
      }

      bool is_consistent()
      {
        if (lcps_.size() != strings_.size()) 
        {
          LOG1 << "lcps.size() = " << lcps_.size() << " != " << strings_.size() << " = strings.size()";
          return false;
        }

        return std::all_of(strings_.begin(), strings_.end(), [this](const String str) -> bool {
            if (str < raw_strings_.data() || str > raw_strings_.data() + raw_strings_.size())
            return false;
            if (str == raw_strings_.data())
            return true;
            return *(str - 1) == 0;
            });
      }

    protected:
      static constexpr size_t approx_string_length = 10;
      std::unique_ptr<std::vector<Char>> raw_strings_;
      std::vector<String> strings_;
      std::vector<size_t> lcps_;

      void update_strings()
      {
        strings_ = InitPolicy<StringSet>::init_strings(*raw_strings_);
      }
  };

  using StringLcpContainerUChar = StringLcpContainer<UCharStringSet>;
}
