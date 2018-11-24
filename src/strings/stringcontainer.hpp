#pragma once

#include <cassert>
#include <stdint.h>
#include <memory>
#include <numa.h>
#include "stringset.hpp"
#include <iostream>
#include <tlx/logger.hpp>
#include <algorithm>
#include "strings/stringtools.hpp"

namespace dss_schimek {

  template <typename CharType>
    class StringLcpContainer_
    {
      public:
        using String = CharType*;
        StringLcpContainer_() = default;
        StringLcpContainer_(std::vector<CharType>&& raw_strings) : 
          raw_strings_(std::move(raw_strings))
      {
        update_strings(); 
        lcps_.resize(size(), 0);
      }

        explicit StringLcpContainer_(std::vector<CharType>&& raw_strings, std::vector<size_t>&& lcp) : 
          raw_strings_(std::move(raw_strings)) {

            update_strings();
            lcps_ = std::move(lcp);
          } 

        String operator[] (size_t i) { return strings_[i]; }
        String front() { return strings_.front(); } 
        String back() { return strings_.back(); }
        String* strings() { return strings_.data();}
        size_t size() const { return strings_.size(); }
        size_t char_size() const { return raw_strings_.size(); }
        std::vector<size_t>& lcps() { return lcps_; }
        size_t* lcp_array() { return lcps_.data(); }
        std::vector<CharType>& raw_strings() { return raw_strings_; }

        void set(std::vector<CharType>&& raw_strings) { raw_strings_ = std::move(raw_strings); }
        void set(std::vector<size_t>&& lcps) { lcps_ = std::move(lcps); }

        void update(std::vector<CharType>&& raw_strings)
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
        std::vector<CharType> raw_strings_;
        std::vector<String> strings_;
        std::vector<size_t> lcps_;

        void update_strings()
        {
          strings_.clear();
          size_t approx_string_size = raw_strings_.size() / approx_string_length;
          strings_.reserve(approx_string_size);
          for (size_t i = 0; i < raw_strings_.size(); ++i)
          {
            strings_.emplace_back(raw_strings_.data() + i);
            while(raw_strings_[i] != 0) ++i;
          }
        }
    };


  template <typename CharType>
    class StringLcpContainer
    {
      public:
        using String = CharType*;
        StringLcpContainer() : raw_strings_(std::make_unique<std::vector<CharType>>()),
          strings_(),
          lcps_() {}
        StringLcpContainer(std::vector<CharType>&& raw_strings) : 
          raw_strings_(std::make_unique<std::vector<CharType>>(std::move(raw_strings)))
      {
        update_strings(); 
        lcps_.resize(size(), 0);
      }

        explicit StringLcpContainer(std::vector<CharType>&& raw_strings, std::vector<size_t>&& lcp) : 
          raw_strings_(std::make_unique<std::vector<CharType>>(std::move(raw_strings))) {

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
        size_t* lcp_array() { return lcps_.data(); }
        std::vector<CharType>& raw_strings() { return *raw_strings_; }
        

        void set(std::vector<CharType>&& raw_strings) { *raw_strings_ = std::move(raw_strings); }
        void set(std::vector<String>&& strings) { strings_ = std::move(strings);  }
        void set(std::vector<size_t>&& lcps) { lcps_ = std::move(lcps); }

        void update(std::vector<CharType>&& raw_strings)
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
        std::unique_ptr<std::vector<CharType>> raw_strings_;
        std::vector<String> strings_;
        std::vector<size_t> lcps_;

        void update_strings()
        {
          strings_.clear();
          size_t approx_string_size = raw_strings_->size() / approx_string_length;
          strings_.reserve(approx_string_size);
          for (size_t i = 0; i < raw_strings_->size(); ++i)
          {
            strings_.emplace_back(raw_strings_->data() + i);
            while(raw_strings_->operator[](i) != 0) ++i;
          }
        }
    };

  using StringLcpContainerUChar = StringLcpContainer<unsigned char>;
}
