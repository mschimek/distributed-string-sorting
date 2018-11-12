/*******************************************************************************
 * src/tools/stringptr.hpp
 *
 * StringPtr, StringPtrOut, and NoLcpCalc specializations. Encapsulate string
 * and shadow array pointers.
 *
 * Additionally: LcpStringPtr encapsulates string and lcp arrays, which may be
 * interleaved or separate.
 *
 * StringLcpPtr               -> (string,lcp,size)
 * StringLcpCachePtr          -> (string,lcp,charcache,size)
 *
 * StringShadowPtr            -> (string,shadow,size,flip)
 * StringShadowOutPtr         -> (string,shadow,output,size,flip)
 * StringShadowLcpPtr         -> (string,shadow=lcp,size,flip)
 * StringShadowLcpOutPtr      -> (string,shadow=lcp,output,size,flip)
 * StringShadowLcpCacheOutPtr -> (string,shadow=lcp,charcache,output,size,flip)
 *
 *******************************************************************************
 * Copyright (C) 2013-2014 Timo Bingmann <tb@panthema.net>
 * Copyright (C) 2013-2014 Andreas Eberle <email@andreas-eberle.com>
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program.  If not, see <http://www.gnu.org/licenses/>.
 ******************************************************************************/

#ifndef PSS_SRC_TOOLS_STRINGPTR_HEADER
#define PSS_SRC_TOOLS_STRINGPTR_HEADER

#include <cassert>
#include <stdint.h>
#include <numa.h>
#include "stringset.hpp"
#include <iostream>

namespace dss_schimek {

  typedef uintptr_t lcp_t;

  /******************************************************************************/

  template <typename CharType>
    class LcpStringContainer
    {
      public:
        LcpStringContainer() = default;
        LcpStringContainer(std::vector<CharType>&& raw_strings) : 
          raw_strings(std::move(raw_strings))
        {
          update(); 
        }
        explicit LcpStringContainer(std::vector<CharType>&& raw_strings, std::vector<size_t>&& lcp) : 
          raw_strings(std::move(raw_strings)) {

              update();
              lcp_values = std::move(lcp);
              std::cout << "lcp_values: size: " << lcp_values.size() << " size: " << strings.size() << std::endl;
              }
        void print() const
        {
          for (size_t i = 0; i < raw_strings.size(); ++i)
          {
            LOG1 << static_cast<uint>(raw_strings[i]) << " " << raw_strings[i];
          }
          //auto charIterator = raw_strings.begin();
          //while(charIterator < raw_strings.end())
          //{
          //    LOG1 << "[" << i++ << "] = " << &(*charIterator)
          //     << " lcp: " << lcp_values[i];
          //    while(charIterator < raw_strings.end() && *charIterator != 0) ++charIterator;
          //    if(charIterator < raw_strings.end()) 
          //      ++charIterator;
          //}
        }
        void make_raw_strings_contiguous()
        {

          std::vector<CharType> tmp(raw_strings.size(), 0);
          size_t offset = 0;
          for (size_t i = 0; i < get_size(); ++i)
          {
            
            auto& charIterator = strings[i];
            bool eos = false;
            bool is_set = false;
            while(!eos) {
              tmp.at(offset) = *(charIterator); 
              eos = *charIterator == 0 ? true : false;
              if(!is_set)
              {
                 is_set = true;
                strings[i] = tmp.data() + offset;
              }
              ++charIterator;
              ++offset;
            } 
            strings[i] = tmp.data() + offset;
          }
          raw_strings = std::move(tmp);
        }
        size_t get_size() const { return size; }
        void update()
        {
          strings.clear();
          strings.reserve(raw_strings.size() / approx_string_length);
          for (auto char_iterator = raw_strings.begin(); char_iterator < raw_strings.end(); ++char_iterator)
          {
            strings.push_back(raw_strings.data() + (char_iterator - raw_strings.begin()));
            while(*char_iterator != 0)
            {
              ++char_iterator;
              assert(char_iterator < raw-strings.end()); //remove for performance
            }
            size = strings.size();
            lcp_values.resize(size, 0);
         }
        }

        void update(std::vector<CharType>&& raw_data)
        {
          raw_strings = std::move(raw_data); 
          update();
        }

        static constexpr size_t approx_string_length = 10;
        size_t size = 0;
        std::vector<CharType> raw_strings;
        std::vector<CharType*> strings;
        std::vector<size_t> lcp_values;
    };

  using LcpStringContainerUChar = LcpStringContainer<unsigned char>;

  template <typename CharType>
    class LcpStringPtr 
    : public dss_schimek::StringSetBase<LcpStringPtr<CharType>,
    dss_schimek::GenericCharStringSetTraits<CharType> >
  {
    public:
      typedef dss_schimek::GenericCharStringSetTraits<CharType> Traits;

      typedef typename Traits::Char Char;
      typedef typename Traits::String String;
      typedef typename Traits::Iterator Iterator;
      typedef typename Traits::CharIterator CharIterator;
      typedef typename std::tuple<Iterator, size_t*, size_t> Container;
      typedef size_t LcpValue;
      typedef size_t* LcpIterator;
      
      LcpStringPtr() : begin_(nullptr), end_(nullptr), lcp_begin(nullptr) {}

      LcpStringPtr(Iterator begin, Iterator end, LcpIterator lcp_begin) 
        : begin_(begin), end_(end), lcp_begin(lcp_begin)
      {}

      LcpStringPtr(const LcpStringContainer<CharType>& lcp_string_container)
      {
        begin_ = const_cast<Iterator>(lcp_string_container.strings.data());
        end_ =  const_cast<Iterator>(begin_ + lcp_string_container.size);
        lcp_begin = const_cast<size_t*>(lcp_string_container.lcp_values.data());
      }

      explicit LcpStringPtr(const Container& c)
        : begin_(std::get<0>(c)), end_(std::get<0>(c) + std::get<2>(c)), lcp_begin(std::get<1>(c)) {}

      size_t size() const { return end_ - begin_; }

      Iterator begin() const { return begin_; }

      Iterator end() const { return end_; }

      String& operator [] (Iterator i) const { return *i; }

      static Container allocate(size_t n)
      {
        return std::make_tuple(new String[n], new size_t[n], n);
      }
      static void deallocate(Container& c)
      {
        delete[] std::get<0>(c); delete[] std::get<1>(c);
        std::get<0>(c) = NULL; std::get<1>(c) = NULL;
      }
      CharIterator get_chars(const String& s, size_t depth) const
      { return s + depth; }

      bool is_end(const String&, const CharIterator& i) const
      { return (*i == 0); }

      std::string get_string(const String& s, size_t depth = 0) const
      { return std::string(reinterpret_cast<const char*>(s) + depth); }

      LcpStringPtr sub(Iterator begin, Iterator end) const
      {
        return LcpStringPtr(begin, end, lcp_begin + (begin - begin_));
      }

      size_t get_lcp(size_t pos) const
      {
        return *(lcp_begin + pos);
      }

      void set_lcp(size_t pos, size_t value) const
      {
        assert(pos >= (end_ -begin_));
        *(lcp_begin + pos) = value;
      }
      void print() const
      {
        size_t i = 0;
        for (Iterator pi = begin(); pi != end(); ++pi)
        {
          LOG1 << "[" << i++ << "] = " << *pi
            << " = " << get_string(*pi, 0) << " lcp: " << get_lcp(pi - begin_);
        }
      }

      std::vector<CharType> get_raw_strings()
      {
        std::vector<CharType> rawStrings;
        rawStrings.reserve(end_ - begin_);
        for(auto curStr = begin_; curStr < end_; ++curStr) {
          CharIterator charIterator = get_chars(*curStr, 0);
          while(!is_end(*curStr, charIterator))  
            rawStrings.push_back(*(charIterator++));
          rawStrings.push_back(0); // end of str
        }
        return rawStrings;
      }

      std::vector<size_t> get_lcp_vector()
      {
        return std::vector<size_t>(lcp_begin, lcp_begin + size());
      }

      size_t get_string_length(const String& str) const
      {
        CharIterator charIt = get_chars(str, 0);
        size_t length = 0;
        while(!is_end(str, charIt)) ++length;
        return length;
      }

    protected:
      Iterator begin_, end_;
      LcpIterator lcp_begin;
  };
  
  using LcpStringPtrUChar = LcpStringPtr<unsigned char>;
}
#endif
/******************************************************************************/
