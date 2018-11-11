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

namespace dss_schimek {

  typedef uintptr_t lcp_t;

  /******************************************************************************/


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

      LcpStringPtr(Iterator begin, Iterator end, LcpIterator lcp_begin) 
        : begin_(begin), end_(end), lcp_begin(lcp_begin)
      {}
      
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
        if (pos >= (end_ - begin_))
          std::cout << " !!!!!!!! pos too big !!!!!!!" << std::endl;
        *(lcp_begin + pos) = value;
      }
      void print() const
      {
        size_t i = 0;
        for (Iterator pi = begin(); pi != end(); ++pi)
        {
          LOG1 << "[" << i++ << "] = " << *pi
            << " = " << get_string(*pi, 0);
        }
      }

    protected:
      Iterator begin_, end_;
      LcpIterator lcp_begin;
  };
}
#endif
/******************************************************************************/
