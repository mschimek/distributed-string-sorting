/******************************************************************************
 * tlx/sort/strings/insertion_sort.hpp
 *
 * Base insertion string sort. This is an internal implementation header, see
 * tlx/sort/strings.hpp for public front-end functions.
 *
 * Part of tlx - http://panthema.net/tlx
 *
 * Copyright (C) 2015-2018 Timo Bingmann <tb@panthema.net>
 *
 * All rights reserved. Published under the Boost Software License, Version 1.0
 ******************************************************************************/

#ifndef TLX_SORT_STRINGS_INSERTION_SORT_HEADER__
#define TLX_SORT_STRINGS_INSERTION_SORT_HEADER__

#include <iostream>
#include <tlx/define/likely.hpp>
#include "../strings/stringset.hpp"

  namespace mysorter {

    /******************************************************************************/

    //! Generic insertion sort for abstract string sets. This method only requires
    //! O(1) additional memory for sorting n strings, but runs in time O(nD).
    template <typename StringSet>
      static inline void insertion_sort(
          StringSet& ss, size_t depth, size_t /* memory */) {
        typedef typename StringSet::Iterator Iterator;
        typedef typename StringSet::String String;
        typedef typename StringSet::CharIterator CharIterator;

        // this stores the begin iterator and size n, making the loops faster
        const Iterator begin = ss.begin();
        Iterator j;
        size_t n = ss.size();
        for (Iterator i = begin + 1; TLX_UNLIKELY(--n != 0); ++i)
        {
          String tmp = std::move(ss[i]);
          j = i;

          while (TLX_LIKELY(j != begin))
          {
            CharIterator s = ss.get_chars(ss[j - 1], depth);
            CharIterator t = ss.get_chars(tmp, depth);

            while (TLX_LIKELY(ss.is_equal(ss[j - 1], s, tmp, t)))
              ++s, ++t;

            if (TLX_UNLIKELY(ss.is_leq(ss[j - 1], s, tmp, t))) {
              break;
            }

            ss[j] = std::move(ss[j - 1]);
            --j;
          }

          ss[j] = std::move(tmp);
        }
        n = ss.size();
        size_t index = 1;
        ss.set_lcp(0, 0);
        for (Iterator i = begin + 1; TLX_UNLIKELY(--n != 0); ++i)
        {
          CharIterator s = ss.get_chars(ss[i], 0);
          CharIterator t = ss.get_chars(ss[i - 1], 0);
          size_t value = 0;
          while ( ss.is_equal(ss[i], s, ss[i - 1], t))
            ++s, ++t, ++value;

          std::cout << "n: " << n << " value:  " << value << std::endl;
          ss.set_lcp(index++, value);
        }
      }

    template <typename StringSet, typename LcpIterator>
      static inline void insertion_sort_lcp(
          const StringSet& ss, size_t depth, size_t /* memory */, LcpIterator lcp) {
        typedef typename StringSet::Iterator Iterator;
        typedef typename StringSet::String String;
        typedef typename StringSet::CharIterator CharIterator;

        // this stores the begin iterator and size n, making the loops faster
        const Iterator begin = ss.begin();
        size_t n = ss.size();

        if (n <= 1) return;

        for (size_t j = 0; j < n - 1; ++j)
        {
          // insert strings[j] into sorted strings[0..j-1]

          String new_str = std::move(ss[begin + j]);
          size_t new_lcp = depth; // start with LCP depth

          size_t i = j;
          while (i > 0)
          {
            size_t prev_lcp = new_lcp;

            String cur_str = std::move(ss[begin + i - 1]);
            size_t cur_lcp = lcp[i];

            if (cur_lcp < new_lcp)
            {
              // CASE 1: lcp goes down -> insert string

              // move comparison string back
              ss[begin + i - 1] = std::move(cur_str);
              break;
            }
            else if (cur_lcp == new_lcp)
            {
              // CASE 2: compare more characters

              CharIterator c1 = ss.get_chars(new_str, new_lcp);
              CharIterator c2 = ss.get_chars(cur_str, new_lcp);

              while (ss.is_equal(new_str, c1, cur_str, c2))
                ++c1, ++c2, ++new_lcp;

              // if (new_str >= curr_str) -> insert string
              if (ss.is_leq(cur_str, c2, new_str, c1))
              {
                // update lcp of prev (smaller string) with inserted string
                lcp[i] = new_lcp;
                // lcp of inserted string with next string
                new_lcp = prev_lcp;

                // move comparison string back
                ss[begin + i - 1] = std::move(cur_str);
                break;
              }
            }
            // else (cur_lcp > new_lcp), CASE 3: nothing to do

            ss[begin + i] = std::move(cur_str);

            lcp[i + 1] = cur_lcp;

            --i;
          }

          ss[begin + i] = std::move(new_str);
          lcp[i + 1] = new_lcp;
        }

        // last loop specialized with checks for out-of-bound access to lcp.
        {
          size_t j = n - 1;

          // insert ssings[j] into sorted ssings[0..j-1]

          String new_ss = std::move(ss[begin + j]);
          size_t new_lcp = depth; // start with LCP depth

          size_t i = j;
          while (i > 0)
          {
            size_t prev_lcp = new_lcp;

            String cur_ss = std::move(ss[begin + i - 1]);
            size_t cur_lcp = lcp[i];

            if (cur_lcp < new_lcp)
            {
              // CASE 1: lcp goes down -> insert ssing

              // move comparison ssing back
              ss[begin + i - 1] = std::move(cur_ss);
              break;
            }
            else if (cur_lcp == new_lcp)
            {
              // CASE 2: compare more characters

              CharIterator c1 = ss.get_chars(new_ss, new_lcp);
              CharIterator c2 = ss.get_chars(cur_ss, new_lcp);

              while (ss.is_equal(new_ss, c1, cur_ss, c2))
                ++c1, ++c2, ++new_lcp;

              // if (new_ss >= curr_ss) -> insert ssing
              if (ss.is_leq(cur_ss, c2, new_ss, c1))
              {
                // update lcp of prev (smaller ssing) with inserted ssing
                lcp[i] = new_lcp;
                // lcp of inserted ssing with next ssing
                new_lcp = prev_lcp;

                // move comparison ssing back
                ss[begin + i - 1] = std::move(cur_ss);
                break;
              }
            }
            // else (cur_lcp > new_lcp), CASE 3: nothing to do

            ss[begin + i] = std::move(cur_ss);

            if (i + 1 < n) // check out-of-bounds copy
              lcp[i + 1] = cur_lcp;

            --i;
          }

          ss[begin + i] = std::move(new_ss);

          if (i + 1 < n) { // check out-of-bounds save
            lcp[i + 1] = new_lcp;
          }
        }
        /*std::cout << "final inside insertion sort" << std::endl;
          for (auto it = lcp; it < lcp + ss.size(); ++it)
          {
          std::cout << *it << " ";
          }
          std::cout << std::endl;*/

      }

    /******************************************************************************/

  } // namespace sort_strings_detail

#endif // !TLX_SORT_STRINGS_INSERTION_SORT_HEADER

