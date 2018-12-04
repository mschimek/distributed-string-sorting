/*******************************************************************************
 * tests/sort_strings_test.cpp
 *
 * String sorting test program
 *
 * Part of tlx - http://panthema.net/tlx
 *
 * Copyright (C) 2015-2018 Timo Bingmann <tb@panthema.net>
 *
 * All rights reserved. Published under the Boost Software License, Version 1.0
 ******************************************************************************/
#include <iostream>
#include "sorter/local/strings/insertion_sort_unified.hpp"
#include "sorter/local/strings/multikey_quicksort_unified.hpp"
#include "sorter/local/strings/radix_sort_unified.hpp"
#include "strings/stringset.hpp" 
#include "strings/stringptr.hpp"

#include <tlx/logger.hpp>

#include <chrono>
#include <random>

#if TLX_MORE_TESTS
static const bool tlx_more_tests = true;
#else
static const bool tlx_more_tests = false;
#endif

static const size_t seed = 1234567;

template <typename Set>
using StringSorter_ = void (*)(const dss_schimek::StringLcpPtr<Set>&, size_t, size_t );
template <typename Set>
using StringSorter = void (*)(const Set&, size_t, size_t);
template <typename Set, typename LcpIterator>
using StringLCPSorter = void (*)(const Set&, size_t, size_t, LcpIterator);

template<typename StringSet, typename LCPContainer>
bool check_lcps(const StringSet& sorted_strings, const LCPContainer& lcp_to_be_tested)
{
  typedef typename StringSet::Iterator Iterator;
  typedef typename StringSet::String String;
  typedef typename StringSet::CharIterator CharIterator;

  const Iterator begin = sorted_strings.begin();
  size_t n = sorted_strings.size();

  if(sorted_strings.size() > 0 && lcp_to_be_tested[0] != 0)
  {
    std::cout << "i: " << 0 << " lcp_value: " << lcp_to_be_tested[0] << std::endl;
    return false;
  }
  
  for(Iterator i = begin + 1; --n > 0; ++i)
  {
    CharIterator s = sorted_strings.get_chars(sorted_strings[i - 1], 0);
    CharIterator t = sorted_strings.get_chars(sorted_strings[i], 0);
    size_t actual_lcp = 0;
    while (TLX_LIKELY(sorted_strings.is_equal(sorted_strings[i - 1], s, sorted_strings[i], t)))
      ++s, ++t, ++actual_lcp;

    //std::cout << "i: " << i - begin << " actual_lcp: " << actual_lcp << " computed lcp " << lcp_to_be_tested[i - begin] << std::endl;

      
    if (actual_lcp != lcp_to_be_tested[i - begin])
     return false;
  }
  return true;
}

template <typename Random, typename Iterator>
void fill_random(Random& rng, const std::string& letters,
                 Iterator begin, Iterator end) {
    for (Iterator i = begin; i != end; ++i)
        *i = letters[(rng() / 100) % letters.size()];
}

//! Returns number of seconds since the epoch, high resolution.
static inline double timestamp() {
    return static_cast<double>(
        std::chrono::duration_cast<std::chrono::microseconds>(
            std::chrono::steady_clock::now().time_since_epoch()).count()) / 1e6;
}

template <typename StringSet, StringSorter_<StringSet> sorter>
void TestUCharString(const char* name,
                     const size_t num_strings, const size_t num_chars,
                     const std::string& letters) {

    std::default_random_engine rng(seed);

    LOG1 << "Running " << name << " on " << num_strings
         << " uint8_t* strings";

    // array of string pointers
    uint8_t** cstrings = new uint8_t*[num_strings];

    // generate random strings of length num_chars
    for (size_t i = 0; i < num_strings; ++i)
    {
        size_t slen = num_chars + (rng() >> 8) % (num_chars / 4);

        cstrings[i] = new uint8_t[slen + 1];
        fill_random(rng, letters, cstrings[i], cstrings[i] + slen);
        cstrings[i][slen] = 0;
    }

    // run sorting algorithm
    double ts1 = timestamp();

    dss_schimek::UCharStringSet ss(cstrings, cstrings + num_strings);
    dss_schimek::LcpType* lcps = new dss_schimek::LcpType[ss.size()];
    if (ss.size() > 0)
      lcps[0] = 0;
    dss_schimek::StringLcpPtr strptr(ss, lcps);
    sorter(strptr, /* depth */ 0, /* memory */ 0);
    if (0) ss.print();

    double ts2 = timestamp();
    LOG1 << "sorting took " << ts2 - ts1 << " seconds";

    // check result
    if (!ss.check_order()) {
        LOG1 << "Result is not sorted!";
        abort();
    }

    if (!check_lcps(ss, lcps))
    {
      LOG1 << "lcp value not correct!";
      //for (size_t i = 0; i < ss.size(); ++i)
      //  std::cout << i << "  " << lcp[i] << std::endl; 
      std::cout << std::endl;
      abort();
    }


    // free memory.
    for (size_t i = 0; i < num_strings; ++i)
        delete[] cstrings[i];

    delete[] cstrings;
}

static const char* letters_alnum
    = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"
      "\xE0\xE1\xE2\xE3\xE4\xE5\xE6\xE7\xE8\xE9\xEA\xEB\xEC\xED\xEE\xEF";

#define run_lcp_tests(func)                                          \
TestUCharString<dss_schimek::UCharStringSet, func>(         \
        #func, num_strings, 16, letters_alnum);                      \

void test_all_lcp(const size_t num_strings) {
  if (num_strings <= 1024) {
   run_lcp_tests(dss_schimek::insertion_sort);
  }
  letters_alnum = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"
      "\xE0\xE1\xE2\xE3\xE4\xE5\xE6\xE7\xE8\xE9\xEA\xEB\xEC\xED\xEE\xEF";

  std::cout << "original test string: " << letters_alnum << std::endl;
   run_lcp_tests(dss_schimek::multikey_quicksort);
   run_lcp_tests(dss_schimek::radixsort_CE0);
   run_lcp_tests(dss_schimek::radixsort_CE2);
   run_lcp_tests(dss_schimek::radixsort_CE3);
  letters_alnum = "a";
  std::cout << "test string: " << letters_alnum << std::endl;
   run_lcp_tests(dss_schimek::multikey_quicksort);
   run_lcp_tests(dss_schimek::radixsort_CE0);
   run_lcp_tests(dss_schimek::radixsort_CE2);
   run_lcp_tests(dss_schimek::radixsort_CE3);

  letters_alnum = "abc";
  std::cout << "test string: " << letters_alnum << std::endl;
   run_lcp_tests(dss_schimek::multikey_quicksort);
   run_lcp_tests(dss_schimek::radixsort_CE0);
   run_lcp_tests(dss_schimek::radixsort_CE2);
   run_lcp_tests(dss_schimek::radixsort_CE3);
 letters_alnum = "acbef";
  std::cout << "test string: " << letters_alnum << std::endl;
   run_lcp_tests(dss_schimek::multikey_quicksort);
   run_lcp_tests(dss_schimek::radixsort_CE0);
   run_lcp_tests(dss_schimek::radixsort_CE2);
   run_lcp_tests(dss_schimek::radixsort_CE3);

}
int main() {
    //test_all(16);
    //test_all(256);
    //test_all(65550);
    //if (tlx_more_tests) {
    //    test_all(1024 * 1024);
    //    // test_all(16 * 1024 * 1024);
    //}
    //test_all_lcp(16);
    //test_all_lcp(65550);
    //test_all(65550);
//    TestUCharString<dss_schimek::UCharStringSet, dss_schimek::insertion_sort>("hallo", 56, 16, letters_alnum);

// array of string pointers
  
   test_all_lcp(16);
    test_all_lcp(256);
    test_all_lcp(65550);
    test_all_lcp(655550);
    test_all_lcp(6555550);
    
    //std::vector<std::string> strings = {"aaabbabbbbbabcbbaac", "aaabbbabbcaacbbcb", "aaabbccaccaaccbc", "aaacabaabcbbaacba", "aaacabcbacaccccaa", "aaacccbacabcccbca", "aababcabbaacacbcab"} ;
    //StdStringSet ss(strings.data(), strings.data() + strings.size());
    //std::vector<uint> lcps(strings.size());
    //insertion_sort_lcp(ss, 2, 0, lcps.begin());
    //std::cout << "final lcp:" << std::endl;
    //for(const auto&i : lcps)
    //  std::cout << i << " ";
    //std::cout << std::endl;
    return 0;
}

/******************************************************************************/
