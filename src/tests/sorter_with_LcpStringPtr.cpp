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
//#include <tlx/sort/strings/insertion_sort.hpp>
//#include <tlx/sort/strings/insertion_sort_lcp.hpp>
//#include <tlx/sort/strings/multikey_quicksort.hpp>
//#include <tlx/sort/strings/multikey_quicksort_lcp.hpp>
//#include <tlx/sort/strings/radix_sort.hpp>
//#include <tlx/sort/strings/radix_sort_lcp.hpp>

#include <tlx/sort/strings/insertion_sort_unified.hpp>
#include <tlx/sort/strings/multikey_quicksort_unified.hpp>
#include <tlx/sort/strings/radix_sort_unified.hpp>

#include <tlx/sort/strings.hpp>

#include <tlx/logger.hpp>

#include <chrono>
#include <random>

#include "../strings/stringtools.hpp"
#include "../strings/stringset.hpp"

#if TLX_MORE_TESTS
static const bool tlx_more_tests = true;
#else
static const bool tlx_more_tests = false;
#endif

using namespace tlx::sort_strings_detail;

static const size_t seed = 1234567;

template <typename Set>
using StringSorter = void (*)(const Set&, size_t, size_t);
template <typename Set>
using StringSorterUnified = void (*)(const Set&, size_t, size_t);
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
    return false;
  
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

template<typename StringSet>
bool check_lcps(const StringSet& sorted_strings)
{
  typedef typename StringSet::Iterator Iterator;
  typedef typename StringSet::String String;
  typedef typename StringSet::CharIterator CharIterator;

  const Iterator begin = sorted_strings.begin();
  size_t n = sorted_strings.size();

  if (n > 0 && sorted_strings.get_lcp(0) != 0)
  {

    std::cout << "i: " << 0 << " computed lcp " << sorted_strings.get_lcp(0) << std::endl;
    return false;
  }
  std::size_t index = 1; 
  bool correct = true;
  for(Iterator i = begin + 1; --n > 0; ++i)
  {
    CharIterator s = sorted_strings.get_chars(sorted_strings[i - 1], 0);
    CharIterator t = sorted_strings.get_chars(sorted_strings[i], 0);
    size_t actual_lcp = 0;
    while (TLX_LIKELY(sorted_strings.is_equal(sorted_strings[i - 1], s, sorted_strings[i], t)))
      ++s, ++t, ++actual_lcp;

    //std::cout << "i: " << i - begin << " actual_lcp: " << actual_lcp << " computed lcp " << sorted_strings.get_lcp(index) << std::endl;

      
    if (actual_lcp != sorted_strings.get_lcp(index))
     return false;
    ++index;
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

template <typename StringSet, StringSorterUnified<StringSet> sorter>
void TestLcpStringPtr(const char* name,
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
    size_t* lcp_values = new size_t[num_strings]; 
    if (num_strings)
      lcp_values[0] = 0;
    // run sorting algorithm
    double ts1 = timestamp();

    dss_schimek::LcpStringPtr<unsigned char> ss(cstrings, cstrings + num_strings, lcp_values);
    sorter(ss, /* depth */ 0, /* memory */ 0);
    double ts2 = timestamp();
    LOG1 << "sorting took " << ts2 - ts1 << " seconds";

    // check result
    if (!ss.check_order()) {
        LOG1 << "Result is not sorted!";
        abort();
    }
    
    if (!check_lcps(ss))
    {
      LOG1 << "lcp value not correct!";
      /*for (size_t i = 0; i < ss.size(); ++i)
        std::cout << i << "  " << ss.get_lcp(i) << std::endl; 
      std::cout << std::endl;*/
      abort();
    }

    // free memory.
    for (size_t i = 0; i < num_strings; ++i)
        delete[] cstrings[i];

    delete[] cstrings;
    delete[] lcp_values;
}

template <typename StringSet, StringSorter<StringSet> sorter>
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

    UCharStringSet ss(cstrings, cstrings + num_strings);
    sorter(ss, /* depth */ 0, /* memory */ 0);
    if (0) ss.print();

    double ts2 = timestamp();
    LOG1 << "sorting took " << ts2 - ts1 << " seconds";

    // check result
    if (!ss.check_order()) {
        LOG1 << "Result is not sorted!";
        abort();
    }
    

    // free memory.
    for (size_t i = 0; i < num_strings; ++i)
        delete[] cstrings[i];

    delete[] cstrings;
}

template <typename StringSet, typename LCPContainer, StringLCPSorter<StringSet, typename LCPContainer::iterator> sorter>
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

    UCharStringSet ss(cstrings, cstrings + num_strings);
    LCPContainer lcp(num_strings);

    sorter(ss, /* depth */ 0, /* memory */ 0, lcp.begin());
    if (0) ss.print();

    double ts2 = timestamp();
    LOG1 << "sorting took " << ts2 - ts1 << " seconds";
    // check result

    if (!ss.check_order()) {
        LOG1 << "Result is not sorted!";
        abort();
    }
    if (!check_lcps(ss, lcp))
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

template <typename StringSet, StringSorter<StringSet> sorter>
void TestVectorStdString(const char* name,
                         const size_t num_strings, const size_t num_chars,
                         const std::string& letters) {

    std::default_random_engine rng(seed);

    LOG1 << "Running " << name << " on " << num_strings
         << " std::vector<std::string> strings";

    // vector of std::string objects
    std::vector<std::string> strings(num_strings);

    // generate random strings of length num_chars
    for (size_t i = 0; i < num_strings; ++i)
    {
        size_t slen = num_chars + (rng() >> 8) % (num_chars / 4);

        strings[i].resize(slen);
        fill_random(rng, letters, strings[i].begin(), strings[i].end());
    }

    // run sorting algorithm
    double ts1 = timestamp();

    StdStringSet ss(strings.data(), strings.data() + strings.size());
    sorter(ss, /* depth */ 0, /* memory */ 0);
    if (0) ss.print();

    double ts2 = timestamp();
    LOG1 << "sorting took " << ts2 - ts1 << " seconds";

    // check result
    if (!ss.check_order()) {
        LOG1 << "Result is not sorted!";
        abort();
    }
}
/***************************** TEST ***************************************/

template <typename StringSet, 
          bool compute_lcp, 
          typename std::conditional<compute_lcp, 
                                    StringLCPSorter<StringSet, typename std::vector<int>::iterator>,  
                                    StringSorter<StringSet>>::type sorter>
void TestVectorStdString_test(const char* name,
                         const size_t num_strings, const size_t num_chars,
                         const std::string& letters) {

    std::default_random_engine rng(seed);

    LOG1 << "Running " << name << " on " << num_strings
         << " std::vector<std::string> strings";

    // vector of std::string objects
    std::vector<std::string> strings(num_strings);

    // generate random strings of length num_chars
    for (size_t i = 0; i < num_strings; ++i)
    {
        size_t slen = num_chars + (rng() >> 8) % (num_chars / 4);

        strings[i].resize(slen);
        fill_random(rng, letters, strings[i].begin(), strings[i].end());
    }

    // run sorting algorithm
    double ts1 = timestamp();

    StdStringSet ss(strings.data(), strings.data() + strings.size());

    if constexpr (compute_lcp) {
      std::vector<int> lcp(ss.size());
      sorter(ss, /* depth */ 0, /* memory */ 0, lcp.begin());
      if (!check_lcps(ss, lcp))
      {
        LOG1 << "lcp value not correct!";
        abort();
      }
    } else
      sorter(ss, 0, 0);
    if (0) ss.print();

    double ts2 = timestamp();
    LOG1 << "sorting took " << ts2 - ts1 << " seconds";

    // check result
    if (!ss.check_order()) {
        LOG1 << "Result is not sorted!";
        abort();
    }
    
}
/*
 *
 *
 *++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

template <typename StringSet, typename LCPContainer, StringLCPSorter<StringSet, typename LCPContainer::iterator> sorter>
void TestVectorStdString(const char* name,
                         const size_t num_strings, const size_t num_chars,
                         const std::string& letters) {

    std::default_random_engine rng(seed);

    LOG1 << "Running " << name << " on " << num_strings
         << " std::vector<std::string> strings";

    // vector of std::string objects
    std::vector<std::string> strings(num_strings);

    // generate random strings of length num_chars
    for (size_t i = 0; i < num_strings; ++i)
    {
        size_t slen = num_chars + (rng() >> 8) % (num_chars / 4);

        strings[i].resize(slen);
        fill_random(rng, letters, strings[i].begin(), strings[i].end());
    }

    // run sorting algorithm
    double ts1 = timestamp();

    StdStringSet ss(strings.data(), strings.data() + strings.size());
    LCPContainer lcp(ss.size());
    sorter(ss, /* depth */ 0, /* memory */ 0, lcp.begin());
    if (0) ss.print();

    double ts2 = timestamp();
    LOG1 << "sorting took " << ts2 - ts1 << " seconds";

    // check result
    if (!ss.check_order()) {
        LOG1 << "Result is not sorted!";
        abort();
    }
    if (!check_lcps(ss, lcp))
    {
      LOG1 << "lcp value not correct!";
      abort();
    }
}

template <typename StringSet, StringSorter<StringSet> sorter>
void TestUPtrStdString(const char* name,
                       const size_t num_strings, const size_t num_chars,
                       const std::string& letters) {

    std::default_random_engine rng(seed);

    LOG1 << "Running " << name << " on " << num_strings
         << " std::vector<std::unique_ptr<std::string>> strings";

    // vector of pointers to std::string objects
    typedef std::unique_ptr<std::string> unique_ptr_string;
    std::vector<unique_ptr_string> strings(num_strings);

    // generate random strings of length num_chars
    for (size_t i = 0; i < num_strings; ++i)
    {
        size_t slen = num_chars + (rng() >> 8) % (num_chars / 4);

        strings[i] = unique_ptr_string(new std::string(slen, 0));
        fill_random(rng, letters, strings[i]->begin(), strings[i]->end());
    }

    // run sorting algorithm
    double ts1 = timestamp();

    UPtrStdStringSet ss(strings.data(), strings.data() + strings.size());
    sorter(ss, /* depth */ 0, /* memory */ 0);
    if (0) ss.print();

    double ts2 = timestamp();
    LOG1 << "sorting took " << ts2 - ts1 << " seconds";

    // check result
    if (!ss.check_order()) {
        LOG1 << "Result is not sorted!";
        abort();
    }
}
template <typename StringSet, typename LCPContainer, StringLCPSorter<StringSet, typename LCPContainer::iterator> sorter>
void TestUPtrStdString(const char* name,
                       const size_t num_strings, const size_t num_chars,
                       const std::string& letters) {

    std::default_random_engine rng(seed);

    LOG1 << "Running " << name << " on " << num_strings
         << " std::vector<std::unique_ptr<std::string>> strings";

    // vector of pointers to std::string objects
    typedef std::unique_ptr<std::string> unique_ptr_string;
    std::vector<unique_ptr_string> strings(num_strings);

    // generate random strings of length num_chars
    for (size_t i = 0; i < num_strings; ++i)
    {
        size_t slen = num_chars + (rng() >> 8) % (num_chars / 4);

        strings[i] = unique_ptr_string(new std::string(slen, 0));
        fill_random(rng, letters, strings[i]->begin(), strings[i]->end());
    }

    // run sorting algorithm
    double ts1 = timestamp();

    UPtrStdStringSet ss(strings.data(), strings.data() + strings.size());
    LCPContainer lcp(ss.size());
    sorter(ss, /* depth */ 0, /* memory */ 0, lcp.begin());

    if (0) ss.print();

    double ts2 = timestamp();
    LOG1 << "sorting took " << ts2 - ts1 << " seconds";

    // check result
    if (!ss.check_order()) {
        LOG1 << "Result is not sorted!";
        abort();
    }
    
    if (!check_lcps(ss, lcp))
    {
      LOG1 << "lcp value not correct!";
      abort();
    }
}

template <typename StringSet, StringSorter<StringSet> sorter>
void TestStringSuffixString(const char* name, const size_t num_chars,
                            const std::string& letters) {

    std::default_random_engine rng(seed);

    LOG1 << "Running " << name << " on " << num_chars
         << " std::vector<size_t> suffixes";

    // std::string text object
    std::string text(num_chars, 0);
    fill_random(rng, letters, text.begin(), text.end());

    std::vector<size_t> suffixarray;
    auto ss = StringSuffixSet::Initialize(text, suffixarray);

    // run sorting algorithm
    double ts1 = timestamp();

    sorter(ss, /* depth */ 0, /* memory */ 0);
    if (0) ss.print();

    double ts2 = timestamp();
    LOG1 << "sorting took " << ts2 - ts1 << " seconds";

    // check result
    if (!ss.check_order()) {
        LOG1 << "Result is not sorted!";
        abort();
    }
    
    
}

template <typename StringSet, typename LCPContainer, StringLCPSorter<StringSet, typename LCPContainer::iterator> sorter>
void TestStringSuffixString(const char* name, const size_t num_chars,
                            const std::string& letters) {

    std::default_random_engine rng(seed);

    LOG1 << "Running " << name << " on " << num_chars
         << " std::vector<size_t> suffixes";

    // std::string text object
    std::string text(num_chars, 0);
    fill_random(rng, letters, text.begin(), text.end());

    std::vector<size_t> suffixarray;
    auto ss = StringSuffixSet::Initialize(text, suffixarray);

    // run sorting algorithm
    double ts1 = timestamp();

    LCPContainer lcp(ss.size());
    sorter(ss, /* depth */ 0, /* memory */ 0, lcp.begin());

    if (0) ss.print();

    double ts2 = timestamp();
    LOG1 << "sorting took " << ts2 - ts1 << " seconds";

    // check result
    if (!ss.check_order()) {
        LOG1 << "Result is not sorted!";
        abort();
    }
    
    if (!check_lcps(ss, lcp))
    {
      LOG1 << "lcp value not correct!";
      abort();
    }
}
void TestFrontend(const size_t num_strings, const size_t num_chars,
                  const std::string& letters) {

    std::default_random_engine rng(seed);

    LOG1 << "Running sort_strings() on " << num_strings
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
    {
        double ts1 = timestamp();

        tlx::sort_strings(cstrings, num_strings, /* memory */ 0);

        double ts2 = timestamp();
        LOG1 << "sorting took " << ts2 - ts1 << " seconds";

        // check result
        if (!UCharStringSet(cstrings, cstrings + num_strings).check_order()) {
            LOG1 << "Result is not sorted!";
            abort();
        }
    }

    // array of const string pointers
    const uint8_t** ccstrings = new const uint8_t*[num_strings];

    for (size_t i = 0; i < num_strings; ++i)
        ccstrings[i] = cstrings[i];

    // run sorting algorithm
    {
        double ts1 = timestamp();

        tlx::sort_strings(ccstrings, num_strings, /* memory */ 0);

        double ts2 = timestamp();
        LOG1 << "sorting took " << ts2 - ts1 << " seconds";

        // check result
        if (!CUCharStringSet(ccstrings, ccstrings + num_strings).check_order()) {
            LOG1 << "Result is not sorted!";
            abort();
        }
    }

    // free memory.
    for (size_t i = 0; i < num_strings; ++i)
        delete[] cstrings[i];

    delete[] cstrings;
}

static const char* letters_alnum
    = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"
      "\xE0\xE1\xE2\xE3\xE4\xE5\xE6\xE7\xE8\xE9\xEA\xEB\xEC\xED\xEE\xEF";

// use macro because one cannot pass template functions as template parameters:
#define run_tests(func)                            \
    TestUCharString<UCharStringSet, mysorter::func>(         \
        #func, num_strings, 16, letters_alnum);    \
    TestVectorStdString<StdStringSet, mysorter::func>(       \
        #func, num_strings, 16, letters_alnum);    \
    TestUPtrStdString<UPtrStdStringSet, mysorter::func>(     \
        #func, num_strings, 16, letters_alnum);    \
    TestStringSuffixString<StringSuffixSet, mysorter::func>( \
        #func, num_strings, letters_alnum);        \
/*
#define run_lcp_tests(func)                                          \
 TestUCharString<UCharStringSet, std::vector<int>, func>(         \
        #func, num_strings, 16, letters_alnum);                      \
    TestVectorStdString<StdStringSet, std::vector<int>, func>(       \
        #func, num_strings, 16, letters_alnum);                      \
    TestUPtrStdString<UPtrStdStringSet, std::vector<int>, func>(     \
        #func, num_strings, 16, letters_alnum);                      \
    TestStringSuffixString<StringSuffixSet, std::vector<int>, func>( \
        #func, num_strings, letters_alnum);                          \
*/
#define run_unified_tests(func)     \
    TestLcpStringPtr<dss_schimek::LcpStringPtr<unsigned char>, func>( \
        #func, num_strings, 16u, letters_alnum); 

void test_all(const size_t num_strings) {
    if (num_strings <= 1024) {
      run_tests(insertion_sort);
    }
    run_tests(multikey_quicksort);
    run_tests(radixsort_CE0);
    run_tests(radixsort_CE2);
    run_tests(radixsort_CE3);
    run_tests(radixsort_CI2);
    run_tests(radixsort_CI3);

    TestFrontend(num_strings, 16, letters_alnum);
}
void test_all_lcp_unified(const size_t num_strings) {
  if (num_strings <= 1024) {
   run_unified_tests(mysorter::insertion_sort);
  }
  letters_alnum = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"
      "\xE0\xE1\xE2\xE3\xE4\xE5\xE6\xE7\xE8\xE9\xEA\xEB\xEC\xED\xEE\xEF";

  std::cout << "original test string: " << letters_alnum << std::endl;
  run_unified_tests(mysorter::multikey_quicksort);
  run_unified_tests(mysorter::radixsort_CE0);
  run_unified_tests(mysorter::radixsort_CE2);
  run_unified_tests(mysorter::radixsort_CE3);
  run_unified_tests(mysorter::radixsort_CI2);
  run_unified_tests(mysorter::radixsort_CI3);

  letters_alnum = "a";
  std::cout << "test string: " << letters_alnum << std::endl;

  run_unified_tests(mysorter::multikey_quicksort);
  run_unified_tests(mysorter::radixsort_CE0);
  run_unified_tests(mysorter::radixsort_CE2);
  run_unified_tests(mysorter::radixsort_CE3);
  run_unified_tests(mysorter::radixsort_CI2);
  run_unified_tests(mysorter::radixsort_CI3);
  letters_alnum = "abc";
  std::cout << "test string: " << letters_alnum << std::endl;
    
  run_unified_tests(mysorter::multikey_quicksort);
  run_unified_tests(mysorter::radixsort_CE0);
  run_unified_tests(mysorter::radixsort_CE2);
  run_unified_tests(mysorter::radixsort_CE3);
  run_unified_tests(mysorter::radixsort_CI2);
  run_unified_tests(mysorter::radixsort_CI3);

  letters_alnum = "abcde";
  std::cout << "test string: " << letters_alnum << std::endl;

  run_unified_tests(mysorter::multikey_quicksort);
  run_unified_tests(mysorter::radixsort_CE0);
  run_unified_tests(mysorter::radixsort_CE2);
  run_unified_tests(mysorter::radixsort_CE3);
  run_unified_tests(mysorter::radixsort_CI2);
  run_unified_tests(mysorter::radixsort_CI3);

}
/*void test_all_lcp(const size_t num_strings) {
  //if (num_strings <= 1024) {
  // run_lcp_tests(insertion_sort);
  //}
  letters_alnum = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"
      "\xE0\xE1\xE2\xE3\xE4\xE5\xE6\xE7\xE8\xE9\xEA\xEB\xEC\xED\xEE\xEF";

  std::cout << "original test string: " << letters_alnum << std::endl;
  //run_lcp_tests(multikey_quicksort);
  //run_lcp_tests(radixsort_CE0);
  //run_lcp_tests(radixsort_CE2);
  //run_lcp_tests(radixsort_CE3);
  //run_lcp_tests(radixsort_CI2);
  //run_lcp_tests(radixsort_CI3);

  letters_alnum = "a";
  std::cout << "test string: " << letters_alnum << std::endl;
  ////run_lcp_tests(multikey_quicksort);
  //run_lcp_tests(radixsort_CE0);
  //run_lcp_tests(radixsort_CE2);
  //run_lcp_tests(radixsort_CE3);
  //run_lcp_tests(radixsort_CI2);
  run_lcp_tests(radixsort_CI3);

  letters_alnum = "abc";
  std::cout << "test string: " << letters_alnum << std::endl;
  //run_lcp_tests(multikey_quicksort);
  //run_lcp_tests(radixsort_CE0);
  //run_lcp_tests(radixsort_CE2);
  //run_lcp_tests(radixsort_CE3);
  //run_lcp_tests(radixsort_CI2);
  //run_lcp_tests(radixsort_CI3);
  
}*/
int main() {
    test_all(16);
    test_all(256);
    test_all(65550);
    //if (tlx_more_tests) {
    //    test_all(1024 * 1024);
    //    // test_all(16 * 1024 * 1024);
    //}
    test_all_lcp_unified(16);
    test_all_lcp_unified(256);
    test_all_lcp_unified(65550);
    test_all_lcp_unified(655550);
    test_all_lcp_unified(6555550);
    //test_all_lcp(16);
    //test_all_lcp(65550);
    //test_all(65550);
    //test_all_lcp(16);
    //test_all_lcp(256);
    //test_all_lcp(65550);
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
