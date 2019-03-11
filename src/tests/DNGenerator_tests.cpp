#include <algorithm>
#include <vector>

#include <tlx/die/core.hpp>
#include <tlx/sort/strings/radix_sort.hpp>

#include "strings/stringset.hpp"
#include "strings/stringtools.hpp"

#include "util/measuringTool.hpp"
#include "util/random_string_generator.hpp"

#include <functional>

#include "mpi/environment.hpp"

namespace dss_schimek {
  namespace tests {

    
    template <typename StringSet>
    size_t computeD(StringSet sortedStringSet) {
      using String = typename StringSet::String;
      using CharIt = typename StringSet::CharIterator;
      if (sortedStringSet.size() < 2)
        return 0;
      size_t D = 0;
      std::vector<size_t> lcps(sortedStringSet.size(), 0);
      auto begin = sortedStringSet.begin();
      for (size_t i = 1; i < sortedStringSet.size(); ++i) {
        String prevString = sortedStringSet[begin + i - 1];
        String curString = sortedStringSet[begin + i];
        CharIt prevChars = sortedStringSet.get_chars(prevString, 0);
        CharIt curChars = sortedStringSet.get_chars(curString, 0);
        lcps[i] = calc_lcp(prevChars, curChars);
      }

      for (size_t i = 0; i + 1 < sortedStringSet.size(); ++i) {
        size_t prevLcp = lcps[i];
        size_t nextLcp = lcps[i + 1];
        D += std::max(prevLcp, nextLcp) + 1;
      }
      D += lcps.back() + 1;
      return D;
    }

    void DNRatioGenerator_test(const size_t size, const size_t stringLength, const double dToNRatio, const double epsilon = 0.01) {
      using namespace dss_schimek;
      using StringSet = UCharLengthStringSet;
      using String = StringSet::String;
      const bool verbose = false;

      if (verbose) 
        std::cout << "size: " << size << " stringLength " << stringLength << " dToNratio: " << dToNRatio << std::endl;

      dss_schimek::DNRatioGenerator<StringSet> generator(size, stringLength, dToNRatio);
      auto stringPtr = generator.make_string_lcp_ptr();
      StringSet ss = stringPtr.active();

      tlx::sort_strings_detail::radixsort_CI3(stringPtr, 0, 0);
      const size_t D = computeD(ss);
      const double acutalDToNRatio = static_cast<double>(D) / ( size * stringLength);

      const bool allStringsHaveCorrectLength = 
        std::all_of(ss.begin(), ss.end(), [&](const String& str) {
            return ss.get_length(str) == stringLength;
            });

      if (verbose) {
        std::cout << "D: " << D << " N: " << size*stringLength << std::endl;
        std::cout << "D: " << static_cast<double>(D) / (size*stringLength) << std::endl;
      }

  
      // +++++++++++
      // tests
      // +++++++++++
      tlx_die_unless(ss.size() == size);
      tlx_die_unless(allStringsHaveCorrectLength);
      tlx_die_unless(std::abs(acutalDToNRatio - dToNRatio) < epsilon);
    }
  }
 }

int main() {
  using namespace dss_schimek::tests;
  const std::vector<size_t> sizes = {1000, 10000, 100000, 1000000};
  const std::vector<size_t> stringLengths = {100,250, 500};
  const double epsilon = 0.01;
  const double relaxedEpsilon = 0.075; // TODO calculate exact bounds
  for (const auto& size : sizes) {
    std::cout << "start tests with size: " << size << " ";
    for (const auto& stringLength : stringLengths) {
      std::cout << " stringLength: " << stringLength << std::endl;
      DNRatioGenerator_test(size, stringLength, 0.0,  relaxedEpsilon);
      DNRatioGenerator_test(size, stringLength, 0.2,  epsilon);
      DNRatioGenerator_test(size, stringLength, 0.25, epsilon);
      DNRatioGenerator_test(size, stringLength, 0.45, epsilon);
      DNRatioGenerator_test(size, stringLength, 0.76, epsilon);
      DNRatioGenerator_test(size, stringLength, 0.95, epsilon);
      DNRatioGenerator_test(size, stringLength, 1.0,  epsilon);
    }
  }
  std::cout << "completed tests successfully" << std::endl;
}
