#include "sorter/distributed/bloomfilter.hpp"
#include <iostream>

#include "util/random_string_generator.hpp"
#include "util/timer.hpp"
#include "mpi/synchron.hpp"
#include <map>

#include <tlx/cmdline_parser.hpp>


int main() {
  using namespace dss_schimek;
  using StringSet = UCharLengthStringSet;
  using StringLcpPtr = typename tlx::sort_strings_detail::StringLcpPtr<StringSet, size_t>;

  dss_schimek::SkewedRandomStringLcpContainer<StringSet> generatedContainer(44);
  StringLcpPtr rand_string_ptr = 
    generatedContainer.make_string_lcp_ptr();

  tlx::sort_strings_detail::radixsort_CI3(rand_string_ptr, 0, 0);
  


  //for (size_t i = 0; i < results.size(); ++i) {
  //  std::cout << "rank: " << i << std::endl;
  //  for (size_t j = 0; j < results[i].size(); ++j) 
  //    std::cout << "\t\t" <<  j << " " << results[i][j] << std::endl;
  }
}
