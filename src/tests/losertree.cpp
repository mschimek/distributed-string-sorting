#include "merge/bingmann-lcp_losertree.hpp"
#include "strings/stringtools.hpp"
#include "strings/stringptr.hpp"
#include "strings/stringset.hpp"
#include "util/random_string_generator.hpp"
#include "sorter/local/strings/insertion_sort_unified.hpp"

int main() {
  using namespace dss_schimek;
  size_t num_strings = 10;
  RandomStringLcpContainer<unsigned char> rand_container(num_strings);
  UCharStringSet ss(rand_container.strings(), rand_container.strings() + rand_container.size()); 
  StringLcpPtr<UCharStringSet> rand_string_ptr(ss, rand_container.lcp_array());

  insertion_sort(rand_string_ptr, 0, 0);
  std::cout << "sorted sequence " << std::endl;
  ss.print();

  std::vector<std::pair<size_t, size_t>> pairs;
  pairs.emplace_back(0,3);
  pairs.emplace_back(3,1);
  pairs.emplace_back(6,4);
  pairs.emplace_back(0,0);
  pairs.emplace_back(0,0);
  stringtools::LcpStringPtr lcp_string_ptr(rand_container.strings(), 
      rand_container.lcp_array(), 
      rand_container.size());
  bingmann::LcpStringLoserTree<5> lt(lcp_string_ptr, pairs.data());
  std::vector<unsigned char*> output_strings(num_strings, 0);
  std::vector<size_t> output_lcp(num_strings, 0);
  stringtools::LcpStringPtr output_string_ptr(output_strings.data(), 
      output_lcp.data(),
      num_strings);
  lt.writeElementsToStream(output_string_ptr, num_strings);

  dss_schimek::print_str(output_strings.data(), output_strings.size());


  
}
