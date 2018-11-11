//#include "../external/dsss/string_sorting/distributed/merge_sort.hpp" 
//#include "../external/dsss/util/string_set.hpp"
//#include "util/random_string_generator.hpp"
//
//#include "sequential/bingmann-mkqs.hpp"
//#include "sequential/bingmann-radix_sort.hpp"
//
//#include "string_sorting/distributed/merge_sort.hpp"

#include "utils/draw_from_distribution.hpp"
//#include "strings/stringtools.hpp"
#include <tlx/sort/strings/insertion_sort.hpp>
#include <tlx/sort/strings/multikey_quicksort.hpp>
#include "sorter/isort.hpp"
/*void sample_sort_kurpicz(dsss::string_set& local_string_set) 
{            
  dsss::sample_sort::sample_sort<bingmann::bingmann_msd_CE0>(local_string_set);
}

dsss::random_string_set create_random_string_set(const size_t size, 
  const size_t min_length, const size_t max_length)
{
  return dsss::random_string_set(size, min_length, max_length);
}*/

int main()
{
  const char* str_a = "aadef";
  const char* str_b = "abc";
  const char* str_c = "abeeeeeee";
  const char* str_d = "abeeeeefffff";
  std::vector<const char*> strings{str_a, str_b, str_c, str_d}; 
  dss_schimek::GenericCharStringSet<const char> ss(std::make_pair(strings.data(), 2));
  size_t* lcp_values = new size_t[4];
  //stringtools::LcpStringPtrSet<const char> ss(strings.data(), strings.data() + 4, lcp_values);
  std::cout << ss[ss.begin() + 0] << " " << ss[ss.begin() + 1] << std::endl;
  mysorter::insertion_sort(ss, 0, 0);
  for (size_t i = 0; i < ss.size(); ++i)
  std::cout << ss[ss.begin() + i] << " " ;
  std::cout << std::endl;
  for(size_t i = 0; i < ss.size(); ++i)
    std::cout << ss.get_lcp(i) << " ";
  std::cout << std::endl;
  constexpr size_t size = 200;
  constexpr size_t min_string_length = 10;
  constexpr size_t max_string_length = 20;
}
