#include "mpi/is_sorted.hpp"
#include "mpi/environment.hpp"
#include "strings/stringptr.hpp"
#include "strings/stringset.hpp"
#include "util/random_string_generator.hpp"
#include <sorter/local/strings/insertion_sort_unified.hpp>

using namespace dss_schimek;
int main()
{
  dsss::mpi::environment env;

  if(env.rank() != 0) {
  PrefixNumberStringLcpContainer<unsigned char> rand_container(10, env.rank() + 65);
  UCharStringSet ss(rand_container.strings(), rand_container.strings() + rand_container.size()); 
  StringLcpPtr<UCharStringSet> rand_string_ptr(ss, rand_container.lcp_array());
  insertion_sort(rand_string_ptr, 0, 0);

  std::cout << "is sorted: " << is_sorted(rand_string_ptr) << std::endl;
  }
  else
  {
    PrefixNumberStringLcpContainer<unsigned char> rand_container(0, env.rank() + 65);
  UCharStringSet ss(rand_container.strings(), rand_container.strings() + rand_container.size()); 
  StringLcpPtr<UCharStringSet> rand_string_ptr(ss, rand_container.lcp_array());
  std::cout << ss.size() << std::endl;
  insertion_sort(rand_string_ptr, 0, 0);

  
  std::cout << "is sorted: " << is_sorted(rand_string_ptr) << std::endl;
  }
  env.finalize();

}
