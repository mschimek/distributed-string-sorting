#include "sorter/distributed/lcp_aware_number_string_splitting.hpp"

#include "util/random_string_generator.hpp"
#include "mpi/synchron.hpp"

int main() {
  using namespace dss_schimek;

  dsss::mpi::environment env;
  RandomStringLcpContainer<unsigned char> rand_container(10);
  UCharStringSet ss(rand_container.strings(), rand_container.strings() + rand_container.size()); 
  StringLcpPtr<UCharStringSet> rand_string_ptr(ss, rand_container.lcp_array());

  dss_schimek::mpi::execute_in_order([&](){
      for (volatile size_t i = 0; i < 1000000; ++i);
      std::cout << "rank " << env.rank() << std::endl;
      ss.print();
      });
  env.barrier();
  StringLcpContainer<unsigned char> sorted_string_cont = 
    merge_sort(rand_string_ptr, std::move(rand_container));


  UCharStringSet sorted_strings_set(sorted_string_cont.strings(),
      sorted_string_cont.strings() + sorted_string_cont.size());
  StringLcpPtr<UCharStringSet> strptr_res(sorted_strings_set, sorted_string_cont.lcp_array());
  const bool is_complete_and_sorted = dss_schimek::is_complete_and_sorted(strptr_res,
      rand_container.char_size(),
      sorted_string_cont.char_size(),
      rand_container.size(),
      sorted_string_cont.size()); 

  dss_schimek::mpi::execute_in_order([&](){
      for (volatile size_t i = 0; i < 1000000; ++i);
      std::cout << "rank " << env.rank() << std::endl;
      sorted_strings_set.print();
      });

  if (env.rank() == 0)
    std::cout << "res: " << is_complete_and_sorted << std::endl;

  env.finalize();
}
