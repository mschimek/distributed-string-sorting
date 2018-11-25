#include "sorter/distributed/lcp_aware_number_string_splitting.hpp"

#include "util/random_string_generator.hpp"
#include "mpi/synchron.hpp"

int main() {
  using namespace dss_schimek;

  dsss::mpi::environment env;
  RandomStringLcpContainer<UCharStringSet> rand_container(10);
  StringLcpPtr<UCharStringSet> rand_string_ptr = 
    rand_container.make_string_lcp_ptr();

  dss_schimek::mpi::execute_in_order([&](){
      for (volatile size_t i = 0; i < 1000000; ++i);
      std::cout << "rank " << env.rank() << std::endl;
      rand_string_ptr.active().print();
      });
  env.barrier();

  StringLcpContainer<UCharStringSet> sorted_string_cont = 
    merge_sort(rand_string_ptr, std::move(rand_container));


  StringLcpPtr<UCharStringSet> sorted_strptr = 
    sorted_string_cont.make_string_lcp_ptr();

  const bool is_complete_and_sorted = dss_schimek::is_complete_and_sorted(sorted_strptr,
      rand_container.char_size(),
      sorted_string_cont.char_size(),
      rand_container.size(),
      sorted_string_cont.size()); 

  dss_schimek::mpi::execute_in_order([&](){
      for (volatile size_t i = 0; i < 1000000; ++i);
      std::cout << "rank " << env.rank() << std::endl;
      sorted_strptr.active().print();
      });

  if (env.rank() == 0)
    std::cout << "res: " << is_complete_and_sorted << std::endl;

  env.finalize();
}
