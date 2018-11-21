#include "sorter/distributed/lcp_aware_number_string_splitting.hpp"

#include "util/random_string_generator.hpp"
#include "mpi/synchron.hpp"

int main() {
  using namespace dss_schimek;
  
  dsss::mpi::environment env;
  RandomStringLcpContainer<unsigned char> rand_container(500);
  UCharStringSet ss(rand_container.strings(), rand_container.strings() + rand_container.size()); 
  StringLcpPtr<UCharStringSet> rand_string_ptr(ss, rand_container.lcp_array());

  dss_schimek::mpi::execute_in_order([&](){
      std::cout << "rank: " << env.rank() << std::endl;
      ss.print();
      }, env);
  env.barrier();
  dss_schimek::merge_sort(rand_string_ptr, rand_container);

  env.finalize();
}
