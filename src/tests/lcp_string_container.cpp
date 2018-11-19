#include "strings/stringptr.hpp"
#include "strings/stringset.hpp"
#include "util/random_string_generator.hpp"
#include <sorter/local/strings/insertion_sort_unified.hpp>
#include "mpi/alltoall.hpp"

using namespace dss_schimek;
int main() {
  dsss::mpi::environment env;

  PrefixNumberStringLcpContainer<unsigned char> rand_container(10, env.rank() + 65);
  UCharStringSet ss(rand_container.strings(), rand_container.strings() + rand_container.size()); 
  StringLcpPtr<UCharStringSet> rand_string_ptr(ss, rand_container.lcp_array());
  rand_string_ptr.active().print();
  
  dss_schimek::insertion_sort(rand_string_ptr, 0, 0);
  for(size_t i = 0; i < env.size(); ++i)
  {
    env.barrier();
    if (i == env.rank())
    {
      std::cout << "sorted rank: " << i << std::endl;
      rand_string_ptr.active().print();
    }
  }
  
  env.barrier();
  StringLcpContainer<unsigned char> recv;
  
  if (env.rank() == 0)
    recv = dsss::mpi::alltoallv(rand_container, {1, 1, 1, 1, 1}, env);
  else
    recv = dsss::mpi::alltoallv(rand_container, {1, 0, 0, 0 , 0}, env);

  //UCharStringSet ss_(recv.strings(), recv.strings() + recv.size());
  //StringLcpPtr<UCharStringSet> receivedSet(ss_, recv.lcp_array());
  //for(size_t i = 0; i < env.size(); ++i)
  //{
  //  env.barrier();
  //  if (i == env.rank())
  //  {
  //    std::cout << "rank: " << i << std::endl;
  //    receivedSet.active().print();
  //  }
  //}
  env.finalize();
  
}
