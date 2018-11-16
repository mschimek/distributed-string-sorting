#include "strings/stringptr.hpp"
#include "util/random_string_generator.hpp"
#include <tlx/sort/strings/insertion_sort_unified.hpp>
#include "mpi/alltoall.hpp"

using namespace dss_schimek;
int main() {
  dsss::mpi::environment env;

  PrefixNumberLcpStringContainer<unsigned char> rand_container(10, env.rank() + 65);
  LcpStringPtr<unsigned char> rand_string_set(rand_container);
  rand_string_set.print();
  
  mysorter::insertion_sort(rand_string_set, 0, 0);
  //for(size_t i = 0; i < env.size(); ++i)
  //{
  //  env.barrier();
  //  if (i == env.rank())
  //  {
  //    std::cout << "sorted rank: " << i << std::endl;
  //    rand_string_set.print();
  //  }
  //}
  
  env.barrier();
  LcpStringContainer<unsigned char> recv;
  
  if (env.rank() == 0)
    recv = dsss::mpi::alltoallv(rand_container, {1, 1, 1, 1, 1}, env);
  else
    recv = dsss::mpi::alltoallv(rand_container, {10, 0, 0, 0 , 0}, env);

  LcpStringPtr<unsigned char> receivedSet(recv);
  for(size_t i = 0; i < env.size(); ++i)
  {
    env.barrier();
    if (i == env.rank())
    {
      std::cout << "rank: " << i << std::endl;
      receivedSet.print();
    }
  }
  env.finalize();
  
}
