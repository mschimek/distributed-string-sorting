#include "util/random_string_generator.hpp"
#include "mpi/environment.hpp"
#include "mpi/synchron.hpp"
#include "mpi/alltoall.hpp"
#include "util/timer.hpp"
#include "sorter/local/strings/multikey_quicksort_unified.hpp"

#include <iostream>

int main() {
  using namespace dss_schimek;
  dsss::mpi::environment env;
  size_t numOfStrings = 100;

  using StringSet = UCharLengthStringSet;
  using MPISendRoutine = dsss::mpi::AllToAllvSmall;
  using Timer = Timer;
  using ByteEncoder = EmptyLcpByteEncoderMemCpy;
  
  dss_schimek::SkewedRandomStringLcpContainer<StringSet> rand_container(numOfStrings);
  StringLcpPtr local_string_ptr = rand_container.make_string_lcp_ptr();
  dss_schimek::multikey_quicksort(local_string_ptr, 0, 0);
  StringSet ss = rand_container.make_string_set();

  mpi::execute_in_order([&]() {
      std::cout << "\nrank " << env.rank() << std::endl;
      rand_container.make_string_set().print();
      
      env.barrier();
      for (volatile size_t i = 0; i < 1000000; ++i);
      env.barrier();
      });

  dsss::mpi::AllToAllStringImpl<StringSet, MPISendRoutine, ByteEncoder, Timer> alltoall;
  Timer timer("");
  std::vector<size_t> sendCounts = {5, 5};
  StringLcpContainer<StringSet> recvCont = alltoall.alltoallv(rand_container, sendCounts, timer);

  mpi::execute_in_order([&]() {
      env.barrier();
      std::cout << "rank " << env.rank() << std::endl;
      recvCont.make_string_set().print();
      });
  
  env.finalize();
}
