#include "sorter/distributed/lcp_aware_number_string_splitting.hpp"

#include "util/random_string_generator.hpp"
#include "mpi/synchron.hpp"

int main() {
  using namespace dss_schimek;

  using StringSet = UCharLengthStringSet;

  dsss::mpi::environment env;
  constexpr size_t size = 20000000;
  RandomStringLcpContainer<StringSet> rand_container(size);
  StringLcpPtr<StringSet> rand_string_ptr = 
    rand_container.make_string_lcp_ptr();

  //dss_schimek::mpi::execute_in_order([&](){
  //    for (volatile size_t i = 0; i < 1000000; ++i);
  //    std::cout << "rank " << env.rank() << std::endl;
  //    rand_string_ptr.active().print();
  //    });

  env.barrier();
  double start_time = MPI_Wtime();
  DistributedMergeSort<StringLcpPtr<StringSet>> sorter;
  StringLcpContainer<StringSet> sorted_string_cont = 
    sorter.sort(rand_string_ptr, std::move(rand_container));

  double end_time = MPI_Wtime();
  double res = end_time - start_time;
  double overall_res = dsss::mpi::allreduce_max(res);

  RandomStringLcpContainer<StringSet> rand_container_2(size);
  StringLcpPtr<StringSet> rand_string_ptr_2 = 
    rand_container_2.make_string_lcp_ptr();

  double start = MPI_Wtime();
  dss_schimek::multikey_quicksort(rand_string_ptr_2, 0, 0);
  double end = MPI_Wtime();

  if (env.rank() == 0)
  std::cout << "local sort " << end - start << std::endl;

  StringLcpPtr<StringSet> sorted_strptr = 
    sorted_string_cont.make_string_lcp_ptr();

  const bool is_complete_and_sorted = dss_schimek::is_complete_and_sorted(sorted_strptr,
      rand_container.char_size(),
      sorted_string_cont.char_size(),
      rand_container.size(),
      sorted_string_cont.size()); 

  //dss_schimek::mpi::execute_in_order([&](){
  //    for (volatile size_t i = 0; i < 1000000; ++i);
  //    std::cout << "rank " << env.rank() << std::endl;
  //    sorted_strptr.active().print();
  //    });

  if (env.rank() == 0) {
    std::cout << "res: " << is_complete_and_sorted << std::endl;
    std::cout << "time in seconds: " << overall_res << std::endl;
  }
  env.finalize();
}
