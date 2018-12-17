#include <vector>
#include "strings/stringtools.hpp"
#include "strings/stringptr.hpp"
#include "strings/stringset.hpp"
#include "util/random_string_generator.hpp"
#include "mpi/alltoall.hpp"
#include "mpi/synchron.hpp"
#include "sorter/local/strings/insertion_sort_unified.hpp"

int main() {
  using namespace dss_schimek;
  using StringSet = UCharLengthStringSet;
  dsss::mpi::environment env;
  size_t numStrings = 10;

  RandomStringLcpContainer<StringSet> randContainer(numStrings);
  dss_schimek::StringLcpPtr strptr = randContainer.make_string_lcp_ptr();
  insertion_sort(strptr, 0, 0);

  dss_schimek::mpi::execute_in_order([&] () {
      std::cout << "rank: " << env.rank() << std::endl;
      strptr.active().print();
      });

  
  std::vector<size_t> sendCounts{3,7};
  dsss::mpi::alltoallv_2(randContainer, sendCounts);
  


  env.finalize();
}
