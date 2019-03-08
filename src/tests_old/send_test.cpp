#include <vector>
#include "strings/stringtools.hpp"
#include "strings/stringptr.hpp"
#include "strings/stringset.hpp"
#include "util/random_string_generator.hpp"
#include "util/timer.hpp"
#include "mpi/alltoall.hpp"
#include "mpi/synchron.hpp"
#include "mpi/byte_encoder.hpp"
#include "sorter/local/strings/insertion_sort_unified.hpp"

  using namespace dss_schimek;
  using StringSet = UCharLengthStringSet;

 template<typename Lhs, typename Rhs>
 bool isEqual(Lhs& lhs, Rhs& rhs) {
        return(lhs.raw_strings() == rhs.raw_strings()) && (lhs.lcps() == rhs.lcps());
      }
int main() {
  using namespace dss_schimek;
  using namespace dss_schimek::mpi;
  using namespace dsss::mpi;
  using StringSet = UCharLengthStringSet;
  dsss::mpi::environment env;
  size_t numStrings = 100;

  RandomStringLcpContainer<StringSet> randContainer(numStrings);
  dss_schimek::StringLcpPtr strptr = randContainer.make_string_lcp_ptr();
  insertion_sort(strptr, 0, 0);

  dss_schimek::mpi::execute_in_order([&] () {
      std::cout << "rank: " << env.rank() << std::endl;
      strptr.active().print();
      });

  
  Timer timer1(" ");
  Timer timer2(" ");
  Timer timer3(" ");
  Timer timer4(" ");
  EmptyTimer emptyTimer;
  std::vector<size_t> sendCounts{2, 2, 95, 1};
  std::cout << "call 0" << std::endl;
  auto recvContainer2 = 
    dsss::mpi::alltoallv<StringSet, 
                         AllToAllvCombined<AllToAllvSmall>, 
                         SequentialDelayedByteEncoder,
                         EmptyTimer>(randContainer, sendCounts, emptyTimer);
  std::cout << "call 1" << std::endl;
  auto recvContainer3 = 
    dsss::mpi::alltoallv<StringSet,
                         AllToAllvDirectMessages, 
                         SequentialByteEncoder,
                         Timer>(randContainer, sendCounts, timer2);
  std::cout << "call 2" << std::endl;
  auto recvContainer4 = 
    dsss::mpi::alltoallv<StringSet,
                         AllToAllvSmall,
                         InterleavedByteEncoder,
                         Timer>(randContainer, sendCounts, timer3);
  std::cout << "call 3" << std::endl;
  auto recvContainerRef = 
    dsss::mpi::alltoallv<StringSet,
                         AllToAllvSmall,
                         EmptyByteEncoder,
                         Timer>(randContainer, sendCounts, timer4);
  
  //if (!(recvContainerRef ==  recvContainer2)) {
  //  std::cout << "alltoall failed! at rank " << env.rank() << std::endl;
  //  //std::abort();
  //}
  if (!(recvContainerRef ==  recvContainer3)) {
    std::cout << "alltoall failed! at rank " << env.rank() << std::endl;
    //std::abort();
  }
  if (!(recvContainerRef ==  recvContainer4)) {
    std::cout << "alltoall failed! at rank " << env.rank() << std::endl;
    //std::abort();
  }
 std::stringstream buffer;
  timer1.writeToStream(buffer);
  if (env.rank() == 0) {
    std::cout << buffer.str() << std::endl;
  }
  //dss_schimek::mpi::execute_in_order([&] () {
  //    std::cout << "rank: " << env.rank() << std::endl;
  //    auto recvStringSet = recvContainer3.make_string_set();
  //    recvStringSet.print();
  //    });



  env.finalize();
}
