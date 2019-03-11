#include <vector>
#include <algorithm>
#include <functional>

#include <tlx/siphash.hpp>
#include <tlx/sort/strings/radix_sort.hpp>

#include "strings/stringcontainer.hpp"
#include "strings/stringset.hpp"
#include "sorter/distributed/bloomfilter.hpp"
#include "util/random_string_generator.hpp"
#include "util/measuringTool.hpp"

namespace dss_schimek {
  namespace tests {
    void GolombEncodingPipeline_test() {
      using Result = std::vector<std::vector<size_t>>;
      dsss::mpi::environment env;
      const size_t sizePerPE = 100;

      // All PE get same messages, namely {0, 1, 2, ... , sizePerPE - 1} from each PE
      std::vector<size_t> hashes;
      std::vector<size_t> sample(sizePerPE);
      std::iota(sample.begin(), sample.end(), 0);
      for (size_t i = 0; i < env.size(); ++i)
        std::copy_n(sample.begin(), sample.size(), std::back_inserter(hashes));
      std::vector<size_t> sendCounts(env.size(), sizePerPE);

      Result recvData = AllToAllHashValuesPipeline::alltoallv(hashes, sendCounts);

      tlx_die_unless(recvData.size() == env.size());
      for (size_t i = 0; i < recvData.size(); ++i) {
        const auto& curVec = recvData[i];
        std::cout << "rank: " << env.rank() << " curVec size: " << curVec.size() <<  std::endl;
        for (const auto& elem : curVec) {
          std::cout << "rank: " << env.rank() << " partner:  " << i <<  " elem: " << elem << std::endl;
        }

        if (!(curVec == sample)) {
          std::cout << "rank: " << env.rank() << " partner:  " << i <<  " failed" << std::endl;
        }
        tlx_die_unless(curVec == sample);
      }
    }    
  }
}

int main() {
  using namespace dss_schimek::tests;
  dsss::mpi::environment env;
  std::cout << "start tests" << std::endl;
  GolombEncodingPipeline_test(); 
  std::cout << "tests completed successfully" << std::endl;
  env.finalize();
}
