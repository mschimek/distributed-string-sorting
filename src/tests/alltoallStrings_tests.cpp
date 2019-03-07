#include <iostream>
#include <algorithm>

#include <tlx/die/core.hpp>

#include "mpi/environment.hpp"
#include "mpi/alltoall.hpp"

#include "util/measuringTool.hpp"

#include "strings/stringset.hpp"
#include "strings/stringptr.hpp"
#include "strings/stringcontainer.hpp"

namespace dss_schimek {
  namespace tests {
    std::vector<unsigned char> generateRawStrings(const size_t sizePerPE, size_t commonPrefixLength, size_t differentSuffixLength) {
      dsss::mpi::environment env;

      std::vector<unsigned char> rawStrings;

      for (size_t peIndex = 0; peIndex < env.size(); ++peIndex) {
        for (size_t i = 1; i < sizePerPE + 1; ++i) {
          for (size_t j = 0; j < commonPrefixLength; ++j)
            rawStrings.push_back(env.rank() + 1);
          rawStrings.push_back(i);
          for (size_t j = 0; j < differentSuffixLength; ++j)
            rawStrings.push_back(env.rank() + 1);
          rawStrings.push_back(0);
        }
      }
      return rawStrings;
    }

    std::vector<unsigned char> getExpectedRawStringsForPrefixCompression(size_t sizePerPE, size_t commonPrefixLength, size_t differentSuffixLength) {
      dsss::mpi::environment env;
      std::vector<unsigned char> expectedRawStrings;
      for (size_t peIndex = 0; peIndex < env.size(); ++peIndex) {
        for (size_t i = 0; i < sizePerPE; ++i) {
          std::fill_n(std::back_inserter(expectedRawStrings), commonPrefixLength, peIndex + 1);
          expectedRawStrings.push_back(i + 1);
          std::fill_n(std::back_inserter(expectedRawStrings), differentSuffixLength, peIndex + 1);
          expectedRawStrings.push_back(0);
        }
      }
      return expectedRawStrings;
    }
    
    std::vector<unsigned char> getExpectedRawStrings(size_t sizePerPE, size_t commonPrefixLength, size_t differentSuffixLength) {
      dsss::mpi::environment env;
      std::vector<unsigned char> expectedRawStrings;
      for (size_t peIndex = 0; peIndex < env.size(); ++peIndex) {
        for (size_t i = 0; i < sizePerPE; ++i) {
          std::fill_n(std::back_inserter(expectedRawStrings), commonPrefixLength, peIndex + 1);
          expectedRawStrings.push_back(i + 1);
          std::fill_n(std::back_inserter(expectedRawStrings), differentSuffixLength, peIndex + 1);
          expectedRawStrings.push_back(0);
        }
      }
      return expectedRawStrings;
    }

    std::vector<unsigned char> getExpectedRawStringsForPrefixDoubling(size_t sizePerPE, size_t commonPrefixLength) {
      dsss::mpi::environment env;
      std::vector<unsigned char> expectedRawStrings;
      for (size_t peIndex = 0; peIndex < env.size(); ++peIndex) {
        std::fill_n(std::back_inserter(expectedRawStrings), commonPrefixLength, peIndex + 1);
        for (size_t i = 0; i < sizePerPE; ++i) {
          expectedRawStrings.push_back(i + 1);
          expectedRawStrings.push_back(0);
        }
      }
      return expectedRawStrings;
    }
    
    std::vector<unsigned char> getExpectedRawStrings_DuplicatesForPrefixDoubling(size_t sizePerPE, size_t commonPrefixLength) {
      dsss::mpi::environment env;
      std::vector<unsigned char> expectedRawStrings;
      for (size_t peIndex = 0; peIndex < env.size(); ++peIndex) {
        std::fill_n(std::back_inserter(expectedRawStrings), commonPrefixLength, peIndex + 1);
        std::fill_n(std::back_inserter(expectedRawStrings), sizePerPE, 0);
      }
      return expectedRawStrings;
    }

    std::vector<size_t> getExpectedLcpValues(size_t sizePerPE, size_t commonPrefixLength) {
      dsss::mpi::environment env;
      std::vector<size_t> expectedLcps;
      for (size_t peIndex = 0; peIndex < env.size(); ++peIndex) {
        expectedLcps.push_back(0);
        std::fill_n(std::back_inserter(expectedLcps), sizePerPE - 1, commonPrefixLength);
      }
      return expectedLcps;
    }

    template <typename StringSet>
    struct Data {
      dss_schimek::StringLcpContainer<StringSet> container;
      std::vector<size_t> sendCounts;
      size_t size;
      size_t sizePerPE;
      size_t commonPrefixLength;
    };

    template <typename StringSet>
      Data<StringSet> commonSetup() {
        using namespace dsss::mpi;
        using StringLcpPtr = typename tlx::sort_strings_detail::StringLcpPtr<StringSet, size_t>;
        dsss::mpi::environment env;

        const size_t commonPrefixLength = 10;
        const size_t differentSuffixLength = 20;
        const size_t sizePerPE = 5;
        const size_t size = sizePerPE * env.size();
        StringLcpContainer<StringSet> container(generateRawStrings(sizePerPE, commonPrefixLength, differentSuffixLength));
        StringLcpPtr localStringPtr = container.make_string_lcp_ptr();

        std::vector<size_t> sendCounts(env.size(), sizePerPE);

        // Set lcp values
        for (size_t i = 0; i < size; ++i)
          localStringPtr.set_lcp(i, commonPrefixLength);
        localStringPtr.set_lcp(0, 0);

      return {std::move(container), std::move(sendCounts), size, sizePerPE, commonPrefixLength};
    }
    
    void AllToAllStringImplPrefixDoubling_test(bool simulateDuplicates) {
      using namespace dsss::mpi;
      using StringSet = UCharLengthStringSet;
      using StringLcpPtr = typename tlx::sort_strings_detail::StringLcpPtr<StringSet, size_t>;
      using AllToAllv = dsss::mpi::AllToAllStringImplPrefixDoubling<StringLcpPtr, dsss::mpi::AllToAllvSmall>;
      AllToAllv sender;

      dsss::mpi::environment env;


      auto data = commonSetup<StringSet>();
      StringLcpContainer<StringSet>& container = data.container; 
      StringLcpPtr localStringPtr = container.make_string_lcp_ptr();
      std::vector<size_t>& sendCounts = data.sendCounts;
      size_t commonPrefixLength = data.commonPrefixLength;
      size_t size = data.size;
      size_t sizePerPE = data.sizePerPE;

      // all PEs get same amount of strings
      // distinguishing prefix size is reduced by 1 if one wants to simulate that all locally generated strings are equal
      size_t distinguishingPrefixLength = simulateDuplicates ? commonPrefixLength : commonPrefixLength + 1;
      std::vector<size_t> distinguishingPrefix(size, distinguishingPrefixLength);

      // exchange strings
      auto recvContainer = sender.alltoallv(localStringPtr, sendCounts, distinguishingPrefix);

      
      const auto& recvLcps = recvContainer.lcps();
      const auto& recvRawStrings = recvContainer.raw_strings();

      // get expected raw strings
      std::vector<unsigned char> expectedRawStrings;
      if (simulateDuplicates)
        expectedRawStrings = getExpectedRawStrings_DuplicatesForPrefixDoubling(sizePerPE, commonPrefixLength);
      else
        expectedRawStrings = getExpectedRawStringsForPrefixDoubling(sizePerPE, commonPrefixLength);

      // check whether lcp values are correct, first lcp value from each PE must be 0, the rest equal to commonPrefixLength
      const std::vector<size_t>& expectedLcps = getExpectedLcpValues(sizePerPE, commonPrefixLength);

      tlx_die_unless(recvLcps == expectedLcps);
      tlx_die_unless(recvRawStrings == expectedRawStrings);
    }
  }
}

int main() {
  using dss_schimek::measurement::MeasuringTool;
  MeasuringTool& measuringTool = MeasuringTool::measuringTool();
  dsss::mpi::environment env;

  bool simulateDuplicates = false;
  dss_schimek::tests::AllToAllStringImplPrefixDoubling_test(simulateDuplicates);

  measuringTool.reset(); // TODO add disable feature in measuringTool

  simulateDuplicates = true;
  dss_schimek::tests::AllToAllStringImplPrefixDoubling_test(simulateDuplicates);

  env.finalize();
}
