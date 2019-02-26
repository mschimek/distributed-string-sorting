#pragma once

#include <algorithm>
#include <numeric>
#include <type_traits>
#include <random>

#include "strings/stringptr.hpp"
#include "strings/stringset.hpp"
#include "strings/stringtools.hpp"

//#include "sorter/local/strings/insertion_sort_unified.hpp"
//#include "sorter/local/strings/multikey_quicksort_unified.hpp"
//#include "sorter/local/strings/radix_sort_unified.hpp"
#include "sorter/distributed/bloomfilter.hpp"
#include "encoding/golomb_encoding.hpp"

#include "mpi/alltoall.hpp"
#include "mpi/allgather.hpp"
#include "mpi/synchron.hpp"
#include "mpi/is_sorted.hpp"
#include "mpi/byte_encoder.hpp"

#include "merge/stringtools.hpp"
#include "merge/bingmann-lcp_losertree.hpp"

#include "util/timer.hpp"
#include <tlx/sort/strings/radix_sort.hpp>
#include <tlx/sort/strings/string_ptr.hpp>
#include <tlx/algorithm/multiway_merge.hpp>



namespace dss_schimek {

  struct DuplicateWitness {
    size_t stringIndex;
    size_t PEIndex;
    size_t depth;
    size_t witnessStringIndex;
    size_t witnessPEIndex;
    DuplicateWitness(size_t stringIndex, size_t PEIndex, size_t depth, size_t witnessStringIndex, size_t witnessPEIndex)
      : stringIndex(stringIndex), PEIndex(PEIndex), depth(depth), witnessStringIndex(witnessStringIndex), witnessPEIndex(witnessPEIndex) 
      {}

    friend std::ostream& operator<< (std::ostream& stream, const DuplicateWitness& duplicateWitness) {
      return stream << "[" << duplicateWitness.stringIndex 
                    << ", " << duplicateWitness.PEIndex 
                    << ", " << duplicateWitness.depth
                    << ", " << duplicateWitness.witnessStringIndex 
                    << " " << duplicateWitness.witnessPEIndex 
                    << "]";
    }
  };

  template <typename StringSet>
    class Tracker {
      using String = typename StringSet::String;
      using CharIt = typename StringSet::CharIterator;
      using StringLcpPtr = tlx::sort_strings_detail::StringLcpPtr<StringSet, size_t>;
      using Container = dss_schimek::StringLcpContainer<StringSet>;

      public:

      void printDuplicateWitness(const size_t stringIndex, const size_t PEIndex) {
        for (const DuplicateWitness elem : overallWitnessesPerPE[PEIndex])
          if (elem.stringIndex == stringIndex) {
            StringSet ownStringSet = stringSets[PEIndex];
            String ownString = ownStringSet[ownStringSet.begin() + stringIndex];
            auto ownChars = ownStringSet.get_chars(ownString, 0);
            StringSet witnessStringSet = stringSets[elem.witnessPEIndex];
            String witnessString = witnessStringSet[witnessStringSet.begin() + elem.witnessStringIndex];
            auto witnessChars = witnessStringSet.get_chars(witnessString, 0);
            std::cout << elem << std::endl;
            std::cout << ownChars << " " << hash(ownChars, elem.depth, filterSize) << std::endl; 
            std::cout << witnessChars << " " << hash(witnessChars, elem.depth, filterSize) << std::endl; 
          }
      }

      std::vector<size_t>& getResultsOf(const size_t rank) {
        return resultsPerPE[rank];
      }

      std::vector<size_t>& getCandidatesOf(const size_t rank) {
        return candidatesPerPE[rank];
      }

      void init(StringLcpPtr localStringData) {
        allgatherStrings(localStringData);
        initCandidates();
        initResults();
      }

      void startIteration(const size_t depth) {
        this->depth = depth;
        generateHashTriples();
        detectDuplicates();
        writeResults();
        updateForNextRound();
      }


      Tracker(const size_t filterSize) : env{}, depth(0), filterSize(filterSize), 
        overallWitnessesPerPE(env.size()),
        curItWitnessesPerPE(env.size()),
        resultsPerPE(env.size()),
        candidatesPerPE(env.size()),
        eosCandidatesPerPE(env.size()) {}

      private: 
      Container container;

      dsss::mpi::environment env;
      size_t depth;
      const size_t filterSize;
      std::vector<StringSet> stringSets;
      std::vector<HashTriple> allHashTriples;
      std::vector<std::vector<DuplicateWitness>> overallWitnessesPerPE;
      std::vector<std::vector<DuplicateWitness>> curItWitnessesPerPE;
      std::vector<std::vector<size_t>> resultsPerPE;

      std::vector<std::vector<size_t>> candidatesPerPE;
      std::vector<std::vector<size_t>> eosCandidatesPerPE;


      std::vector<size_t> getDuplicatesOf(const size_t rank) {
        std::vector<size_t> duplicates;
        for (const auto& witness : curItWitnessesPerPE[rank])
          duplicates.push_back(witness.stringIndex);
        return duplicates;
      }

      void initCandidates() {
        for (size_t curRank = 0; curRank < env.size(); ++curRank)
          for(size_t i = 0; i < stringSets[curRank].size(); ++i) 
            candidatesPerPE[curRank].push_back(i);
      }
      void initResults() {
        for (size_t curRank = 0; curRank < env.size(); ++curRank)
          for(size_t i = 0; i < stringSets[curRank].size(); ++i) 
            resultsPerPE[curRank].push_back(0);

      }

      void updateForNextRound() {
        for (size_t curRank = 0; curRank < env.size(); ++curRank) {
          candidatesPerPE[curRank] = getDuplicatesOf(curRank);
          eosCandidatesPerPE[curRank].clear();
          std::copy(curItWitnessesPerPE[curRank].begin(), curItWitnessesPerPE[curRank].end(), std::back_inserter(overallWitnessesPerPE[curRank]));
          curItWitnessesPerPE[curRank].clear();
        }
        allHashTriples.clear();
      }
      inline size_t hash(CharIt str, const size_t maxDepth, const size_t m) {
        size_t hash = 5381;
        size_t c = 0, i = 0;

        while ((c = *str++) && i < maxDepth) {
          hash = ((hash << 5) + hash) + c * 33; /* hash * 33 + c */
          ++i;
        }
        return hash % m;
      }

      void writeResults() {
        for (size_t curRank = 0; curRank < env.size(); ++curRank) {
          StringSet curSet = stringSets[curRank];
          for (size_t curCandidate : candidatesPerPE[curRank])
            resultsPerPE[curRank][curCandidate] = depth; 

          for (size_t curEOSCandidate : eosCandidatesPerPE[curRank]) {
            String str = curSet[curSet.begin() + curEOSCandidate];
            size_t length = curSet.get_length(str);
            resultsPerPE[curRank][curEOSCandidate] = length;
          }
        }
      }

      void detectDuplicates() {
        dsss::mpi::environment env;

        if (allHashTriples.empty())
          return;

        std::stable_sort(allHashTriples.begin(), allHashTriples.end());
        HashTriple prevHashTriple = allHashTriples[0];
        bool duplicate = false;

        for(size_t i = 1; i < allHashTriples.size(); ++i) {
          const HashTriple curHashTriple = allHashTriples[i];
          if (prevHashTriple.hashValue == curHashTriple.hashValue) {
            curItWitnessesPerPE[prevHashTriple.PEIndex].emplace_back(
                prevHashTriple.stringIndex,
                prevHashTriple.PEIndex,
                depth,
                curHashTriple.stringIndex,
                curHashTriple.PEIndex);  
            duplicate = true;
          } else if (duplicate) {
            curItWitnessesPerPE[prevHashTriple.PEIndex].emplace_back(
                prevHashTriple.stringIndex, 
                prevHashTriple.PEIndex, 
                depth,
                allHashTriples[i - 2].stringIndex, 
                allHashTriples[i - 2].PEIndex);  
            duplicate = false;
          }
          prevHashTriple = curHashTriple;
        }
        if (duplicate) {
          const size_t secondToLastIndex = allHashTriples.size() - 2;
          curItWitnessesPerPE[prevHashTriple.PEIndex].emplace_back(
              prevHashTriple.stringIndex, 
              prevHashTriple.PEIndex, 
              depth,
              allHashTriples[secondToLastIndex].stringIndex, 
              allHashTriples[secondToLastIndex].PEIndex);  
        }
      }

      void generateHashTriples() {
        dsss::mpi::environment env;
        for (size_t curRank = 0; curRank < env.size(); ++curRank) {
          StringSet curSet = stringSets[curRank];
          for (size_t i = 0; i < candidatesPerPE[curRank].size(); ++i) {
            size_t curCandidate = candidatesPerPE[curRank][i];
            String curString = curSet[curSet.begin() + curCandidate];
            const size_t length = curSet.get_length(curString);
            if (depth > length) {
              eosCandidatesPerPE[curRank].push_back(curCandidate);
            } else {
              const size_t curHash = hash(curSet.get_chars(curString, 0), depth, filterSize);
              allHashTriples.emplace_back(curHash, curCandidate, curRank);
            }
          } 
        }
      }

      void splitCandidates(const std::vector<size_t>& candidatesSizes, const std::vector<size_t>& candidates) {
        auto curIt = candidates.begin();
        for (size_t curRank = 0; curRank < candidatesSizes.size(); ++curRank) {
          const size_t curCandidatesSize = candidatesSizes[curRank];
          std::copy(curIt, curIt + curCandidatesSize, std::back_inserter(candidatesPerPE[curRank]));
          curIt += curCandidatesSize;
        }
      }

      void splitResults(const std::vector<size_t>& resultsSizes, const std::vector<size_t>& results) {
        auto curIt = results.begin();
        for (size_t curRank = 0; curRank < resultsSizes.size(); ++curRank) {
          const size_t curResultsSize = resultsSizes[curRank];
          std::copy(curIt, curIt + curResultsSize, std::back_inserter(resultsPerPE[curRank]));
          curIt += curResultsSize;
        }
      }

      void splitContainer(const std::vector<size_t>& stringSetSizes) {
        StringSet ss = container.make_string_set();
        std::vector<size_t> stringSetBounds;

        auto begin = ss.begin();
        size_t curStringSetBound = 0;
        for (size_t curRank = 0; curRank < stringSetSizes.size(); ++curRank) {
          StringSet newStringSet(begin + curStringSetBound, begin + curStringSetBound + stringSetSizes[curRank]);
          stringSets.push_back(newStringSet);
          curStringSetBound += stringSetSizes[curRank];
        }
      }

      void allgatherStrings(StringLcpPtr strptr) {
        StringSet ss = strptr.active();
        std::vector<unsigned char> send_buffer;

        for (size_t j = 0; j < ss.size(); ++j) {
          String str = ss[ss.begin() + j];
          size_t string_length = ss.get_length(str) + 1; 
          std::copy_n(ss.get_chars(str, 0), string_length, std::back_inserter(send_buffer));
        }

        size_t numStrings = ss.size();
        std::vector<size_t> recvCounts = dsss::mpi::allgather(numStrings);
        std::vector<unsigned char> recvBuffer = dsss::mpi::allgatherv(send_buffer);
        container.update(std::move(recvBuffer));
        splitContainer(recvCounts);
      }

      void allgatherCandidates(std::vector<size_t>& candidates) {
        size_t candidatesSize = candidates.size();
        std::vector<size_t> recvCounts = dsss::mpi::allgather(candidatesSize);
        std::vector<size_t> recvCandidates = dsss::mpi::allgatherv(candidates);
        splitCandidates(recvCounts, recvCandidates);
      }

      void allgatherResults(std::vector<size_t>& results) {
        size_t resultsSize = results.size();
        std::vector<size_t> recvCounts = dsss::mpi::allgather(resultsSize);
        std::vector<size_t> recvResults = dsss::mpi::allgatherv(results);
        splitResults(recvCounts, recvResults);
      }
    };
}
