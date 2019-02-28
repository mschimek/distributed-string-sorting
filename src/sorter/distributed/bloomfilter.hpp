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


  struct Duplicate {
    size_t index;
    bool hasReachedEOS;
    Duplicate(size_t index, bool hasReachedEOS) : index(index), hasReachedEOS(hasReachedEOS) {}
  };

  struct HashTriple {
    size_t hashValue;
    size_t stringIndex; 
    size_t PEIndex;
    HashTriple() = default;

    HashTriple(size_t hashValue, size_t stringIndex, size_t PEIndex) : hashValue(hashValue), stringIndex(stringIndex), PEIndex(PEIndex) {}

    bool operator<(const HashTriple& rhs) const {
      return hashValue < rhs.hashValue;
    }

    friend std::ostream& operator<< (std::ostream& stream, const HashTriple& hashTriple) {
      return stream << "[" << hashTriple.hashValue << ", " << hashTriple.stringIndex << ", "  << hashTriple.PEIndex << "]";
    }
    operator std::string() const { return "(" + std::to_string(hashValue) + ", " + std::to_string(stringIndex) + ", " +  std::to_string(PEIndex) + ")"; }
  };
  
  struct StringTriple {
    const unsigned char* string;
    size_t stringIndex; 
    size_t PEIndex;
    StringTriple() = default;

    StringTriple(const unsigned char* string_, size_t stringIndex, size_t PEIndex) : string(string_), stringIndex(stringIndex), PEIndex(PEIndex) {}

    bool operator<(const StringTriple& rhs) const {
      size_t i = 0;
      while (string[i] != 0 && string[i] == rhs.string[i])
        ++i;
      return string[i] < rhs.string[i];
    }

    friend std::ostream& operator<< (std::ostream& stream, const StringTriple& stringTriple) {
      return stream << "[" << stringTriple.string << ", " << stringTriple.stringIndex << ", " << stringTriple.PEIndex << "]" << std::endl;
    }
  };

  struct HashTripleComp {
    bool operator() (const HashTriple& lhs, const HashTriple& rhs) { 
      const bool firstCrit = lhs.hashValue < rhs.hashValue;
      const bool secondCrit = !firstCrit && lhs.stringIndex < rhs.stringIndex;
      return firstCrit || secondCrit;
    }
    bool operator() (const size_t& lhs, const HashTriple& rhs) { return lhs < rhs.hashValue; }
  };

  struct HashStringIndex {
    size_t hashValue;
    size_t stringIndex;
    bool isLocalDuplicate = false;
    bool isLocalDuplicateButSendAnyway = false;
    HashStringIndex(const size_t hashValue, const size_t stringIndex, bool isLocalDuplicate, bool isLocalDuplicateButSendAnyway)
      : hashValue(hashValue), stringIndex(stringIndex), isLocalDuplicate(isLocalDuplicate), isLocalDuplicateButSendAnyway(isLocalDuplicateButSendAnyway) {}
    HashStringIndex(const size_t hashValue, const size_t stringIndex) : hashValue(hashValue), stringIndex(stringIndex), isLocalDuplicate(false) {}

    bool operator< (const HashStringIndex& rhs) const {
      return hashValue < rhs.hashValue;
    } 

    friend std::ostream& operator<< (std::ostream& stream, const HashStringIndex& hashStringIndex) {
      return stream << "[" << hashStringIndex.hashValue << ", " << hashStringIndex.stringIndex << ", localDup: " << hashStringIndex.isLocalDuplicate << 
                              ", sendAnyway: " << hashStringIndex.isLocalDuplicateButSendAnyway << "]";
    } 
  };

  struct HashPEIndex {
    size_t hashValue;
    size_t PEIndex;
    HashPEIndex() : hashValue(0), PEIndex(0) {}
    HashPEIndex(const size_t hashValue, const size_t PEIndex) : hashValue(hashValue), PEIndex(PEIndex) {}

    bool operator< (const HashPEIndex& rhs) const {
      return hashValue < rhs.hashValue;
    } 

    friend std::ostream& operator<< (std::ostream& stream, const HashPEIndex& hashPEIndex) {
      return stream << "[" << hashPEIndex.hashValue << ", " << hashPEIndex.PEIndex << "]";
    } 
  };

  struct AllToAllHashesNaive {
    template <typename DataType>
      static inline std::vector<DataType> allToAllv(const std::vector<DataType>& sendData, const std::vector<size_t>& intervalSizes) {
        return dsss::mpi::AllToAllvSmall::alltoallv(sendData.data, intervalSizes);
      }
  };

  std::vector<size_t> computeIntervalSizes(const std::vector<size_t>& hashes, const size_t bloomFilterSize,
      dsss::mpi::environment env = dsss::mpi::environment()) {
    std::vector<size_t> indices;
    indices.reserve(env.size());
    auto curPosWithinVector = hashes.begin(); 
    size_t upperPartitionLimit = bloomFilterSize / env.size();
    for (size_t i = 0; i < env.size(); ++i) {
      const size_t upperPartitionLimit = (i + 1) * (bloomFilterSize / env.size()) - 1;
      auto pos = std::upper_bound(curPosWithinVector, hashes.end(), upperPartitionLimit);
      indices.push_back(pos - curPosWithinVector);
      curPosWithinVector = pos;
    }
    return indices;
  }

  struct RecvData {
    std::vector<size_t> data;
    std::vector<size_t> intervalSizes;
    std::vector<size_t> globalOffsets;
    RecvData(std::vector<size_t>&& data, std::vector<size_t>&& intervalSizes, std::vector<size_t>&& globalOffsets)
      : data(std::move(data)), intervalSizes(std::move(intervalSizes)), globalOffsets(std::move(globalOffsets)) 
    {}
  };

  template<typename SendPolicy>
  struct SendOnlyHashesToFilter : private SendPolicy {
    using SendType = size_t;

    static inline std::vector<SendType> extractSendValues(const std::vector<HashStringIndex>& hashStringIndices) {
    dsss::mpi::environment env;
      std::vector<size_t> hashValues;
      hashValues.reserve(hashStringIndices.size());
      for (const auto& hashStringIndex : hashStringIndices) {
         hashValues.push_back(hashStringIndex.hashValue);
      }

      //dss_schimek::mpi::execute_in_order([&](){
      //    std::cout << "extract send values rank: " << env.rank() << std::endl;
      //  for (const auto& elem: hashValues) {
      //    std::cout << elem << std::endl;
      //  }
      //    });
      return hashValues;
    }

    static inline RecvData sendToFilter(const std::vector<HashStringIndex>& hashes, size_t bloomfilterSize) {
      std::vector<size_t> sendValues = extractSendValues(hashes);
      std::vector<size_t> intervalSizes = computeIntervalSizes(sendValues, bloomfilterSize);
      std::vector<size_t> offsets;
      offsets.reserve(intervalSizes.size());
      offsets.push_back(0);
      std::partial_sum(intervalSizes.begin(), intervalSizes.end() - 1, std::back_inserter(offsets));
      dsss::mpi::environment env;
      offsets = dsss::mpi::alltoall(offsets);

      std::vector<size_t> recvIntervalSizes = dsss::mpi::alltoall(intervalSizes);
      std::vector<size_t> result = SendPolicy::alltoallv(sendValues.data(), intervalSizes);
      return RecvData(std::move(result), std::move(recvIntervalSizes), std::move(offsets));
    }

    static inline std::vector<HashPEIndex> addPEIndex(const RecvData& recvData) {
      std::vector<HashPEIndex> hashesPEIndex;
      hashesPEIndex.reserve(recvData.data.size());

      size_t curPE = 0;
      size_t curBoundary = recvData.intervalSizes[0];
      for (size_t i = 0; i < recvData.data.size(); ++i) {
        while (i == curBoundary) 
          curBoundary += recvData.intervalSizes[++curPE];
        hashesPEIndex.emplace_back(recvData.data[i], curPE);
      }  
      return hashesPEIndex;
    }
  };

  std::vector<size_t> getEncoding(const std::vector<size_t>& values, const std::vector<size_t>& intervalSizes) {
      std::vector<size_t> encodedValues; 
      encodedValues.reserve(values.size()); // should i do this? 
      auto outIt = std::back_inserter(encodedValues);
      size_t offset = 0;
      auto inputIterator = values.begin();
      for (const size_t curIntervalSize : intervalSizes) {
        getDeltaEncoding(inputIterator, inputIterator + curIntervalSize, outIt);
        inputIterator += curIntervalSize;
      }
      return encodedValues;
    }

  struct FindDuplicates {
    using DataType = HashPEIndex;

    static inline std::vector<size_t> findDuplicates(std::vector<HashPEIndex>& hashPEIndices, const RecvData& recvData, Timer& timer, size_t curIteration) {
      using ConstIterator = std::vector<HashPEIndex>::const_iterator;
      using Iterator = std::vector<HashPEIndex>::iterator;
      using IteratorPair = std::pair<Iterator, Iterator>;
      dsss::mpi::environment env;

      timer.add(std::string("bloomfilter_recvHashValues"), curIteration, hashPEIndices.size());
      env.barrier();
      timer.start(std::string("bloomfilter_findDuplicatesOverallIntern"), curIteration);
      timer.localStart(std::string("bloomfilter_findDuplicatesSetup"), curIteration);
      std::vector<IteratorPair> iteratorPairs;
      size_t elementsToMerge = std::accumulate(recvData.intervalSizes.begin(), recvData.intervalSizes.end(), 0);
      std::vector<HashPEIndex> mergedElements(elementsToMerge);
      auto outputIt = std::back_inserter(mergedElements);
      Iterator it = hashPEIndices.begin(); 

      for (size_t i = 0; i < recvData.intervalSizes.size(); ++i) {
       iteratorPairs.emplace_back(it, it + recvData.intervalSizes[i]);
       it += recvData.intervalSizes[i];
      }
      timer.end(std::string("bloomfilter_findDuplicatesSetup"), curIteration);
      //
      //++++++++++++++++++++++++++++++++++++++
      //
      timer.localStart(std::string("bloomfilter_findDuplicatesMerge"), curIteration);
      tlx::multiway_merge(iteratorPairs.begin(), iteratorPairs.end(), mergedElements.begin(), elementsToMerge);
      timer.end(std::string("bloomfilter_findDuplicatesMerge"), curIteration);
      //
      //++++++++++++++++++++++++++++++++++++++
      //

      timer.localStart(std::string("bloomfilter_findDuplicatesFind"), curIteration);
      std::vector<std::vector<size_t>> result_sets(recvData.intervalSizes.size());
      std::vector<size_t> counters(recvData.intervalSizes.size(), 0);

      HashPEIndex prevHashTriple = mergedElements.empty() ? HashPEIndex{0, 0} : mergedElements[0];
      bool duplicate = false;


      for(size_t i = 1; i < mergedElements.size(); ++i) {
        const HashPEIndex curHashTriple = mergedElements[i];
        if (prevHashTriple.hashValue == curHashTriple.hashValue) {
          result_sets[prevHashTriple.PEIndex].push_back(counters[prevHashTriple.PEIndex]++);
          duplicate = true;
        } else if (duplicate) {
          result_sets[prevHashTriple.PEIndex].push_back(counters[prevHashTriple.PEIndex]++);
          duplicate = false;
        } else {
          ++counters[prevHashTriple.PEIndex];
        }
        prevHashTriple = curHashTriple;
      }

      if (duplicate)
          result_sets[prevHashTriple.PEIndex].push_back(counters[prevHashTriple.PEIndex]++);

      std::vector<size_t> sendBuffer;
      sendBuffer.reserve(elementsToMerge);
      std::vector<size_t> sendCounts_;
        // write back to PEs
        for (size_t i = 0; i < result_sets.size(); ++i) {
          sendCounts_.push_back(result_sets[i].size());
          for (size_t j = 0; j < result_sets[i].size(); ++j) {
            sendBuffer.push_back(result_sets[i][j] + recvData.globalOffsets[i] );
          }
        }
      timer.end(std::string("bloomfilter_findDuplicatesFind"), curIteration);
      //
      //++++++++++++++++++++++++++++
      //
      size_t totalNumSendDuplicates = std::accumulate(sendCounts_.begin(), sendCounts_.end(), 0);
      timer.add(std::string("bloomfilter_findDuplicatesNumberSendDups"), curIteration, totalNumSendDuplicates);
      timer.start(std::string("bloomfilter_findDuplicatesSendDups"), curIteration);
      int mpiSmallTypes = 0;
      if (totalNumSendDuplicates > 0)
        mpiSmallTypes = 1;
      bool dupsToSend = (0 != dsss::mpi::allreduce_max(mpiSmallTypes));
      std::vector<size_t> duplicates;
      if (dupsToSend)
       duplicates = dsss::mpi::AllToAllvSmall::alltoallv(sendBuffer.data(), sendCounts_);
      timer.end(std::string("bloomfilter_findDuplicatesSendDups"), curIteration);
      timer.end(std::string("bloomfilter_findDuplicatesOverallIntern"), curIteration);

      return duplicates; 
    }
   
    
    void setUniqueElements(std::vector<bool>& duplicates, std::vector<size_t>& depth, const size_t curDepth, const std::vector<HashStringIndex>& originalMapping) {
      for (size_t i = 0; i < duplicates.size(); ++i)
        if (!duplicates[i]) {
          const size_t stringIndex = originalMapping[i].stringIndex;
          depth[stringIndex] = curDepth;
        }
    }
//TODO add && reference for localDuplicates
    std::vector<size_t> getIndicesOfDuplicates(std::vector<size_t>& localDuplicates, std::vector<size_t>& remoteDuplicates, const std::vector<HashStringIndex>& originalMapping) {
      std::vector<size_t> indicesOfAllDuplicates(localDuplicates);
      dsss::mpi::environment env;
      //indicesOfAllDuplicates.reserve(localDuplicates.size() + remoteDuplicates.size());
      dss_schimek::mpi::execute_in_order([&](){
          std::cout << "rank: " << env.rank() << std::endl;
          for (size_t i = 0; i < remoteDuplicates.size(); ++i) {
          const size_t curIndex = remoteDuplicates[i];
          bool isAlsoLocalDuplicate = originalMapping[curIndex].isLocalDuplicateButSendAnyway;
          std::cout << i << " curIndex " << curIndex << " isLocalDuplicateButSendAnyway " << isAlsoLocalDuplicate << std::endl;
          if (!isAlsoLocalDuplicate) {
          const size_t stringIndex = originalMapping[curIndex].stringIndex;
          indicesOfAllDuplicates.push_back(stringIndex);
          }
          }  
          });

      return indicesOfAllDuplicates;
    }
  };

  struct SendHashesAndPEIndexToFilter {
    using SendType = HashStringIndex;
  };

  class AllToAllHashValuesNaive {
    public:
      static inline std::vector<HashTriple> allToAllHashTriples(std::vector<HashTriple> hashTriples, std::vector<size_t> intervalSizes) {
        return dsss::mpi::AllToAllvSmall::alltoallv(hashTriples.data(), intervalSizes);
      }
  };

  class AllToAllHashValuesPipeline {
    public:
      static inline std::vector<HashTriple> allToAllHashTriples(std::vector<HashTriple> hashTriples, std::vector<size_t> intervalSizes) {
        dsss::mpi::environment env;
        const bool PESizeIsEven = env.size() % 2 == 0;
        const size_t modulator = PESizeIsEven ? env.size() - 1 : env.size();
        for (size_t j = 0; j < env.size(); ++j) {
          size_t idlePE = (env.size() / 2 * j) % (env.size() - 1);
          if (PESizeIsEven && env.rank() == env.size() - 1) {
            //exchange with PE idle
          } else {
            if (PESizeIsEven && env.rank() == idlePE) {
              //exchange with PE env.size() - 1
            } else {
              if (PESizeIsEven) {
                // exchange with PE   ((j - i) % env.size() - 1
              } else {
                // exchange with PE (( j - i) % env.size() 
              }
            }
          }

        }
        return std::vector<HashTriple>();
      }
  };


  template <typename StringSet, typename AllToAllHashValuePolicy, typename FindDuplicatesPolicy, template<typename> typename SendPolicy>
    class BloomFilter : private AllToAllHashValuePolicy,
                        private FindDuplicatesPolicy, 
                        private SendPolicy<dsss::mpi::AllToAllvSmall> {

      using String = typename StringSet::String;
      using Iterator = typename StringSet::Iterator;
      using CharIt = typename StringSet::CharIterator;
      using StringLcpPtr = typename tlx::sort_strings_detail::StringLcpPtr<StringSet, size_t>;

      HashTripleComp comp;
      dsss::mpi::environment env;
      public:

      const size_t bloomFilterSize = std::numeric_limits<uint32_t>::max();

      std::vector<size_t> computeIntervals(std::vector<HashTriple>& hashes) {
        std::vector<size_t> indices;
        indices.reserve(env.size());
        typename std::vector<HashTriple>::iterator curPosWithinVector = hashes.begin(); 
        size_t upperPartitionLimit = bloomFilterSize / env.size();
        for (size_t i = 0; i < env.size(); ++i) {
          size_t upperPartitionLimit = (i + 1) * (bloomFilterSize / env.size()) - 1;
          auto pos = std::upper_bound(curPosWithinVector, hashes.end(), upperPartitionLimit, comp);
          indices.emplace_back(pos - curPosWithinVector);
          curPosWithinVector = pos;
        }
        return indices;
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

      std::vector<size_t> getOwnDuplicates(std::vector<HashTriple>& hashTriples) {
        dsss::mpi::environment env;
        std::vector<size_t> duplicateIndices;

        if (hashTriples.empty())
          return duplicateIndices;

        std::stable_sort(hashTriples.begin(), hashTriples.end());
        HashTriple prevHashTriple = hashTriples[0];
        bool duplicate = false;

        for(size_t i = 1; i < hashTriples.size(); ++i) {
          const HashTriple curHashTriple = hashTriples[i];
          if (prevHashTriple.hashValue == curHashTriple.hashValue) {
            if (prevHashTriple.PEIndex == env.rank())
              duplicateIndices.push_back(prevHashTriple.stringIndex);
            duplicate = true;
          } else if (duplicate) {
            if (prevHashTriple.PEIndex == env.rank())
              duplicateIndices.push_back(prevHashTriple.stringIndex);
            duplicate = false;
          }
          prevHashTriple = curHashTriple;
        }
        if (duplicate)
          if (prevHashTriple.PEIndex == env.rank())
            duplicateIndices.push_back(prevHashTriple.stringIndex);

        return duplicateIndices;
      }

      void computeExactDistPrefixLengths(std::vector<StringTriple>& stringTriples, std::vector<size_t>& distinguishingPrefixLength) {
        dsss::mpi::environment env;
        if (stringTriples.empty())
          return;

        std::stable_sort(stringTriples.begin(), stringTriples.end());
        StringTriple prevHashTriple = stringTriples[0];
        bool duplicate = false;

        for (size_t i = 1; i < stringTriples.size(); ++i) {
          const StringTriple prevStringTriple = stringTriples[i - 1];
          const StringTriple curStringTriple = stringTriples[i];

          const unsigned char * s1 = prevStringTriple.string;
          const unsigned char * s2 = curStringTriple.string;

          const size_t distValuePrevCur = 1 + dss_schimek::calc_lcp(s1, s2); 
          if (prevStringTriple.PEIndex == env.rank()) {
            const size_t oldValue = distinguishingPrefixLength[prevStringTriple.stringIndex];
            distinguishingPrefixLength[prevStringTriple.stringIndex] = distValuePrevCur > oldValue ? distValuePrevCur : oldValue;
          }
          if (curStringTriple.PEIndex == env.rank()) {
            const size_t oldValue = distinguishingPrefixLength[curStringTriple.stringIndex];
            distinguishingPrefixLength[curStringTriple.stringIndex] = distValuePrevCur > oldValue ? distValuePrevCur : oldValue;
          }
        }
      }

      struct ContainerSizesIndices {
        using Container = dss_schimek::StringLcpContainer<StringSet>;
        Container container;
        std::vector<size_t> intervalSizes;
        std::vector<size_t> stringIndices;
        ContainerSizesIndices(Container&& container, std::vector<size_t>&& intervalSizes, std::vector<size_t>&& stringIndices) 
          : container(std::move(container)), intervalSizes(std::move(intervalSizes)), stringIndices(std::move(stringIndices)) {}
      };

      
      
      std::vector<StringTriple> generateStringTriples(ContainerSizesIndices& containerSizesIndices, const size_t depth) {
        const std::vector<size_t>& intervalSizes = containerSizesIndices.intervalSizes;
        const StringSet globalSet = containerSizesIndices.container.make_string_set();
        std::vector<size_t>& stringIndices = containerSizesIndices.stringIndices;

        size_t totalNumSentStrings = std::accumulate(intervalSizes.begin(), intervalSizes.end(), 0);

        std::vector<StringTriple> stringTriples;
        if (totalNumSentStrings == 0) 
          return stringTriples;

        stringTriples.reserve(totalNumSentStrings);
        size_t curOffset = 0; 
        auto begin = globalSet.begin();

        for (size_t curRank = 0; curRank < env.size(); ++curRank) {
          for (size_t i = 0; i < intervalSizes[curRank]; ++i) {
            String curString = globalSet[begin + curOffset + i];
            stringTriples.emplace_back(globalSet.get_chars(curString, 0), stringIndices[curOffset + i], curRank);
          } 
          curOffset += intervalSizes[curRank];
        }
        return stringTriples;
      }

      ContainerSizesIndices allgatherStrings(StringLcpPtr strptr,
          std::vector<size_t>& candidates) {

        StringSet ss = strptr.active();
        std::vector<unsigned char> send_buffer;

        for (size_t j = 0; j < candidates.size(); ++j) {
          String str = ss[ss.begin() + candidates[j]];
          size_t string_length = ss.get_length(str) + 1; 
          std::copy_n(ss.get_chars(str, 0), string_length, std::back_inserter(send_buffer));
        }
        size_t numStrings = candidates.size();

        std::vector<size_t> recvCounts = dsss::mpi::allgather(numStrings);
        std::vector<size_t> stringIndices = dsss::mpi::allgatherv(candidates);
        std::vector<unsigned char> recvBuffer = dsss::mpi::allgatherv(send_buffer);
        return ContainerSizesIndices(std::move(recvBuffer), std::move(recvCounts), std::move(stringIndices));
      }

      template <typename T>
        struct GeneratedHashStructuresEOSCandidates {
          std::vector<T> data;
          std::vector<size_t> eosCandidates;
          GeneratedHashStructuresEOSCandidates(std::vector<T>&& data, std::vector<size_t>&& eosCandidates) 
            : data(std::move(data)), eosCandidates(std::move(eosCandidates)) {}
        };


      GeneratedHashStructuresEOSCandidates<HashTriple> generateHashTriples(ContainerSizesIndices& containerSizesIndices, const size_t depth) {
        const std::vector<size_t>& intervalSizes = containerSizesIndices.intervalSizes;
        const StringSet globalSet = containerSizesIndices.container.make_string_set();
        std::vector<size_t>& stringIndices = containerSizesIndices.stringIndices;
        std::vector<size_t> eosCandidates;

        size_t totalNumSentStrings = std::accumulate(intervalSizes.begin(), intervalSizes.end(), 0);

        std::vector<HashTriple> hashTriples;
        if (totalNumSentStrings == 0) 
          return GeneratedHashStructuresEOSCandidates<HashTriple>(std::move(hashTriples), std::move(eosCandidates));

        hashTriples.reserve(totalNumSentStrings);
        size_t curOffset = 0; 
        auto begin = globalSet.begin();

        for (size_t curRank = 0; curRank < env.size(); ++curRank) {
          for (size_t i = 0; i < intervalSizes[curRank]; ++i) {
            String curString = globalSet[begin + curOffset + i];
            const size_t length = globalSet.get_length(curString);
            if (depth > length) {
              if (curRank == env.rank())
                eosCandidates.push_back(stringIndices[curOffset + i]);
            } else {
              const size_t curHash = hash(globalSet.get_chars(curString, 0), depth, bloomFilterSize);
              hashTriples.emplace_back(curHash, stringIndices[curOffset + i], curRank);
            }
          } 
          curOffset += intervalSizes[curRank];
        }
        return GeneratedHashStructuresEOSCandidates<HashTriple>(std::move(hashTriples), std::move(eosCandidates));
      }

      GeneratedHashStructuresEOSCandidates<HashStringIndex> generateHashStringIndices(StringSet ss, const std::vector<size_t>& candidates, const size_t depth) {
        std::vector<HashStringIndex> hashStringIndices;
        std::vector<size_t> eosCandidates;
        hashStringIndices.reserve(candidates.size());
        const Iterator begin = ss.begin();

        for (const size_t curCandidate : candidates) {
          String curString = ss[begin + curCandidate];
          const size_t length = ss.get_length(curString);
          if (depth > length) {
            eosCandidates.push_back(curCandidate);
          } else {
            const size_t curHash = hash(ss.get_chars(curString, 0), depth, bloomFilterSize);
            hashStringIndices.emplace_back(curHash, curCandidate);
          }
        }
        return GeneratedHashStructuresEOSCandidates(std::move(hashStringIndices), std::move(eosCandidates));
      }

      std::vector<size_t> getIndicesOfLocalDuplicates(std::vector<HashStringIndex>& hashStringIndices) {
        std::vector<size_t> indicesOfLocalDuplicates;
        if (hashStringIndices.empty())
          return indicesOfLocalDuplicates;

        for (size_t i = 0; i < hashStringIndices.size() - 1;) {
          HashStringIndex& pivotHashStringIndex = hashStringIndices[i];
          size_t j = i + 1;
          HashStringIndex& curHashStringIndex = hashStringIndices[j];
          if (curHashStringIndex.hashValue == pivotHashStringIndex.hashValue) {
            indicesOfLocalDuplicates.push_back(pivotHashStringIndex.stringIndex);
            indicesOfLocalDuplicates.push_back(curHashStringIndex.stringIndex);

            pivotHashStringIndex.isLocalDuplicate = true;
            pivotHashStringIndex.isLocalDuplicateButSendAnyway = true;
            curHashStringIndex.isLocalDuplicate = true;
            ++j;
            while (j < hashStringIndices.size() && hashStringIndices[j].hashValue == pivotHashStringIndex.hashValue) {
              hashStringIndices[j].isLocalDuplicate = true;
              indicesOfLocalDuplicates.push_back(hashStringIndices[j].stringIndex);
              ++j;
            }
          }
          i = j;
        }
        return indicesOfLocalDuplicates;
      }
      
      void setDepth(StringLcpPtr strptr, const size_t depth, const std::vector<size_t>& candidates, 
                    const std::vector<size_t>& eosCandidates, std::vector<size_t>& results) {

        // eosCandidates is subset of candidates whose length is <= depth
        StringSet ss = strptr.active();
        for (const size_t curCandidate : candidates)
            results[curCandidate] = depth;

        for (const size_t curEOSCandidate : eosCandidates) {
          String str = ss[ss.begin() + curEOSCandidate];
          size_t length = ss.get_length(str);
            results[curEOSCandidate] = length;
        }
      }

      void filter_exact(StringLcpPtr strptr, const size_t depth, 
          std::vector<size_t>& candidates, std::vector<size_t>& results) {
        ContainerSizesIndices containerSizesIndices = allgatherStrings(strptr, candidates);
        std::vector<StringTriple> globalStringTriples = generateStringTriples(containerSizesIndices, depth);
        computeExactDistPrefixLengths(globalStringTriples, results);
      }

      

      std::vector<size_t> filter_simple(StringLcpPtr strptr, const size_t depth, 
          std::vector<size_t>& candidates, std::vector<size_t>& results) {

        ContainerSizesIndices containerSizesIndices = allgatherStrings(strptr, candidates);
        GeneratedHashStructuresEOSCandidates<HashTriple> generatedData = generateHashTriples(containerSizesIndices, depth);
        std::vector<HashTriple>& globalHashTriples = generatedData.data;
        const std::vector<size_t>& eosCandidates = generatedData.eosCandidates;
        std::vector<size_t> ownDuplicateIndices = getOwnDuplicates(globalHashTriples);
        setDepth(strptr, depth, candidates, eosCandidates, results);

        return ownDuplicateIndices;
      } 

      std::vector<size_t> filter(StringLcpPtr strptr, const size_t depth, const std::vector<size_t>& candidates, std::vector<size_t>& results, Timer& timer, size_t curIteration) {
        dsss::mpi::environment env;

        timer.start(std::string("bloomfilter_generateHashStringIndices"), curIteration);
        GeneratedHashStructuresEOSCandidates<HashStringIndex> hashStringIndicesEOSCandidates = generateHashStringIndices(strptr.active(), candidates, depth);

        std::vector<HashStringIndex>& hashStringIndices = hashStringIndicesEOSCandidates.data;
        const std::vector<size_t>& eosCandidates = hashStringIndicesEOSCandidates.eosCandidates;
        timer.end(std::string("bloomfilter_generateHashStringIndices"), curIteration);

        timer.start(std::string("bloomfilter_sortHashStringIndices"), curIteration);
        std::sort(hashStringIndices.begin(), hashStringIndices.end());
        timer.end(std::string("bloomfilter_sortHashStringIndices"), curIteration);
        

        std::vector<size_t> indicesOfLocalDuplicates = getIndicesOfLocalDuplicates(hashStringIndices);
        std::vector<HashStringIndex> reducedHashStringIndices;
        reducedHashStringIndices.reserve(hashStringIndices.size());
        std::copy_if(hashStringIndices.begin(), hashStringIndices.end(), std::back_inserter(reducedHashStringIndices), [&](const HashStringIndex& v) {
            return !v.isLocalDuplicate || v.isLocalDuplicateButSendAnyway;
            }); dss_schimek::mpi::execute_in_order([&]() { std::cout << "local dups: rank: " << env.rank() << std::endl;
            for (const auto& elem : indicesOfLocalDuplicates)
              std::cout << elem << std::endl;
            });
dss_schimek::mpi::execute_in_order([&]() {
            std::cout << "hashStringIndices dups: rank: " << env.rank() << std::endl;
            for (size_t i = 0; i < hashStringIndices.size() && i << reducedHashStringIndices.size(); ++i)
              std::cout << hashStringIndices[i] << " " << reducedHashStringIndices[i] << std::endl;
            for (volatile size_t i = 0; i < 1000000; ++i) {

            }
            });



        timer.start(std::string("bloomfilter_sendHashStringIndices"), curIteration);
        RecvData recvData = SendPolicy<dsss::mpi::AllToAllvSmall>::sendToFilter(reducedHashStringIndices, bloomFilterSize);
        timer.end(std::string("bloomfilter_sendHashStringIndices"), curIteration);

        timer.start(std::string("bloomfilter_addPEIndex"), curIteration);
        std::vector<HashPEIndex> recvHashPEIndices = SendPolicy<dsss::mpi::AllToAllvSmall>::addPEIndex(recvData);
        timer.end(std::string("bloomfilter_addPEIndex"), curIteration);
        
        //timer.start(std::string("bloomfilter_findDuplicatesOverall"), curIteration);
        std::vector<size_t> indicesOfRemoteDuplicates = FindDuplicatesPolicy::findDuplicates(recvHashPEIndices, recvData, timer, curIteration);
        //timer.end(std::string("bloomfilter_findDuplicatesOverall"), curIteration);
          dss_schimek::mpi::execute_in_order([&]() {
            std::cout << "remote dups: rank: " << env.rank() << std::endl;
            for (const auto& elem : indicesOfRemoteDuplicates)
              std::cout << elem << std::endl;
            });


        timer.start(std::string("bloomfilter_getIndices"), curIteration);
        std::vector<size_t> indicesOfAllDuplicates = FindDuplicatesPolicy::getIndicesOfDuplicates(indicesOfLocalDuplicates, indicesOfRemoteDuplicates, reducedHashStringIndices);
      dss_schimek::mpi::execute_in_order([&]() {
            std::cout << "all dups: rank: " << env.rank() << std::endl;
            for (const auto& elem : indicesOfAllDuplicates)
              std::cout << elem << std::endl;
            });


        std::cout << "hallo " << std::endl;
        timer.end(std::string("bloomfilter_getIndices"), curIteration);

        timer.start(std::string("bloomfilter_setDepth"), curIteration);
        setDepth(strptr, depth, candidates, eosCandidates, results);
        timer.end(std::string("bloomfilter_setDepth"), curIteration);

        return indicesOfAllDuplicates;
      }

      struct RemoveEOSIndices {};
      std::vector<size_t> filter(StringLcpPtr strptr, const size_t depth, const std::vector<size_t>& candidates, std::vector<size_t>& results, RemoveEOSIndices) {
        std::vector<size_t> indicesOfDuplicates = filter(strptr, depth, candidates, results);
        StringSet ss = strptr.active();
        std::vector<size_t> indicesOfDuplicatesWithoutEOS;
        indicesOfDuplicatesWithoutEOS.reserve(indicesOfDuplicates.size());
        for (const size_t curIndex : indicesOfDuplicatesWithoutEOS) {
          String str = ss[ss.begin() + curIndex];
          size_t length = ss.get_length(str) + 1;
          if (depth > length)
            indicesOfDuplicatesWithoutEOS.push_back(curIndex);
        }
        return indicesOfDuplicatesWithoutEOS;
      }

                        };

}
