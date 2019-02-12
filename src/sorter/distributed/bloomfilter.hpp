#pragma once

#include <algorithm>
#include <type_traits>
#include <random>

#include "strings/stringptr.hpp"
#include "strings/stringset.hpp"
#include "strings/stringtools.hpp"

//#include "sorter/local/strings/insertion_sort_unified.hpp"
//#include "sorter/local/strings/multikey_quicksort_unified.hpp"
//#include "sorter/local/strings/radix_sort_unified.hpp"

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
  struct HashStringIndex {
    size_t hashValue;
    size_t stringIndex;
    HashStringIndex(const size_t hashValue, const size_t stringIndex) : hashValue(hashValue), stringIndex(stringIndex) {}

    bool operator< (const HashStringIndex& rhs) const {
      return hashValue < rhs.hashValue;
    } 

    friend std::ostream& operator<< (std::ostream& stream, const HashStringIndex& hashStringIndex) {
      return stream << "[" << hashStringIndex.hashValue << ", " << hashStringIndex.stringIndex << "]";
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


  template<typename SendPolicy>
  struct SendOnlyHashesToFilter : private SendPolicy {
    using SendType = size_t;

    static inline std::vector<SendType> extractSendValues(const std::vector<HashStringIndex>& hashStringIndices) {
      std::vector<size_t> hashValues;
      hashValues.reserve(hashStringIndices.size());
      for (const auto& hashStringIndex : hashStringIndices)
        hashValues.push_back(hashStringIndex.hashValue);
      return hashValues;
    }

    static inline std::pair<std::vector<SendType>, std::vector<size_t>> sendToFilter(const std::vector<HashStringIndex>& hashes, size_t bloomfilterSize) {
      std::vector<size_t> sendValues = extractSendValues(hashes);//TODO return computeIntervalSizes
      std::vector<size_t> intervalSizes = computeIntervalSizes(sendValues, bloomfilterSize);
      std::vector<size_t> recvIntervalSizes = dsss::mpi::alltoall(intervalSizes);
      std::vector<size_t> result = SendPolicy::alltoallv(sendValues.data(), intervalSizes);
      auto p = std::make_pair(result, recvIntervalSizes);
      return p;
    }

    static inline std::vector<HashPEIndex> addPEIndex(const std::vector<size_t>& hashes, const std::vector<size_t>& intervalSizes) {
      std::vector<HashPEIndex> hashesPEIndex;
      hashesPEIndex.reserve(hashes.size());

      //dsss::mpi::environment env;
      //dss_schimek::mpi::execute_in_order([&]() {
      //    std::cout << "intervalSizes  within addPEIndex rank: " << env.rank()  << std::endl;
      //    for (size_t i = 0; i < intervalSizes.size(); ++i)
      //    std::cout << i << " " << intervalSizes[i] << std::endl;
      //    std::cout << "intervalSizes within addPEIndex end" << std::endl;
      //    });

      size_t curPE = 0;
      size_t curBoundary = intervalSizes[0];
      for (size_t i = 0; i < hashes.size(); ++i) {
        while (i == curBoundary) 
          curBoundary += intervalSizes[++curPE];
        hashesPEIndex.emplace_back(hashes[i], curPE);
      }  
      return hashesPEIndex;
    }
  };



  struct FindDuplicates {
    using DataType = HashPEIndex;
    
    static inline std::vector<size_t> findDuplicates(std::vector<HashPEIndex>& hashPEIndices, std::vector<size_t>& intervalSizes) {
      using ConstIterator = std::vector<HashPEIndex>::const_iterator;
      using Iterator = std::vector<HashPEIndex>::iterator;
      using IteratorPair = std::pair<Iterator, Iterator>;
      dsss::mpi::environment env;
      std::vector<IteratorPair> iteratorPairs;

      size_t elementsToMerge = std::accumulate(intervalSizes.begin(), intervalSizes.end(), 0);
      std::vector<HashPEIndex> mergedElements(elementsToMerge);
      auto outputIt = std::back_inserter(mergedElements);
      Iterator it = hashPEIndices.begin(); 
      

      for (size_t i = 0; i < intervalSizes.size(); ++i) {
       iteratorPairs.emplace_back(it, it + intervalSizes[i]);
       it += intervalSizes[i];
      }

      //dss_schimek::mpi::execute_in_order([&]() {
      //    std::cout << "hashPEIndices within findDuplicates rank: " << env.rank()  << std::endl;
      //    std::cout << "elementsToMerge " << elementsToMerge << std::endl;
      //    for (size_t i = 0; i < hashPEIndices.size(); ++i)
      //    std::cout << i << " " << hashPEIndices[i] << std::endl;
      //    std::cout << "hashPEIndices within findDuplicates end" << std::endl;
      //    });


      tlx::multiway_merge(iteratorPairs.begin(), iteratorPairs.end(), mergedElements.begin(), elementsToMerge);

      //dss_schimek::mpi::execute_in_order([&]() {
      //    std::cout << "mergedElements rank: " << env.rank()  << std::endl;
      //    for (size_t i = 0; i < mergedElements.size(); ++i)
      //    std::cout << i << " " << mergedElements[i] << std::endl;
      //    std::cout << "mergedElements end" << std::endl;
      //    });


      std::vector<std::vector<size_t>> result_sets(intervalSizes.size());

      std::vector<size_t> counters(intervalSizes.size(), 0);

      for (auto elem : mergedElements)
        std::cout << elem << std::endl;

      HashPEIndex prevHashTriple = mergedElements.empty() ? HashPEIndex{0, 0} : mergedElements[0];
      bool duplicate = false;


      dss_schimek::mpi::execute_in_order([&]() {
          std::cout << "merge rank: " << env.rank() << std::endl;
      for(size_t i = 1; i < mergedElements.size(); ++i) {
        const HashPEIndex curHashTriple = mergedElements[i];
        std::cout << "prev: " << prevHashTriple << std::endl;
        std::cout << "cur: " << curHashTriple << std::endl;
        for (size_t j = 0; j < counters.size(); ++j) {
        std::cout << counters[j] << " ";
        std::cout << "resultset: " << j;
        for (auto elem : result_sets[j])
          std::cout << elem << " ";
        std::cout << std::endl;
        }
        std::cout << std::endl;
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
      });

      if (duplicate)
          result_sets[prevHashTriple.PEIndex].push_back(counters[prevHashTriple.PEIndex]++);

      std::vector<size_t> sendBuffer;
      sendBuffer.reserve(elementsToMerge);
      std::vector<size_t> sendCounts_;
        // write back to PEs
        for (size_t i = 0; i < result_sets.size(); ++i) {
          sendCounts_.push_back(result_sets[i].size());
          for (size_t j = 0; j < result_sets[i].size(); ++j) {
            sendBuffer.push_back(result_sets[i][j]);
          }
        }

      return dsss::mpi::AllToAllvSmall::alltoallv(sendBuffer.data(), sendCounts_);
    }
    
    void setUniqueElements(std::vector<bool>& duplicates, std::vector<size_t>& depth, const size_t curDepth, const std::vector<HashStringIndex>& originalMapping) {
      for (size_t i = 0; i < duplicates.size(); ++i)
        if (!duplicates[i]) {
          const size_t stringIndex = originalMapping[i].stringIndex;
          depth[stringIndex] = curDepth;
        }
    }

    void getIndicesOfDuplicates(std::vector<size_t>& duplicates, const std::vector<HashStringIndex>& originalMapping) {
      for (size_t i = 0; i < duplicates.size(); ++i) {
        const size_t stringIndex = originalMapping[duplicates[i]].stringIndex;
        duplicates[i] = stringIndex;
      }
    }
  };

  struct SendHashesAndPEIndexToFilter {
    using SendType = HashStringIndex;
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


  struct HashTripleComp {
    bool operator() (const HashTriple& lhs, const HashTriple& rhs) { 
      const bool firstCrit = lhs.hashValue < rhs.hashValue;
      const bool secondCrit = !firstCrit && lhs.stringIndex < rhs.stringIndex;
      return firstCrit || secondCrit;
    }
    bool operator() (const size_t& lhs, const HashTriple& rhs) { return lhs < rhs.hashValue; }
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

      const size_t bloomFilterSize = 100000;

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
          hash = ((hash << 5) + hash) + c; /* hash * 33 + c */
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
      std::vector<std::vector<size_t>> merge(const std::vector<HashTriple>& hashTriples) {
        dsss::mpi::environment env;
        std::vector<std::vector<size_t>> resultSets(env.size());
        if (hashTriples.empty())
          return resultSets;

        HashTriple prevHashTriple = hashTriples[0];
        bool duplicate = false;

        for(size_t i = 1; i < hashTriples.size(); ++i) {
          const HashTriple curHashTriple = hashTriples[i];
          if (prevHashTriple.hashValue == curHashTriple.hashValue) {
            resultSets[prevHashTriple.PEIndex].push_back(prevHashTriple.stringIndex);
            duplicate = true;
          } else if (duplicate) {
            resultSets[prevHashTriple.PEIndex].push_back(prevHashTriple.stringIndex);
            duplicate = false;
          }
          prevHashTriple = curHashTriple;
        }
        if (duplicate)
          resultSets[prevHashTriple.PEIndex].push_back(prevHashTriple.stringIndex);

        return resultSets;
      } 

      std::vector<size_t> getLocalDuplicateIndices(const std::vector<HashTriple>& hashTriples) {

        dsss::mpi::environment env;
        std::vector<std::vector<size_t>> globalDuplicates = merge(hashTriples);
        std::cout << "merge: rank: " << env.rank() << std::endl;
        
        std::vector<size_t> intervalSizes;
        size_t totalCounts = 0;
        for (const auto& vec : globalDuplicates) {
          intervalSizes.push_back(vec.size());
          totalCounts += vec.size();
        } 
        std::cout << "merge: rank: " << env.rank() << " totalCounts: " << totalCounts << std::endl;
        std::vector<size_t> sendBuffer;
        sendBuffer.reserve(totalCounts);

        for (size_t i = 0; i < globalDuplicates.size(); ++i) {
          for (size_t j = 0; j < globalDuplicates[i].size(); ++j) 
            sendBuffer.push_back(globalDuplicates[i][j]);
        }
        return dsss::mpi::AllToAllvSmall::alltoallv(sendBuffer.data(), intervalSizes);
      }

      std::vector<HashTriple> collectOnPE0(std::vector<HashTriple>& hashTriples) {
        std::vector<HashTriple> recvBuffer;
        int32_t localSize = hashTriples.size();
        std::vector<int32_t> intervalSizes = dsss::mpi::allgather(localSize);
        const int32_t totalSize = std::accumulate(intervalSizes.begin(), intervalSizes.end(), 0);
        std::vector<int32_t> sendDisplacements(intervalSizes.size(), 0);
        for (size_t i = 1; i < sendDisplacements.size(); ++i) {
          sendDisplacements[i] = sendDisplacements[i - 1] + intervalSizes[i - 1];
        }
        if (env.rank() == 0)
          recvBuffer.resize(totalSize);

        dss_schimek::mpi::execute_in_order([&]() {
            std::cout << "hashTriples rank: " << env.rank()  << std::endl;
            for (size_t i = 0; i < hashTriples.size(); ++i)
            std::cout << i << " " << hashTriples[i] << std::endl;
            std::cout << "hashTriples end" << std::endl;
            });


        dsss::mpi::data_type_mapper<HashTriple> dtm;
        MPI_Gatherv(
            hashTriples.data(),
            hashTriples.size(),
            dtm.get_mpi_type(),
            recvBuffer.data(),
            intervalSizes.data(),
            sendDisplacements.data(),
            dtm.get_mpi_type(),
            0,
            env.communicator());

        dss_schimek::mpi::execute_in_order([&]() {
            std::cout << "recvBuffer rank: " << env.rank()  << std::endl;
            for (size_t i = 0; i < recvBuffer.size(); ++i)
            std::cout << i << " " << recvBuffer[i] << std::endl;
            std::cout << "recvBuffer end" << std::endl;
            });


        return recvBuffer; 
      }

      std::vector<HashTriple> collectOnPE0(StringLcpPtr strptr,
          std::vector<size_t>& candidates,
          const size_t depth) {

        StringSet ss = strptr.active();
        std::vector<unsigned char> send_buffer;
        std::vector<unsigned char> receive_buffer_char;

        for (size_t j = 0; j < candidates.size(); ++j) {
          String str = ss[ss.begin() + candidates[j]];
          size_t string_length = ss.get_length(str) + 1; 
          std::copy_n(ss.get_chars(str, 0), string_length, std::back_inserter(send_buffer));
        }
        size_t numStrings = candidates.size();

        std::vector<size_t> recvCounts = dsss::mpi::allgather(numStrings);
        std::vector<size_t> stringIndices = dsss::mpi::allgatherv(candidates);
        std::vector<unsigned char> recvBuffer = dsss::mpi::allgatherv(send_buffer);
        size_t totalCount = std::accumulate(recvCounts.begin(), recvCounts.end(), 0);
        dss_schimek::StringLcpContainer<StringSet> container(std::move(recvBuffer));
        StringSet globalSet = container.make_string_set();

        dss_schimek::mpi::execute_in_order([&]() {
            std::cout << "print global Stringset rank: " << env.rank()  << std::endl;
            globalSet.print();
            std::cout << "end global Stringset " << std::endl;
            });

        std::vector<HashTriple> hashTriples;
        if (totalCount == 0) 
          return hashTriples;

        hashTriples.reserve(totalCount);
        size_t curOffset = 0; 
        auto begin = globalSet.begin();

        for (size_t curRank = 0; curRank < env.size(); ++curRank) {
         for (size_t i = 0; i < recvCounts[curRank]; ++i) {
          String curString = globalSet[begin + curOffset + i];
          const size_t curHash = hash(globalSet.get_chars(curString, 0), depth, bloomFilterSize);
          hashTriples.emplace_back(curHash, stringIndices[curOffset + i], curRank);
         } 
         curOffset += recvCounts[curRank];
        }

        dss_schimek::mpi::execute_in_order([&]() {
            std::cout << "print global hashes rank: " << env.rank()  << std::endl;
            for (const auto& hashTriple : hashTriples)
            std::cout << hashTriple << std::endl;
            std::cout << "end global hashes" << std::endl;
            });


        return hashTriples; 
      }

      std::vector<size_t> filter_simple(StringLcpPtr strptr, const size_t depth, 
          std::vector<size_t>& candidates, std::vector<size_t>& results) {

        std::vector<HashTriple> globalHashTriples = collectOnPE0(strptr, candidates, depth);
        std::vector<size_t> ownDuplicateIndices = getOwnDuplicates(globalHashTriples);

        for (size_t i = 0; i < candidates.size(); ++i) {
          results[candidates[i]] = depth;
        }
        return ownDuplicateIndices;
      } 

      std::vector<size_t> filter_naive(StringLcpPtr strptr, const size_t depth, const std::vector<size_t>& candidates, std::vector<size_t>& results) {
        dsss::mpi::environment env;
        std::cout << "filter_naive: rank: " << env.rank() << " depth: " << depth << std::endl;
        std::vector<HashTriple> hashTriples;
        const StringSet ss = strptr.active();
        const Iterator begin = ss.begin();
        for (size_t i = 0; i < candidates.size(); ++i) {
          String curString = ss[begin + candidates[i]];
          const size_t curHash = hash(ss.get_chars(curString, 0), depth, bloomFilterSize);
          hashTriples.emplace_back(curHash, candidates[i], env.rank());
        }

        std::sort(hashTriples.begin(), hashTriples.end());
        std::vector<HashTriple> recvTriples = collectOnPE0(hashTriples);
        std::sort(recvTriples.begin(), recvTriples.end());
        std::vector<size_t> localDuplicateIndices = getLocalDuplicateIndices(recvTriples); 

        dss_schimek::mpi::execute_in_order([&]() {
            std::cout << "localDuplicateIndices rank: " << env.rank()  << std::endl;
            for (size_t i = 0; i < localDuplicateIndices.size(); ++i)
            std::cout << i << " " << localDuplicateIndices[i] << std::endl;
            std::cout << "localDuplicateIndices end" << std::endl;
            });

        for (size_t i = 0; i < candidates.size(); ++i)
          results[candidates[i]] = depth;
        std::cout<< "set depth on rank: " << env.rank() << std::endl;

        return localDuplicateIndices;
      }

      std::vector<size_t> filter_new(StringLcpPtr strptr, const size_t depth, const std::vector<size_t>& candidates, std::vector<size_t>& results) {
        dsss::mpi::environment env;

        dss_schimek::mpi::execute_in_order([&]() {
            std::cout << "candidates rank: " << env.rank()  << std::endl;
            for (size_t i = 0; i < candidates.size(); ++i)
            std::cout << i << " " << candidates[i] << std::endl;
            std::cout << "candidates end" << std::endl;
            });

        std::vector<HashStringIndex> hashStringIndices;
        const StringSet ss = strptr.active();
        const Iterator begin = ss.begin();
        for (size_t i = 0; i < candidates.size(); ++i) {
          String curString = ss[begin + candidates[i]];

          const size_t curHash = hash(ss.get_chars(curString, 0), depth, bloomFilterSize);
          hashStringIndices.emplace_back(curHash, candidates[i]);
        }

        std::sort(hashStringIndices.begin(), hashStringIndices.end());
        std::cout<< "rank: " << env.rank() << " hashes are sorted " << std::endl;


        //auto [recvHashValues, intervalSizes] = SendPolicy<dsss::mpi::AllToAllvSmall>::sendToFilter(hashStringIndices, bloomFilterSize);
        auto p = SendPolicy<dsss::mpi::AllToAllvSmall>::sendToFilter(hashStringIndices, bloomFilterSize);
        std::vector<size_t> recvHashValues = p.first;
        std::vector<size_t> intervalSizes = p.second;

        //dss_schimek::mpi::execute_in_order([&]() {
        //    std::cout << "recvHashValues rank: " << env.rank()  << std::endl;
        //    for (size_t i = 0; i < recvHashValues.size(); ++i)
        //    std::cout << i << " " << recvHashValues[i] << std::endl;
        //    std::cout << "recvHashValues end" << std::endl;
        //    });


        std::vector<HashPEIndex> recvHashPEIndices = SendPolicy<dsss::mpi::AllToAllvSmall>::addPEIndex(recvHashValues, intervalSizes);

        dss_schimek::mpi::execute_in_order([&]() {
            std::cout << "recvHashPEIndices rank: " << env.rank()  << std::endl;
            for (size_t i = 0; i < recvHashPEIndices.size(); ++i)
            std::cout << i << " " << recvHashPEIndices[i] << std::endl;
            std::cout << "recvHashPEIndices end" << std::endl;
            });


    
        std::vector<size_t> indicesOfDuplicates = FindDuplicatesPolicy::findDuplicates(recvHashPEIndices, intervalSizes);
        dss_schimek::mpi::execute_in_order([&]() {
            std::cout << "indicesOfDuplicates1 rank: " << env.rank()  << std::endl;
            for (size_t i = 0; i < indicesOfDuplicates.size(); ++i)
            std::cout << i << " " << indicesOfDuplicates[i] << std::endl;
            std::cout << "indicesOfDuplicates1 end" << std::endl;
            });


        FindDuplicatesPolicy::getIndicesOfDuplicates(indicesOfDuplicates, hashStringIndices);

dss_schimek::mpi::execute_in_order([&]() {
            std::cout << "indicesOfDuplicates rank: " << env.rank()  << std::endl;
            for (size_t i = 0; i < indicesOfDuplicates.size(); ++i)
            std::cout << i << " " << indicesOfDuplicates[i] << std::endl;
            std::cout << "indicesOfDuplicates end" << std::endl;
            });


       
        // assume candidates are sorted
        for (size_t i = 0; i < candidates.size(); ++i) {
          results[candidates[i]] = depth;
        }

        return indicesOfDuplicates;

      }

      //std::vector<size_t> filter(StringLcpPtr strptr, const size_t depth, const std::vector<size_t>& candidates, std::vector<size_t>& results) {
      //  dsss::mpi::environment env;

      //  dss_schimek::mpi::execute_in_order([&]() {
      //      std::cout << "candidates rank: " << env.rank()  << std::endl;
      //      for (size_t i = 0; i < candidates.size(); ++i)
      //      std::cout << i << " " << candidates[i] << std::endl;
      //      std::cout << "candidates end" << std::endl;
      //      });

      //  std::vector<HashTriple> hashes;
      //  const StringSet ss = strptr.active();
      //  const Iterator begin = ss.begin();
      //  for (size_t i = 0; i < candidates.size(); ++i) {
      //    String curString = ss[begin + candidates[i]];

      //    const size_t curHash = hash(ss.get_chars(curString, 0), depth, bloomFilterSize);
      //    hashes.emplace_back(curHash, i, env.rank());
      //  }

      //  std::sort(hashes.begin(), hashes.end(), comp);
      //  std::cout<< "rank: " << env.rank() << " hashed are sorted " << std::endl;
      //  std::vector<size_t> intervalSizes= computeIntervals(hashes);

      //  dss_schimek::mpi::execute_in_order([&]() {
      //      std::cout << "filter boundaries for rank: " << env.rank() << std::endl;
      //      for (size_t i = 0; i < intervalSizes.size(); ++i)
      //      std::cout << i << " " << intervalSizes[i] << std::endl;
      //      });

      //  std::vector<HashTriple> recvHashTriples = AllToAllHashValuePolicy::allToAllHashTriples(hashes, intervalSizes);
      //  dss_schimek::mpi::execute_in_order([&]() {
      //      std::cout << "recv hashes for  rank: " << env.rank() << std::endl;
      //      for (size_t i = 0; i < recvHashTriples.size(); ++i) {
      //      std::string str = recvHashTriples[i];
      //      std::cout << i << " " << str << std::endl;
      //      }
      //      });


      //  std::sort(recvHashTriples.begin(), recvHashTriples.end(), comp);

      //  dss_schimek::mpi::execute_in_order([&]() {
      //      std::cout << "sorted hashes for  rank: " << env.rank() << std::endl;
      //      for (size_t i = 0; i < recvHashTriples.size(); ++i) {
      //      std::string str = recvHashTriples[i];
      //      std::cout << i << " " << str << std::endl;
      //      }
      //      });


      //  // find duplicates
      //  std::vector<std::vector<size_t>> duplicates(env.size(), std::vector<size_t>());
      //  HashTriple prevHashTriple = recvHashTriples.empty() ? HashTriple() : recvHashTriples[0];
      //  bool duplicate = false;

      //  for(size_t i = 1; i < recvHashTriples.size(); ++i) {
      //    const HashTriple curHashTriple = recvHashTriples[i];
      //    if (prevHashTriple.hashValue == curHashTriple.hashValue) {
      //      duplicates[prevHashTriple.PEIndex].push_back(prevHashTriple.stringIndex);
      //      duplicate = true;
      //    } else if (duplicate) {
      //      duplicates[prevHashTriple.PEIndex].push_back(prevHashTriple.stringIndex);
      //      duplicate = false;
      //    }
      //    prevHashTriple = curHashTriple;
      //  }
      //  std::cout << "rank: " << env.rank() << " rearrange indices" << std::endl;
      //  if (duplicate)
      //    duplicates[prevHashTriple.PEIndex].push_back(prevHashTriple.stringIndex);

      //  size_t size = 0;
      //  std::vector<size_t> sendCounts(env.size());
      //  for (size_t i = 0; i < duplicates.size(); ++i) {
      //    size += duplicates[i].size();
      //    sendCounts[i] = duplicates[i].size();
      //  }

      //  std::cout << "rank: " << env.rank() << " write back to sendBuffer " << std::endl;
      //  std::vector<size_t> sendBuffer;
      //  sendBuffer.reserve(size);
      //  // write back to PEs
      //  for (size_t i = 0; i < duplicates.size(); ++i) {
      //    for (size_t j = 0; j < duplicates[i].size(); ++j) {
      //      sendBuffer.push_back(duplicates[i][j]);
      //    }
      //  }

      //  std::vector<size_t> recvDuplicateIndices = dsss::mpi::AllToAllvSmall::alltoallv(sendBuffer.data(), sendCounts); // should I sort recvDuplicateIndices?

      //  dss_schimek::mpi::execute_in_order([&] () { 
      //      std::cout << "recv duplicates rank: " << env.rank() << std::endl;
      //      for (size_t i = 0; i < recvDuplicateIndices.size(); ++i) 
      //      std::cout << i << " " << recvDuplicateIndices[i] << std::endl;
      //      std::cout << "recv duplicates end " << std::endl;
      //      });

      //  for (size_t i = 0; i < candidates.size(); ++i) {
      //    results[candidates[i]] = depth;
      //  }
      //  for (size_t i = 0; i < recvDuplicateIndices.size(); ++i) {
      //    results[recvDuplicateIndices[i]] = 0;
      //  }

      //  return recvDuplicateIndices;
      //} 

      std::vector<size_t> filter(StringLcpPtr strptr) {
        return {};
      }
    };
}
