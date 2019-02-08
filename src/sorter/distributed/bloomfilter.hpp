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
    }

    static inline std::vector<SendType> sendToFilter(const std::vector<HashStringIndex>& hashes, size_t bloomfilterSize) {
      std::vector<size_t> sendValues = extractSendValues(hashes);
      std::vector<size_t> intervalSizes = computeIntervalSizes(sendValues, bloomfilterSize);
      return SendPolicy::allToAllv(sendValues, intervalSizes);
    }

    static inline HashPEIndex addPEIndex(const std::vector<size_t>& hashes, const std::vector<size_t>& intervalSizes) {
      std::vector<HashPEIndex> hashesPEIndex;
      hashesPEIndex.reserve(hashes.size());

      size_t curPE = 0;
      size_t curBoundary = intervalSizes[0];
      for (size_t i = 0; i < hashes.size(); ++i) {
        while (i == curBoundary) 
          curBoundary += intervalSizes[++curPE];
        hashesPEIndex.emplace_back(hashes[i], curPE);
      }  
    }
  };



  struct FindDuplicates {
    using DataType = HashPEIndex;
    
    static inline std::vector<bool> findDuplicates(std::vector<HashPEIndex>& hashPEIndices, std::vector<size_t>& intervalSizes) {
      using ConstIterator = std::vector<HashPEIndex>::const_iterator;
      using Iterator = std::vector<HashPEIndex>::iterator;
      using IteratorPair = std::pair<Iterator, Iterator>;
      std::vector<IteratorPair> iteratorPairs;

      size_t elementsToMerge = std::accumulate(intervalSizes.begin(), intervalSizes.end(), 0);
      std::vector<HashPEIndex> mergedElements(elementsToMerge);
      auto outputIt = std::back_inserter(mergedElements);
      Iterator it = hashPEIndices.begin(); 
      

      for (size_t i = 0; i < intervalSizes.size(); ++i) {
       iteratorPairs.emplace_back(it, it + intervalSizes[i]);
       it += intervalSizes[i];
      }

      tlx::multiway_merge(iteratorPairs.begin(), iteratorPairs.end(), mergedElements.begin(), elementsToMerge);

      std::vector<std::vector<bool>> result_sets;
      result_sets.reserve(intervalSizes.size());
      for (size_t i = 0; i < intervalSizes.size(); ++i) 
        result_sets.push_back(std::vector<bool>(intervalSizes[i], false));

      std::vector<size_t> counters(intervalSizes.size(), 0);

      for (auto elem : mergedElements)
        std::cout << elem << std::endl;

      HashPEIndex prevHashTriple = mergedElements.empty() ? HashPEIndex{0, 0} : mergedElements[0];
      bool duplicate = false;

      for(size_t i = 1; i < mergedElements.size(); ++i) {
        const HashPEIndex curHashTriple = mergedElements[i];
        if (prevHashTriple.hashValue == curHashTriple.hashValue) {
          result_sets[prevHashTriple.PEIndex][counters[prevHashTriple.PEIndex]++] = true;
          duplicate = true;
        } else if (duplicate) {
          result_sets[prevHashTriple.PEIndex][counters[prevHashTriple.PEIndex]++] = true;
          duplicate = false;
        }
        prevHashTriple = curHashTriple;
      }

      if (duplicate)
        result_sets[prevHashTriple.PEIndex][counters[prevHashTriple.PEIndex]++] = true;

      std::vector<bool> sendBuffer;
      sendBuffer.reserve(elementsToMerge);
        // write back to PEs
        for (size_t i = 0; i < result_sets.size(); ++i) {
          for (size_t j = 0; j < result_sets[i].size(); ++j) {
            sendBuffer.push_back(result_sets[i][j]);
          }
        }
      dsss::mpi::AllToAllvSmall::alltoallv(sendBuffer.data(), intervalSizes);
    }
    
    void setUniqueElements(std::vector<bool>& duplicates, std::vector<size_t>& depth, const size_t curDepth, const std::vector<HashStringIndex>& originalMapping) {
      for (size_t i = 0; i < duplicates.size(); ++i)
        if (!duplicates[i]) {
          const size_t stringIndex = originalMapping[i].stringIndex;
          depth[stringIndex] = curDepth;
        }
    }

    std::vector<size_t> getIndicesOfDuplicates(const std::vector<bool>& duplicates, const std::vector<HashStringIndex>& originalMapping) {
      std::vector<size_t> duplicateIndices;
      for (size_t i = 0; i < duplicates.size(); ++i)
        if (!duplicates[i]) {
          const size_t stringIndex = originalMapping[i].stringIndex;
          duplicateIndices.push_back(stringIndex);
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
      }
  };


  template <typename StringSet, typename AllToAllHashValuePolicy>
    class BloomFilter : private AllToAllHashValuePolicy {
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

      std::vector<size_t> filter(StringLcpPtr strptr, size_t depth, const std::vector<size_t>& candidates, std::vector<size_t>& results) {
        dsss::mpi::environment env;

        dss_schimek::mpi::execute_in_order([&]() {
            std::cout << "candidates rank: " << env.rank()  << std::endl;
            for (size_t i = 0; i < candidates.size(); ++i)
            std::cout << i << " " << candidates[i] << std::endl;
            std::cout << "candidates end" << std::endl;
            });

        std::vector<HashTriple> hashes;
        const StringSet ss = strptr.active();
        const Iterator begin = ss.begin();
        for (size_t i = 0; i < candidates.size(); ++i) {
          String curString = ss[begin + candidates[i]];

          const size_t curHash = hash(ss.get_chars(curString, 0), depth, bloomFilterSize);
          hashes.emplace_back(curHash, i, env.rank());
        }

        std::sort(hashes.begin(), hashes.end(), comp);
        std::cout<< "rank: " << env.rank() << " hashed are sorted " << std::endl;
        std::vector<size_t> intervalSizes= computeIntervals(hashes);

        dss_schimek::mpi::execute_in_order([&]() {
            std::cout << "filter boundaries for rank: " << env.rank() << std::endl;
            for (size_t i = 0; i < intervalSizes.size(); ++i)
            std::cout << i << " " << intervalSizes[i] << std::endl;
            });

        std::vector<HashTriple> recvHashTriples = AllToAllHashValuePolicy::allToAllHashTriples(hashes, intervalSizes);
        dss_schimek::mpi::execute_in_order([&]() {
            std::cout << "recv hashes for  rank: " << env.rank() << std::endl;
            for (size_t i = 0; i < recvHashTriples.size(); ++i) {
            std::string str = recvHashTriples[i];
            std::cout << i << " " << str << std::endl;
            }
            });


        std::sort(recvHashTriples.begin(), recvHashTriples.end(), comp);

        dss_schimek::mpi::execute_in_order([&]() {
            std::cout << "sorted hashes for  rank: " << env.rank() << std::endl;
            for (size_t i = 0; i < recvHashTriples.size(); ++i) {
            std::string str = recvHashTriples[i];
            std::cout << i << " " << str << std::endl;
            }
            });


        // find duplicates
        std::vector<std::vector<size_t>> duplicates(env.size(), std::vector<size_t>());
        HashTriple prevHashTriple = recvHashTriples.empty() ? HashTriple() : recvHashTriples[0];
        bool duplicate = false;

        for(size_t i = 1; i < recvHashTriples.size(); ++i) {
          const HashTriple curHashTriple = recvHashTriples[i];
          if (prevHashTriple.hashValue == curHashTriple.hashValue) {
            duplicates[prevHashTriple.PEIndex].push_back(prevHashTriple.stringIndex);
            duplicate = true;
          } else if (duplicate) {
            duplicates[prevHashTriple.PEIndex].push_back(prevHashTriple.stringIndex);
            duplicate = false;
          }
          prevHashTriple = curHashTriple;
        }
        std::cout << "rank: " << env.rank() << " rearrange indices" << std::endl;
        if (duplicate)
          duplicates[prevHashTriple.PEIndex].push_back(prevHashTriple.stringIndex);

        size_t size = 0;
        std::vector<size_t> sendCounts(env.size());
        for (size_t i = 0; i < duplicates.size(); ++i) {
          size += duplicates[i].size();
          sendCounts[i] = duplicates[i].size();
        }

        std::cout << "rank: " << env.rank() << " write back to sendBuffer " << std::endl;
        std::vector<size_t> sendBuffer;
        sendBuffer.reserve(size);
        // write back to PEs
        for (size_t i = 0; i < duplicates.size(); ++i) {
          for (size_t j = 0; j < duplicates[i].size(); ++j) {
            sendBuffer.push_back(duplicates[i][j]);
          }
        }

        std::vector<size_t> recvDuplicateIndices = dsss::mpi::AllToAllvSmall::alltoallv(sendBuffer.data(), sendCounts); // should I sort recvDuplicateIndices?

        dss_schimek::mpi::execute_in_order([&] () { 
            std::cout << "recv duplicates rank: " << env.rank() << std::endl;
            for (size_t i = 0; i < recvDuplicateIndices.size(); ++i) 
            std::cout << i << " " << recvDuplicateIndices[i] << std::endl;
            std::cout << "recv duplicates end " << std::endl;
            });

        for (size_t i = 0; i < candidates.size(); ++i) {
          results[candidates[i]] = depth;
        }
        for (size_t i = 0; i < recvDuplicateIndices.size(); ++i) {
          results[recvDuplicateIndices[i]] = 0;
        }

        return recvDuplicateIndices;
      } 

      std::vector<size_t> filter(StringLcpPtr strptr) {
        return {};
      }
    };
}