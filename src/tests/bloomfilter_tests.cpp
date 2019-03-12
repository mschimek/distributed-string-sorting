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
    template <typename HashPolicy, typename StringSet>
      class BloomfilterTester : HashPolicy{
        using String = typename StringSet::String;


        struct HashesEOSHashes {
          std::vector<size_t> hashes;
          std::vector<size_t> eosHashes; // End Of String hashes, i.e. hashes from strings that
          // shorter than current depth. These strings are hashed
          // til their end
        };
        struct DuplicateUniqueHashes {
          std::vector<size_t> duplicates;
          std::vector<size_t> uniques;
        };
        public:
        BloomfilterTester(size_t bloomFilterSize) : bloomFilterSize(bloomFilterSize) {}

        DuplicateUniqueHashes filter(StringSet ss, const std::vector<size_t>& candidates, const size_t depth) {
          HashesEOSHashes hashContainer = generateHashes(ss, candidates, depth);  
          std::vector<size_t> globalHashes = dsss::mpi::allgatherv(hashContainer.hashes);
          DuplicateUniqueHashes result = getDuplicateUniqueHashes(globalHashes);
          std::copy_n(hashContainer.eosHashes.begin(), hashContainer.eosHashes.size(), std::back_inserter(result.uniques));
          // -> can contain duplicates if eosHash and globalHash has same hash value but this is no problem
          
          return result;
        }

        private:
        size_t bloomFilterSize;

        HashesEOSHashes generateHashes(StringSet ss, const std::vector<size_t>& candidates, const size_t depth) {
          HashesEOSHashes hashContainer;
          std::vector<size_t>& hashes = hashContainer.hashes;
          std::vector<size_t>& eosHashes = hashContainer.eosHashes;

          auto begin = ss.begin();
          for (const size_t curCandidate : candidates) {
            String curString = ss[begin + curCandidate];
            const size_t length = ss.get_length(curString);
            if (depth > length) {
              const size_t hash = HashPolicy::hash(ss.get_chars(curString, 0), length, bloomFilterSize);
              eosHashes.push_back(hash);
            } else {
              const size_t hash = HashPolicy::hash(ss.get_chars(curString, 0), depth, bloomFilterSize);
              hashes.push_back(hash);
            }
          }
          return hashContainer;
        }

        DuplicateUniqueHashes getDuplicateUniqueHashes(std::vector<size_t>& hashes) {
          DuplicateUniqueHashes hashContainer;
          std::vector<size_t>& duplicates = hashContainer.duplicates;
          std::vector<size_t>& uniques = hashContainer.uniques;
          duplicates.reserve(hashes.size());
          uniques.reserve(hashes.size());

          // special cases
          if (hashes.empty())
            return hashContainer;
          if (hashes.size() == 1) {
            uniques.push_back(hashes.front());
            return hashContainer;
          }

          std::sort(hashes.begin(), hashes.end());

          for (size_t i = 0; i + 1 < hashes.size();) {
            const size_t curHash = hashes[i];
            size_t indexFirstDiffHash = i + 1;
            while (indexFirstDiffHash < hashes.size() && curHash == hashes[indexFirstDiffHash])
              ++indexFirstDiffHash; // advance til hashes are different or end of vector

            if (indexFirstDiffHash - i == 1)
              uniques.push_back(curHash);
            else
              duplicates.push_back(curHash);

            i = indexFirstDiffHash;
          }
          if (hashes.back() != hashes[hashes.size() - 2])
            uniques.push_back(hashes.back());

          return hashContainer;
        } 
      };

    std::vector<size_t> getDifference(std::vector<size_t> minuendSet, std::vector<size_t> subtrahendSet) {
      std::vector<size_t> setDifference;
      setDifference.reserve(minuendSet.size());

      std::sort(minuendSet.begin(), minuendSet.end());
      std::sort(subtrahendSet.begin(), subtrahendSet.end());

      std::set_difference(minuendSet.begin(), minuendSet.end(), 
          subtrahendSet.begin(), subtrahendSet.end(),
          std::back_inserter(setDifference));
      return setDifference;
    }

    struct CharacterBasesHash {
      static inline size_t hash(const unsigned char* str, size_t maxDepth, size_t m) {
        size_t hash = 5381;
        size_t c = 0, i = 0;

        while ((c = *str++) && i < maxDepth) {
          hash = ((hash << 5) + hash) + c * 33; /* hash * 33 + c */
          ++i;
        }
        return hash % m;
      }
    };

    template<typename StringSet, typename Hasher>
      std::vector<size_t> getHashValues(StringSet ss, const std::vector<size_t>& indicesToHash, const size_t depth, const size_t bloomFilterSize) {
        std::vector<size_t> hashes;
        hashes.reserve(indicesToHash.size());
        auto begin = ss.begin();
        for (const size_t index : indicesToHash) {
          const auto str = ss[begin + index];
          const auto chars = ss.get_chars(str, 0);
          const size_t stringLength = ss.get_length(str);
          const size_t actualDepth = stringLength < depth ? stringLength : depth;
          const size_t hash = Hasher::hash(chars, actualDepth, bloomFilterSize);
          hashes.push_back(hash);
        }
        return hashes;
      }

    void testResult(std::vector<size_t> globalDuplicateHashes,
        std::vector<size_t> globalUniqueHashes,
        const std::vector<size_t>& hashesDuplicates, 
        const std::vector<size_t>& hashesUniqueElems) {

      std::sort(globalUniqueHashes.begin(), globalUniqueHashes.end());
      const bool foundAllUniqueHashes = std::all_of(hashesUniqueElems.begin(), hashesUniqueElems.end(), [&] (const size_t hash) {
          return std::binary_search(globalUniqueHashes.begin(), globalUniqueHashes.end(), hash);
          });

      // verify that all duplicate hashes found by the bloomfilter are also found by the tester
      std::sort(globalDuplicateHashes.begin(), globalDuplicateHashes.end());
      const bool foundAllDuplicateHashes = std::all_of(hashesDuplicates.begin(), hashesDuplicates.end(), [&] (const size_t hash) {
          return std::binary_search(globalDuplicateHashes.begin(), globalDuplicateHashes.end(), hash);
          });

      tlx_die_unless(foundAllUniqueHashes);
      tlx_die_unless(foundAllDuplicateHashes);
    }

    template<typename GolombPolicy>
      void bloomfilter_test(const size_t size, const size_t stringLength, const double dToNRatio) {
        using namespace dss_schimek;
        using namespace dss_schimek::tests;

        using StringSet = UCharLengthStringSet;
        using Generator = DNRatioGenerator<StringSet>;
        using StringLcpPtr = typename tlx::sort_strings_detail::StringLcpPtr<StringSet, size_t>;

        const size_t startDepth = 2;

        Generator container(size, stringLength, dToNRatio);
        StringLcpPtr localStringPtr = container.make_string_lcp_ptr();
        StringSet ss = container.make_string_set();
        const size_t actualSize = ss.size();

        BloomFilter<StringSet, FindDuplicates, SendOnlyHashesToFilter<GolombPolicy>> bloomFilter;
        BloomfilterTester<CharacterBasesHash, StringSet> bftest(bloomFilter.bloomFilterSize);

        std::vector<size_t> results(actualSize, 0);
        std::vector<size_t> duplicates(actualSize);
        std::iota(duplicates.begin(), duplicates.end(), 0);
        duplicates.reserve(actualSize); // all indices are duplicates at the beginning
        //std::generate_n(std::back_inserter(duplicates), actualSize, []() {static size_t i = 0; return i++;});

        for (size_t depth = startDepth; depth < std::numeric_limits<size_t>::max(); depth *= 2) {
          const auto candidates = duplicates;

          // run Bloomfilter and tests
          duplicates = bloomFilter.filter(localStringPtr, depth, candidates, results);
          auto duplicatesUniques = bftest.filter(ss, candidates, depth);

          std::vector<size_t>& globalDuplicateHashes = duplicatesUniques.duplicates;
          std::vector<size_t>& globalUniqueHashes = duplicatesUniques.uniques;


          std::vector<size_t> hashesDuplicates = getHashValues<StringSet, CharacterBasesHash>(ss, duplicates, depth, bloomFilter.bloomFilterSize);
          std::vector<size_t> indicesUniqueElems = getDifference(candidates, duplicates);
          std::vector<size_t> hashesUniqueElems = getHashValues<StringSet, CharacterBasesHash>(ss, indicesUniqueElems, depth, bloomFilter.bloomFilterSize);

          // basic tests
          tlx_die_unless(depth != startDepth || candidates.size() == actualSize); // depth == startDepth -> all indices are candidates
          tlx_die_unless(indicesUniqueElems.size() + duplicates.size() == candidates.size());

          tlx_die_unless(indicesUniqueElems.size() == 0 || hashesUniqueElems.size() > 0); // ther are some uniques -> there are some unique Hashes
          tlx_die_unless(indicesUniqueElems.size() <= hashesUniqueElems.size());

          tlx_die_unless(duplicates.size() == 0 || hashesDuplicates.size() > 0); // ther are some dup -> there are some unique dup
          tlx_die_unless(duplicates.size() <= hashesDuplicates.size());

          // complete test: verify that all unique hashes found by the bloomfilter are also found by the tester
          testResult(globalDuplicateHashes, globalUniqueHashes, hashesDuplicates, hashesUniqueElems);

          bool noMoreCandidates = duplicates.empty();
          const bool allEmpty = dsss::mpi::allreduce_and(noMoreCandidates);
          if (allEmpty)
            break;
        }
      }
  }
}

int main() {
  using namespace dss_schimek;
  using namespace tests;
  using measurement::MeasuringTool;
  MeasuringTool& measuringTool = MeasuringTool::measuringTool();
  measuringTool.disable();

  dsss::mpi::environment env;

  const size_t size = 10000;
  const size_t stringLength = 500;
  std::vector<double> dToNRatios{0.0, 0.2, 0.4, 0.6, 0.8, 1.0};

  for (const double dToNRatio : dToNRatios) {
    std::cout << "test with dToNRatio: " << dToNRatio << std::endl;

    std::cout << "test with " << AllToAllHashesNaive::getName() << std::endl; 
    bloomfilter_test<AllToAllHashesNaive>(size, stringLength, dToNRatio);
    std::cout << "test with " << AllToAllHashesGolomb::getName() << std::endl; 
    bloomfilter_test<AllToAllHashesGolomb>(size, stringLength, dToNRatio);
    std::cout << "test with " << AllToAllHashValuesPipeline::getName() << std::endl; 
    bloomfilter_test<AllToAllHashValuesPipeline>(size, stringLength, dToNRatio);

  }
  std::cout << "tests completed successfully!" << std::endl;
  env.finalize();
  return 0;
}
