#include <vector>
#include <algorithm>
#include <functional>

#include <tlx/siphash.hpp>
#include <tlx/sort/strings/radix_sort.hpp>

#include "strings/stringcontainer.hpp"
#include "strings/stringset.hpp"
#include "sorter/distributed/bloomfilter.hpp"
#include "util/random_string_generator.hpp"

namespace dss_schimek {
  namespace tests {
    template <typename HashFunctor, typename StringSet>
      class BloomfilterTester {
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
        HashFunctor hasher;
        size_t bloomFilterSize = 1000;

        HashesEOSHashes generateHashes(StringSet ss, const std::vector<size_t>& candidates, const size_t depth) {
          HashesEOSHashes hashContainer;
          std::vector<size_t>& hashes = hashContainer.hashes;
          std::vector<size_t>& eosHashes = hashContainer.eosHashes;

          auto begin = ss.begin();
          for (const size_t curCandidate : candidates) {
            String curString = ss[begin + curCandidate];
            const size_t length = ss.get_length(curString);
            if (depth > length) {
              const size_t hash = hasher(ss.get_chars(curString, 0), length, bloomFilterSize);
              eosHashes.push_back(hash);
            } else {
              const size_t hash = hasher(ss.get_chars(curString, 0), depth, bloomFilterSize);
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

        DuplicateUniqueHashes filter(StringSet ss, const std::vector<size_t>& candidates, const size_t depth) {
          HashesEOSHashes hashContainer = generateHashes(ss, candidates, depth);  
          std::cout << "hashes" << std::endl;
          for (size_t i = 0; i < hashContainer.hashes.size(); ++i)
            std::cout << i << " " << hashContainer.hashes[i] << std::endl;
          std::cout << "EOShashes" << std::endl;
          for (size_t i = 0; i < hashContainer.eosHashes.size(); ++i)
            std::cout << i << " " << hashContainer.eosHashes[i] << std::endl;

          std::vector<size_t> globalHashes = dsss::mpi::allgatherv(hashContainer.hashes);
          DuplicateUniqueHashes result = getDuplicateUniqueHashes(globalHashes);
          std::copy_n(hashContainer.eosHashes.begin(), hashContainer.eosHashes.size(), std::back_inserter(result.uniques));
          // -> can contain duplicates if eosHash and globalHash has same hash value but this is no problem
          std::cout << "dups: " << std::endl;
          for (size_t i = 0; i < result.duplicates.size(); ++i)
            std::cout << i << " " << result.duplicates[i] << std::endl;
          std::cout << "uniques:" << std::endl;
          for (size_t i = 0; i < result.uniques.size(); ++i)
            std::cout << i << " " << result.uniques[i] << std::endl;
        return result;
        }
      };
  }
}

std::vector<size_t> getDifference(std::vector<size_t>& superSet, std::vector<size_t>& subSet) {
  std::vector<size_t> differenceSet;
  differenceSet.reserve(superSet.size());
  std::sort(superSet.begin(), superSet.end());
  std::sort(subSet.begin(), subSet.end());
  std::set_difference(superSet.begin(), superSet.end(), 
                      subSet.begin(), subSet.end(),
                      std::back_inserter(differenceSet));
}

int main() {
  dsss::mpi::environment env;
  using namespace dss_schimek;
  using namespace dss_schimek::tests;

  using StringSet = UCharLengthStringSet;
  using Hasher = std::function<size_t(const unsigned char*, size_t, size_t)>;
  using Generator = DNRatioGenerator<StringSet>;

  using StringLcpPtr = typename tlx::sort_strings_detail::StringLcpPtr<StringSet, size_t>;

  using GolombPolicy = AllToAllHashesNaive;
  
  const size_t size = 100;
  const size_t stringLength = 10;
  const double dToNRatio = 0.0;
  const size_t startDepth = 2;

  BloomfilterTester<Hasher, StringSet> bftest;
  BloomFilter<StringSet, FindDuplicates, SendOnlyHashesToFilter<GolombPolicy>> bloomFilter;
  Generator container(size, stringLength, dToNRatio);
  StringLcpPtr localStringPtr = container.make_string_lcp_ptr();
  StringSet ss = container.make_string_set();
  Hasher hasher = [](const unsigned char* str, size_t maxDepth, size_t m) {
    //return tlx::siphash(chars, depth) % bloomFilterSize;
    size_t hash = 5381;
    size_t c = 0, i = 0;

    while ((c = *str++) && i < maxDepth) {
      hash = ((hash << 5) + hash) + c * 33; /* hash * 33 + c */
      ++i;
    }
    return hash % m;
  };
  ss.print();
  bftest.hasher = hasher;
  bftest.bloomFilterSize = bloomFilter.bloomFilterSize;

  std::vector<size_t> candidates_init;
  candidates_init.reserve(size);
  std::generate_n(std::back_inserter(candidates_init), size, [size]() {static size_t i = 0; return i++;});


  std::vector<size_t> results(size, 0);
  std::vector<size_t> candidates = bloomFilter.filter(localStringPtr, startDepth, results);

  bftest.filter(ss, candidates_init, 2);

  env.finalize();
  return 0;
}
