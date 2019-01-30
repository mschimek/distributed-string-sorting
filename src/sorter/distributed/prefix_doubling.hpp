#include <algorithm>
#include <type_traits>
#include <random>

#include "strings/stringptr.hpp"
#include "strings/stringset.hpp"
#include "strings/stringtools.hpp"

#include "sorter/local/strings/insertion_sort_unified.hpp"
#include "sorter/local/strings/multikey_quicksort_unified.hpp"
#include "sorter/local/strings/radix_sort_unified.hpp"

#include "mpi/alltoall.hpp"
#include "mpi/allgather.hpp"
#include "mpi/synchron.hpp"
#include "mpi/is_sorted.hpp"
#include "mpi/byte_encoder.hpp"

#include "merge/stringtools.hpp"
#include "merge/bingmann-lcp_losertree.hpp"

#include "util/timer.hpp"

namespace dss_schimek {

  static constexpr bool debug = false;

  template <typename StringSet>
    class  SampleSplittersNumStringsPolicy
    {
      public:
        static std::string getName() {
          return "NumStrings";
        }
      protected:
        std::vector<typename StringSet::Char> sample_splitters(const StringSet& ss,
            dsss::mpi::environment env = dsss::mpi::environment()) {

          using Char = typename StringSet::Char;
          using String = typename StringSet::String;

          const size_t local_num_strings = ss.size();
          const size_t nr_splitters = std::min<size_t>(env.size() - 1, local_num_strings);
          const size_t splitter_dist = local_num_strings / (nr_splitters + 1);
          std::vector<Char> raw_splitters;

          for (size_t i = 1; i <= nr_splitters; ++i) {
            const String splitter = ss[ss.begin() + i * splitter_dist];
            std::copy_n(ss.get_chars(splitter, 0), ss.get_length(splitter) + 1,
                std::back_inserter(raw_splitters));
          }
          return raw_splitters; 
        }
    };

  template <typename StringSet>
    class  SampleSplittersNumCharsPolicy
    {
      public: 
        static std::string getName() {
          return "NumChars";
        }
      protected:
        std::vector<typename StringSet::Char> sample_splitters(const StringSet& ss,
            dsss::mpi::environment env = dsss::mpi::environment()) {

          using Char = typename StringSet::Char;
          using String = typename StringSet::String;

          const size_t num_chars = std::accumulate(ss.begin(), ss.end(), 0, 
              [&ss](const size_t& sum, const String& str) {
              return sum + ss.get_length(str);
              });

          const size_t local_num_strings = ss.size();
          const size_t nr_splitters = std::min<size_t>(env.size() - 1, local_num_strings);
          const size_t splitter_dist = num_chars / (nr_splitters + 1);
          std::vector<Char> raw_splitters;

          size_t string_index = 0;
          for (size_t i = 1; i <= nr_splitters; ++i) {
            size_t num_chars_seen = 0;
            while (num_chars_seen < splitter_dist) {
              num_chars_seen += ss.get_length(ss[ss.begin() + string_index]);
              ++string_index;
            }

            const String splitter = ss[ss.begin() + string_index - 1];
            std::copy_n(ss.get_chars(splitter, 0), ss.get_length(splitter) + 1,
                std::back_inserter(raw_splitters));
          }
          return raw_splitters; 
        }
    };


  template <typename StringSet>
    StringLcpContainer<StringSet> choose_splitters(
        const StringSet& ss, 
        std::vector<typename StringSet::Char>& all_splitters,
        dsss::mpi::environment env = dsss::mpi::environment())
    {
      using Char = typename StringSet::Char;
      using String = typename StringSet::String;

      StringLcpContainer<StringSet> all_splitters_cont(std::move(all_splitters));
      dss_schimek::StringLcpPtr all_splitters_strptr = all_splitters_cont.make_string_lcp_ptr();
      const StringSet& all_splitters_set = all_splitters_strptr.active();

      insertion_sort(all_splitters_strptr, 0, 0);

      const size_t nr_splitters = std::min<std::size_t>(env.size() - 1, all_splitters_set.size());
      const size_t splitter_dist = all_splitters_set.size() / (nr_splitters + 1);

      std::vector<Char> raw_chosen_splitters;
      for (std::size_t i = 1; i <= nr_splitters; ++i) {
        const auto begin = all_splitters_set.begin();
        const String splitter = all_splitters_set[begin + i * splitter_dist];
        std::copy_n(ss.get_chars(splitter, 0), ss.get_length(splitter) + 1,
            std::back_inserter(raw_chosen_splitters));
      }
      return StringLcpContainer<StringSet>(std::move(raw_chosen_splitters));
    }

  template <typename StringSet>
    inline std::vector<size_t> compute_interval_sizes(const StringSet& ss, 
        const StringSet& splitters,
        dsss::mpi::environment env = dsss::mpi::environment())
    {
      std::vector<size_t> interval_sizes;
      interval_sizes.reserve(splitters.size());

      size_t nr_splitters = std::min<size_t>(env.size() - 1, ss.size());
      size_t splitter_dist = ss.size() / (nr_splitters + 1);
      size_t element_pos = 0;

      for (std::size_t i = 0; i < splitters.size(); ++i) {
        element_pos = (i + 1) * splitter_dist;

        while(element_pos > 0 && !dss_schimek::leq(
              ss.get_chars(ss[ss.begin() + element_pos], 0), 
              splitters.get_chars(splitters[splitters.begin() + i], 0))) 
        { --element_pos; }

        while (element_pos < ss.size() && dss_schimek::leq(
              ss.get_chars(ss[ss.begin() + element_pos], 0),
              splitters.get_chars(splitters[splitters.begin() + i], 0))) 
        { ++element_pos; }

        interval_sizes.emplace_back(element_pos);
      }
      interval_sizes.emplace_back(ss.size());
      for (std::size_t i = interval_sizes.size() - 1; i > 0; --i) {
        interval_sizes[i] -= interval_sizes[i - 1];
      }
      return interval_sizes;
    } 

  template <typename StringSet>
    inline static int binarySearch(const StringSet& ss, typename StringSet::CharIterator elem) {
      using String = typename StringSet::String;
      using CharIt = typename StringSet::CharIterator;

      auto left = ss.begin();
      auto right = ss.end();


      while(left != right) {
        size_t dist = (right - left) / 2;
        String curStr = ss[left + dist];
        size_t curLcp = 0;
        int res = dss_schimek::scmp(ss.get_chars(curStr, 0), elem);
        if (res < 0) {
          left = left + dist + 1;
        } else if (res == 0) {
          return left + dist - ss.begin();
        } else {
          right = left + dist;
        }
      }
      return left -ss.begin();
    }

  template <typename StringSet>
    inline std::vector<size_t> compute_interval_binary(const StringSet& ss, 
        const StringSet& splitters,
        dsss::mpi::environment env = dsss::mpi::environment())
    {
      using String = typename StringSet::String;
      using CharIt = typename StringSet::CharIterator;
      std::vector<size_t> interval_sizes;
      interval_sizes.reserve(splitters.size());

      size_t nr_splitters = std::min<size_t>(env.size() - 1, ss.size());
      size_t splitter_dist = ss.size() / (nr_splitters + 1);
      size_t element_pos = 0;

      for (std::size_t i = 0; i < splitters.size(); ++i) {
        //element_pos = (i + 1) * splitter_dist;

        //while(element_pos > 0 && !dss_schimek::leq(
        //      ss.get_chars(ss[ss.begin() + element_pos], 0), 
        //      splitters.get_chars(splitters[splitters.begin() + i], 0))) 
        //{ --element_pos; }

        //while (element_pos < ss.size() && dss_schimek::leq(
        //      ss.get_chars(ss[ss.begin() + element_pos], 0),
        //      splitters.get_chars(splitters[splitters.begin() + i], 0))) 
        //{ ++element_pos; }

        CharIt splitter = splitters.get_chars(splitters[splitters.begin() + i], 0);
        size_t pos = binarySearch(ss, splitter);
        interval_sizes.emplace_back(pos);
      }
      interval_sizes.emplace_back(ss.size());
      for (std::size_t i = interval_sizes.size() - 1; i > 0; --i) {
        interval_sizes[i] -= interval_sizes[i - 1];
      }
      return interval_sizes;
    }

  static inline void print_interval_sizes(const std::vector<size_t>& sent_interval_sizes,
      const std::vector<size_t>& recv_interval_sizes,
      dsss::mpi::environment env = dsss::mpi::environment())
  {
    constexpr bool print_interval_details = true;
    if constexpr (print_interval_details) {
      for (std::int32_t rank = 0; rank < env.size(); ++rank) {
        if (env.rank() == rank) {
          std::size_t total_size = 0;
          std::cout << "### Sending interval sizes on PE " << rank << std::endl;
          for (const auto is : sent_interval_sizes) {
            total_size += is;
            std::cout << is << ", ";
          }
          std::cout << "Total size: " << total_size << std::endl;
        }
        env.barrier();
      }
      for (std::int32_t rank = 0; rank < env.size(); ++rank) {
        if (env.rank() == rank) {
          std::size_t total_size = 0;
          std::cout << "### Receiving interval sizes on PE " << rank << std::endl;
          for (const auto is : recv_interval_sizes) {
            total_size += is;
            std::cout << is << ", ";
          }
          std::cout << "Total size: " << total_size << std::endl;
        }
        env.barrier();
      }
      if (env.rank() == 0) { std::cout << std::endl; }
    }
  }

  template <typename StringLcpContainer>
    static inline std::vector<std::pair<size_t, size_t>> compute_ranges_and_set_lcp_at_start_of_range(
        StringLcpContainer& recv_string_cont,
        std::vector<size_t>& recv_interval_sizes,
        dsss::mpi::environment env = dsss::mpi::environment())
    {
      std::vector<std::pair<size_t, size_t>> ranges;
      for(size_t i = 0, offset = 0; i < env.size(); ++i)
      {
        if (recv_interval_sizes[i] == 0)
        {
          ranges.emplace_back(0, 0);
          continue;
        }
        *(recv_string_cont.lcp_array() + offset) = 0;
        ranges.emplace_back(offset, recv_interval_sizes[i]);
        offset += recv_interval_sizes[i];
      }

      if constexpr (debug)
        dss_schimek::mpi::execute_in_order([&] () {
            std::cout << "rank: " << env.rank() << " pairs:" << std::endl;
            for(size_t i = 0; i < env.size(); ++i)
            std::cout << i << " " << ranges[i].first << " " << ranges[i].second << std::endl;
            });

      return ranges;
    }
  template<typename AllToAllStringPolicy, size_t K, typename StringSet>
    static inline StringLcpContainer<StringSet> merge(
        dss_schimek::StringLcpContainer<StringSet>&& recv_string_cont,
        const std::vector<std::pair<size_t, size_t>>& ranges,
        const size_t num_recv_elems) {



      std::vector<typename StringSet::String> sorted_string(recv_string_cont.size());
      std::vector<size_t> sorted_lcp(recv_string_cont.size());
      StringSet ss = recv_string_cont.make_string_set();
      //dss_schimek::mpi::execute_in_order([&] () {
      //    dsss::mpi::environment env;
      //    env.barrier();
      //    std::cout << "before merge: \n rank: " << env.rank() << std::endl;
      //    ss.print();
      //    for (size_t i = 0; i < ss.size(); ++i)
      //    std::cout << i << " " << recv_string_cont.lcps()[i] << std::endl;
      //    });
      dss_schimek::StringLcpPtrMergeAdapter<StringSet> mergeAdapter(ss, recv_string_cont.lcp_array());
      dss_schimek::LcpStringLoserTree_<K, StringSet> loser_tree(mergeAdapter, ranges.data());
      StringSet sortedSet(sorted_string.data(), sorted_string.data() + sorted_string.size());
      dss_schimek::StringLcpPtrMergeAdapter out_(sortedSet, sorted_lcp.data());
      //dss_schimek::mpi::execute_in_order([&] () {
      //    dsss::mpi::environment env;
      //    env.barrier();

      //    std::cout << "merge  rank: " << env.rank() << std::endl;
      //loser_tree.writeElementsToStream(out_, num_recv_elems);
      //    });
      std::vector<size_t> oldLcps;
      if (AllToAllStringPolicy::PrefixCompression) {
        loser_tree.writeElementsToStream(out_, num_recv_elems, oldLcps);
      }
      else {
        loser_tree.writeElementsToStream(out_, num_recv_elems);
      }
      StringLcpContainer<StringSet> sorted_string_cont;//(std::move(recv_string_cont));

      sorted_string_cont.set(std::move(recv_string_cont.raw_strings()));
      sorted_string_cont.set(std::move(sorted_string));
      sorted_string_cont.set(std::move(sorted_lcp));
      sorted_string_cont.setSavedLcps(std::move(oldLcps));

      return sorted_string_cont;
    }


  template<typename AllToAllStringPolicy, typename StringLcpContainer>
    static inline StringLcpContainer choose_merge(StringLcpContainer&& recv_string_cont,
        std::vector<std::pair<size_t, size_t>> ranges,
        size_t num_recv_elems,
        dsss::mpi::environment env = dsss::mpi::environment()) {

      switch (env.size()) {
        case 1 :  return merge<AllToAllStringPolicy,1>(std::move(recv_string_cont),
                      ranges,
                      num_recv_elems);
        case 2 : return merge<AllToAllStringPolicy,2>(std::move(recv_string_cont),
                     ranges,
                     num_recv_elems);
        case 4 : return merge<AllToAllStringPolicy,4>(std::move(recv_string_cont),
                     ranges,
                     num_recv_elems);
        case 8 : return merge<AllToAllStringPolicy,8>(std::move(recv_string_cont),
                     ranges,
                     num_recv_elems);
        case 16 : return merge<AllToAllStringPolicy,16>(std::move(recv_string_cont),
                      ranges,
                      num_recv_elems);
        case 32 : return merge<AllToAllStringPolicy,32>(std::move(recv_string_cont),
                      ranges,
                      num_recv_elems);
        case 64 : return merge<AllToAllStringPolicy,64>(std::move(recv_string_cont),
                      ranges,
                      num_recv_elems);
        case 128 : return merge<AllToAllStringPolicy,128>(std::move(recv_string_cont),
                       ranges,
                       num_recv_elems);
        case 264 : return merge<AllToAllStringPolicy,264>(std::move(recv_string_cont),
                       ranges,
                       num_recv_elems);
        case 512 : return merge<AllToAllStringPolicy,512>(std::move(recv_string_cont),
                       ranges,
                       num_recv_elems);
        default : std::cout << "Error in merge: K is not 2^i for i in {0,...,9} " << std::endl; 
                  std::abort();
      }
      return StringLcpContainer();
    }

  template <typename StringSet>
    class BloomFilter { 
      using String = typename StringSet::String;
      using Iterator = typename StringSet::Iterator;
      using CharIt = typename StringSet::CharIterator;
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

      std::vector<size_t> filter(dss_schimek::StringLcpPtr<StringSet> strptr, size_t depth, const std::vector<size_t>& candidates, std::vector<size_t>& results) {
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

        std::vector<HashTriple> recvHashTriples = dsss::mpi::AllToAllvSmall::alltoallv(hashes.data(), intervalSizes);
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

      std::vector<size_t> filter(dss_schimek::StringLcpPtr<StringSet> strptr) {
        return {};
      }
    };

    

  template<typename StringPtr, typename SampleSplittersPolicy, typename AllToAllStringPolicy, typename Timer>
    class DistributedMergeSort : private SampleSplittersPolicy, private AllToAllStringPolicy
  {
    public:
      dss_schimek::StringLcpContainer<typename StringPtr::StringSet>
        sort(StringPtr& local_string_ptr,
            dss_schimek::StringLcpContainer<typename StringPtr::StringSet>&& local_string_container, 
            dss_schimek::Timer& timer,
            dsss::mpi::environment env = dsss::mpi::environment()) {

          constexpr bool debug = false;

          using StringSet = typename StringPtr::StringSet;
          using Char = typename StringSet::Char;
          const StringSet& ss = local_string_ptr.active();
          std::size_t local_n = ss.size();
          
          // sort locally
          timer.start("sort_locally");
          dss_schimek::radixsort_CE3(local_string_ptr, 0, 0);
          timer.end("sort_locally");
          
          
          std::vector<size_t> results(ss.size(), 0);
          std::vector<size_t> candidates(ss.size());
          std::iota(candidates.begin(), candidates.end(), 0);
          BloomFilter<StringSet> bloomFilter;
          if (env.rank() == 0)
            ss.print();
          for (size_t i = 1; i < 10; ++i) {
            env.barrier();
            candidates = bloomFilter.filter(local_string_ptr, i, candidates, results);
            if (env.rank() == 0) {
              std::cout << "results:" << std::endl;
              for (size_t j = 0; j < results.size(); ++j) {
                std::cout << j << " " << results[j] << std::endl;
              }
            }
          }



          // There is only one PE, hence there is no need for distributed sorting 
          if (env.size() == 1)
            return dss_schimek::StringLcpContainer<StringSet>(std::move(local_string_container));

          timer.start("sample_splitters");
          std::vector<Char> raw_splitters = SampleSplittersPolicy::sample_splitters(ss);
          timer.end("sample_splitters");

          asm volatile("" ::: "memory");
          timer.add("allgather_splitters_bytes_sent", raw_splitters.size());
          timer.start("allgather_splitters");
          std::vector<Char> splitters =
            dss_schimek::mpi::allgather_strings(raw_splitters, env);
          timer.end("allgather_splitters");
          asm volatile("" ::: "memory");
          timer.start("choose_splitters");
          dss_schimek::StringLcpContainer chosen_splitters_cont = choose_splitters(ss, splitters);
          timer.end("choose_splitters");


          const StringSet chosen_splitters_set(chosen_splitters_cont.strings(),
              chosen_splitters_cont.strings() + chosen_splitters_cont.size());

          timer.start("compute_interval_sizes");
          std::vector<std::size_t> interval_sizes = compute_interval_binary(ss, chosen_splitters_set);
          std::vector<std::size_t> receiving_interval_sizes = dsss::mpi::alltoall(interval_sizes);
          timer.end("compute_interval_sizes");
          //print_interval_sizes(interval_sizes, receiving_interval_sizes);

          dss_schimek::StringLcpContainer<StringSet> recv_string_cont; 
          if constexpr(std::is_same<Timer, dss_schimek::Timer>::value) {
            timer.start("all_to_all_strings");
            recv_string_cont = 
              AllToAllStringPolicy::alltoallv(local_string_container, interval_sizes, timer);
            timer.end("all_to_all_strings");
          } else {
            timer.start("all_to_all_strings");
            EmptyTimer emptyTimer;
            recv_string_cont = 
              AllToAllStringPolicy::alltoallv(local_string_container, interval_sizes, emptyTimer);
            timer.end("all_to_all_strings");
          }
          timer.add("num_received_chars", recv_string_cont.char_size() - recv_string_cont.size());
          
          size_t num_recv_elems = 
            std::accumulate(receiving_interval_sizes.begin(), receiving_interval_sizes.end(), 0);

          assert(num_recv_elems == recv_string_cont.size());

          timer.start("compute_ranges");
          std::vector<std::pair<size_t, size_t>> ranges = 
            compute_ranges_and_set_lcp_at_start_of_range(recv_string_cont, receiving_interval_sizes);
          timer.end("compute_ranges");

          timer.start("merge_ranges");
          auto sorted_container = choose_merge<AllToAllStringPolicy>(std::move(recv_string_cont), ranges, num_recv_elems);
          timer.end("merge_ranges");
          return sorted_container;
        }
  };

}