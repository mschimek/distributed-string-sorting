#include <algorithm>
#include "strings/stringptr.hpp"
#include "strings/stringset.hpp"
#include "strings/stringtools.hpp"

#include "sorter/local/strings/insertion_sort_unified.hpp"

#include "mpi/alltoall.hpp"
#include "mpi/allgather.hpp"
#include "mpi/synchron.hpp"
#include "mpi/is_sorted.hpp"

#include "merge/stringtools.hpp"
#include "merge/bingmann-lcp_losertree.hpp"

namespace dss_schimek {

  static constexpr bool debug = false;

  template <typename StringSet>
    static inline std::vector<typename StringSet::Char> sample_splitters(const StringSet& ss,
        dsss::mpi::environment env = dsss::mpi::environment()) {

      using Char = typename StringSet::Char;
      using String = typename StringSet::String;

      const size_t local_num_strings = ss.size();
      const size_t nr_splitters = std::min<size_t>(env.size() - 1, local_num_strings);
      const size_t splitter_dist = local_num_strings / (nr_splitters + 1);
      std::vector<Char> raw_splitters;

      for (size_t i = 1; i <= nr_splitters; ++i) {
        const String splitter = ss[ss.begin() + i * splitter_dist];
        std::copy_n(splitter, dss_schimek::string_length(splitter) + 1,
            std::back_inserter(raw_splitters));
      }
      return raw_splitters; 
    }

  template <typename StringLcpContainer, typename StringSet>
    static inline StringLcpContainer choose_splitters(const StringSet& ss, 
        std::vector<typename StringSet::Char>& all_splitters,
        dsss::mpi::environment env = dsss::mpi::environment())
    {
      using Char = typename StringSet::Char;
      using String = typename StringSet::String;

      //TODO find a more convenient way to get a StrPtr out of an raw_string vector
      StringLcpContainer all_splitters_cont(std::move(all_splitters));
      const auto start = all_splitters_cont.strings();
      StringSet all_splitters_set(start, start + all_splitters_cont.size());
      dss_schimek::StringLcpPtr all_splitters_strptr(all_splitters_set,
          all_splitters_cont.lcp_array());

      insertion_sort(all_splitters_strptr, 0, 0);

      const size_t nr_splitters = std::min<std::size_t>(env.size() - 1, all_splitters_set.size());
      const size_t splitter_dist = all_splitters_set.size() / (nr_splitters + 1);

      std::vector<Char> raw_chosen_splitters;
      for (std::size_t i = 1; i <= nr_splitters; ++i) {
        const auto begin = all_splitters_set.begin();
        const String splitter = all_splitters_set[begin + i * splitter_dist];
        std::copy_n(splitter, dss_schimek::string_length(splitter) + 1,
            std::back_inserter(raw_chosen_splitters));
      }

      return StringLcpContainer(std::move(raw_chosen_splitters));
    }

  template <typename StringSet>
    static inline std::vector<size_t> compute_interval_sizes(const StringSet& ss, 
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
              ss[ss.begin() + element_pos], splitters[splitters.begin() + i])) { --element_pos; }
        while (element_pos < ss.size() && dss_schimek::leq(
              ss[ss.begin() + element_pos], splitters[splitters.begin() + i])) { ++element_pos; }
        interval_sizes.emplace_back(element_pos);
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

  template<typename StringPtr>
    static inline void merge_sort(StringPtr& local_string_ptr,
        dss_schimek::StringLcpContainer<unsigned char>& local_string_container, 
        dsss::mpi::environment env = dsss::mpi::environment()) {

      constexpr bool debug = false;

      using StringSet = typename StringPtr::StringSet;
      using Char = typename StringSet::Char;
      const StringSet& ss = local_string_ptr.active();
      std::size_t local_n = ss.size();

      // sort locally
      insertion_sort(local_string_ptr, 0, 0);

      // There is only one PE, hence there is no need for distributed sorting 
      if (env.size() == 1)
        return;


      if constexpr (debug) {
        if (env.rank() == 0) { std::cout << "Begin sampling" << std::endl; }
        env.barrier();
      }

      std::vector<Char> raw_splitters = sample_splitters(ss);
      size_t nr_splitters = std::min<size_t>(env.size() - 1, ss.size());
      size_t splitter_dist = ss.size() / (nr_splitters + 1);
      // Gather all splitters and sort them to determine the final splitters
      std::vector<Char> splitters =
        dss_schimek::mpi::allgather_strings(raw_splitters, env);


      dss_schimek::StringLcpContainer chosen_splitters_cont = 
        choose_splitters<dss_schimek::StringLcpContainer<unsigned char>>(ss, splitters);

      const StringSet chosen_splitters_set(chosen_splitters_cont.strings(),
          chosen_splitters_cont.strings() + chosen_splitters_cont.size());

      // Now we need to split the local strings using the splitters
      // The size is given by the NUMBER of strings, not their lengths. The number
      // of characters that must be send to other PEs is computed in the alltoall-
      // function.
      std::vector<std::size_t> interval_sizes = compute_interval_sizes(ss, chosen_splitters_set);
      std::vector<std::size_t> receiving_interval_sizes = dsss::mpi::alltoall(interval_sizes);
      print_interval_sizes(interval_sizes, receiving_interval_sizes);

      dss_schimek::StringLcpContainer<Char> recv_string_cont = 
        dsss::mpi::alltoallv(local_string_container, interval_sizes);

      size_t num_recv_elems = 
        std::accumulate(receiving_interval_sizes.begin(), receiving_interval_sizes.end(), 0);

      assert(num_recv_elems == recv_string_cont.size());

      std::vector<std::pair<size_t, size_t>> ranges = 
        compute_ranges_and_set_lcp_at_start_of_range(recv_string_cont, receiving_interval_sizes);


      stringtools::LcpStringPtr lt_all_strings(recv_string_cont.strings(),
          recv_string_cont.lcp_array(),
          recv_string_cont.size());

      //TODO refactor
      constexpr size_t KWAY = 16;
      bingmann::LcpStringLoserTree<KWAY> loser_tree_(lt_all_strings, ranges.data());
      std::vector<Char*> sortedStrings_(num_recv_elems);
      std::vector<size_t> sortedLcp_(num_recv_elems);
      stringtools::LcpStringPtr out_(sortedStrings_.data(), sortedLcp_.data(), num_recv_elems);
      loser_tree_.writeElementsToStream(out_, num_recv_elems);
      asm volatile("": : : "memory");
      dss_schimek::UCharStringSet ss_res(sortedStrings_.data(), sortedStrings_.data() + num_recv_elems);
      StringPtr strptr_res(ss_res, sortedLcp_.data());
      const bool is_sorted = dss_schimek::is_sorted(strptr_res);

      std::cout << "res: " << is_sorted << std::endl;
    }
}
