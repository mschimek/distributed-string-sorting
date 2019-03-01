#pragma once

#include "strings/stringcontainer.hpp"

#include "merge/stringptr.hpp"
#include "merge/bingmann-lcp_losertree.hpp"

#include "mpi/environment.hpp"
#include "mpi/synchron.hpp"

#include <tlx/sort/strings/radix_sort.hpp>

namespace dss_schimek {

  template <typename StringSet>
    StringLcpContainer<StringSet> choose_splitters(
        const StringSet& ss, 
        std::vector<typename StringSet::Char>& all_splitters,
        dsss::mpi::environment env = dsss::mpi::environment())
    {
      using Char = typename StringSet::Char;
      using String = typename StringSet::String;

      StringLcpContainer<StringSet> all_splitters_cont(std::move(all_splitters));
      tlx::sort_strings_detail::StringLcpPtr all_splitters_strptr = all_splitters_cont.make_string_lcp_ptr();
      const StringSet& all_splitters_set = all_splitters_strptr.active();

      tlx::sort_strings_detail::radixsort_CI3(all_splitters_strptr, 0, 0);

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
      return ranges;
    }
}
