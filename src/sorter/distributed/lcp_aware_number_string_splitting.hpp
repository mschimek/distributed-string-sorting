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

template <typename StringSet>
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
    std::copy_n(splitter, dss_schimek::string_length(splitter) + 1,
        std::back_inserter(raw_splitters));
  }
  return raw_splitters; 
}



template<typename StringPtr>
static inline void merge_sort(StringPtr& local_string_ptr, dss_schimek::StringLcpContainer<unsigned char>& local_string_container, 
    dsss::mpi::environment env = dsss::mpi::environment()) {

  constexpr bool debug = false;

  using StringSet = typename StringPtr::StringSet;
  using Char = typename StringSet::Char;
  const StringSet& ss = local_string_ptr.active();
  std::size_t local_n = ss.size();

  // sort locally
  insertion_sort(local_string_ptr, 0, 0);

  // There is only one PE, hence there is no need for distributed sorting 
  if (env.size() == 1) {
    return;
  }

  if constexpr (debug) {
    if (env.rank() == 0) { std::cout << "Begin sampling" << std::endl; }
    env.barrier();
  }

  //auto nr_splitters = std::min<std::size_t>(env.size() - 1, local_n);
  //auto splitter_dist = local_n / (nr_splitters + 1);
  //std::vector<Char> raw_splitters;

  //for (std::size_t i = 1; i <= nr_splitters; ++i) {
  //  const auto splitter = ss[ss.begin() + i * splitter_dist];
  //  std::copy_n(splitter, dsss::string_length(splitter) + 1,
  //      std::back_inserter(raw_splitters));
  //}

  //if constexpr (debug) {
  //  if (env.rank() == 0) {
  //    std::cout << "nr_splitters " << nr_splitters << " splitter_dist: " << splitter_dist << std::endl;
  //    std::cout << raw_splitters.data() << std::endl;
  //  }
  //}

  std::vector<Char> raw_splitters = sample_splitters(ss);
  size_t nr_splitters = std::min<size_t>(env.size() - 1, ss.size());
  size_t splitter_dist = ss.size() / (nr_splitters + 1);
  // Gather all splitters and sort them to determine the final splitters
  std::vector<Char> splitters =
    dss_schimek::mpi::allgather_strings(raw_splitters, env);

  if constexpr (debug) {
    if (env.rank() == 0) { std::cout << "Received all splitters char_length: " << splitters.size() << std::endl; 
      dss_schimek::print(splitters);
    }

    env.barrier();
  }

  dss_schimek::StringLcpContainer<Char> splitter_container(std::move(splitters));
  StringSet splitter_set(splitter_container.strings(), splitter_container.strings() + splitter_container.size());
  dss_schimek::StringLcpPtr splitter_strptr(splitter_set, splitter_container.lcp_array());
  insertion_sort(splitter_strptr, 0, 0);

  if (debug) {
    if(env.rank() == 0) {
      splitter_set.print();
    }
  }


  nr_splitters = std::min<std::size_t>(env.size() - 1, splitter_set.size());
  splitter_dist = splitter_set.size() / (nr_splitters + 1);
  raw_splitters.clear();
  if (debug) {
    if (env.rank() == 0)
    {
      std::cout << "splitter_dist " << splitter_dist << std::endl;
    }
  }
  for (std::size_t i = 1; i <= nr_splitters; ++i) {
    const auto splitter = splitter_set[splitter_set.begin() + i * splitter_dist];
    std::copy_n(splitter, dss_schimek::string_length(splitter) + 1,
        std::back_inserter(raw_splitters));
  }

  if (debug){
    if (env.rank() == 0) {
      std::cout << " selected splitter " << std::endl;
      dss_schimek::print(raw_splitters);
    }
  }
  splitter_container = dss_schimek::StringLcpContainer(std::move(raw_splitters));


  if constexpr (debug) {
    if (env.rank() == 0) { std::cout << "global splitters " << splitter_container.size() << std::endl; 
      dss_schimek::print(splitter_container.raw_strings());
    }

    env.barrier();
  }

  // Now we need to split the local strings using the splitters
  // The size is given by the NUMBER of strings, not their lengths. The number
  // of characters that must be send to other PEs is computed in the alltoall-
  // function.
  std::vector<std::size_t> interval_sizes;
  std::size_t element_pos = 0;
  splitter_dist = local_n / (nr_splitters + 1);
  for (std::size_t i = 0; i < splitter_container.size(); ++i) {
    element_pos = (i + 1) * splitter_dist;
    if (debug) {
      if (env.rank() == 0) {
        std::cout << splitter_container.strings()[i] << std::endl;
      }
    }
    while(element_pos > 0 && !dss_schimek::leq(
          ss[ss.begin() + element_pos], splitter_container.strings()[i])) { --element_pos; }
    while (element_pos < local_n && dss_schimek::leq(
          ss[ss.begin() + element_pos], splitter_container.strings()[i])) { ++element_pos; }
    interval_sizes.emplace_back(element_pos);
  }
  interval_sizes.emplace_back(local_n);
  for (std::size_t i = interval_sizes.size() - 1; i > 0; --i) {
    interval_sizes[i] -= interval_sizes[i - 1];
  }
  if constexpr (debug) {
    if (env.rank() == 0) {
      std::cout << "interval sizes: " << std::endl;
      for (size_t i = 0; i < interval_sizes.size(); ++i) {
        std::cout << i << " " << interval_sizes[i] << std::endl;
      }
    }
  }

  constexpr bool print_interval_details = true;
  std::vector<std::size_t> receiving_sizes = dsss::mpi::alltoall(interval_sizes);
  if constexpr (print_interval_details) {
    for (std::int32_t rank = 0; rank < env.size(); ++rank) {
      if (env.rank() == rank) {
        std::size_t total_size = 0;
        std::cout << "### Sending interval sizes on PE " << rank << std::endl;
        for (const auto is : interval_sizes) {
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
        for (const auto is : receiving_sizes) {
          total_size += is;
          std::cout << is << ", ";
        }
        std::cout << "Total size: " << total_size << std::endl;
      }
      env.barrier();
    }
    if (env.rank() == 0) { std::cout << std::endl; }
  }

  if constexpr (debug) {
    if (env.rank() == 0) {
      std::cout << "Send strings to corresponding PEs" << std::endl;
    }
    env.barrier();
  }


  dss_schimek::StringLcpContainer<Char> recv_container = dsss::mpi::alltoallv(local_string_container, interval_sizes);
  env.barrier();
  if constexpr(debug) {
    if (env.rank() == 0) {
      dss_schimek::print(recv_container.raw_strings());
      std::cout << std::endl;
      for (size_t i = 0; i < recv_container.size(); ++i)
        std::cout << i << " " << recv_container.lcp_array()[i] << std::endl;
    }
  }

  //if(env.rank() == 0) {
  //  for (size_t i = 0; i < env.size(); ++i)
  //    std::cout << "receiving " << i << " " << receiving_sizes[i] << std::endl;
  //}

  std::vector<unsigned char> expected_on_PE_0 = {
    'A', '1', 0, 
    'A', '2', 0,
    'A', '3', 0,
    'A', '4', 0,
    'A', '5', 0,
    'B', '1', 0,
    'B', '2', 0
  };
  std::vector<unsigned char> expected_on_PE_1 = {
    'B', '3', 0,
    'B', '4', 0,
    'B', '5', 0,
    'C', '1', 0,
    'C', '2', 0,
  };
  std::vector<unsigned char> expected_on_PE_2 = {
    'C', '3', 0,
    'C', '4', 0,
    'C', '5', 0
  };
  std::vector<size_t> expected_lcp_on_PE_0 = { 0, 1, 1, 1, 1, 0, 1};
  std::vector<size_t> expected_lcp_on_PE_1 = { 1, 1, 1, 0, 1 };
  std::vector<size_t> expected_lcp_on_PE_2 = { 1, 1, 1 }; 
  //switch (env.rank()) {
  // case 0 : 
  //   if (expected_on_PE_0 != recv_container.raw_strings() ||
  //       expected_lcp_on_PE_0 != recv_container.lcps())
  //   {
  //     std::cout << "\t\tPE 0 did not receive correct data: " << std::endl;
  //     std::abort();
  //   }
  //   break;
  // case 1 : 
  //   if (expected_on_PE_1 != recv_container.raw_strings() ||
  //       expected_lcp_on_PE_1 != recv_container.lcps())
  //   {
  //     std::cout << "\t\tPE 1 did not receive correct data: " << std::endl;
  //     std::abort();
  //   }
  //   break;
  // case 2 : 
  //   if (expected_on_PE_2 != recv_container.raw_strings() ||
  //       expected_lcp_on_PE_2 != recv_container.lcps())
  //   {
  //     std::cout << "\t\tPE 2 did not receive correct data: " << std::endl;
  //     std::abort();
  //   }
  //   break;
  //}

  std::vector<std::pair<size_t, size_t>> pairs;
  for(size_t i = 0, offset = 0; i < env.size(); ++i)
  {
    if (receiving_sizes[i] == 0)
    {
      pairs.emplace_back(0, 0);
      continue;
    }
    *(recv_container.lcp_array() + offset) = 0;
    pairs.emplace_back(offset, receiving_sizes[i]);
    offset += receiving_sizes[i];

  }

  dss_schimek::mpi::execute_in_order([&] () {
      std::cout << "rank: " << env.rank() << " pairs:" << std::endl;
      for(size_t i = 0; i < env.size(); ++i)
      std::cout << i << " " << pairs[i].first << " " << pairs[i].second << std::endl;

      });
  int num_recv_elems = std::accumulate(receiving_sizes.begin(), receiving_sizes.end(), 0);
 
  for (size_t i = 0; i < num_recv_elems; ++i)
    *(recv_container.lcp_array() + i) = 0;

  if (!recv_container.is_consistent())
  {
    std::cout << "rank: " << env.rank() << " not consistent: " << std::endl;
    std::abort();
  }

  stringtools::LcpStringPtr total(recv_container.strings(), recv_container.lcp_array(), num_recv_elems);
  //dss_schimek::print_str(total.strings, num_recv_elems);
  const size_t KWAY = 16;


  env.barrier();
  bingmann::LcpStringLoserTree<KWAY> loser_tree_(total, pairs.data());
  //dss_schimek::print_str(total.strings, total.size);
  std::vector<Char*> sortedStrings_(num_recv_elems);
  std::vector<size_t> sortedLcp_(num_recv_elems);
  stringtools::LcpStringPtr out_(sortedStrings_.data(), sortedLcp_.data(), num_recv_elems);
  loser_tree_.writeElementsToStream(out_, num_recv_elems);
  asm volatile("": : : "memory");
  dss_schimek::UCharStringSet ss_res(sortedStrings_.data(), sortedStrings_.data() + num_recv_elems);
  StringPtr strptr_res(ss_res, sortedLcp_.data());
  const bool is_sorted = dss_schimek::is_sorted(strptr_res);

  std::cout << "res: " << is_sorted << std::endl;
  //std::cout << "print final result " << std::endl;
    /*
     std::vector<decltype(local_string_set.cbegin())> string_it(
     env.size(), local_string_set.cbegin());
     std::vector<decltype(local_string_set.cbegin())> end_it(
     env.size(), local_string_set.cbegin() + receiving_sizes[0]);

     for (std::int32_t i = 1; i < env.size(); ++i) {
     string_it[i] = string_it[i - 1] + receiving_sizes[i - 1];
     end_it[i] = end_it[i - 1] + receiving_sizes[i];
     }

     if constexpr (debug) {
     if (env.rank() == 0) {
     std::cout << "Sort received strings" << std::endl;
     }
     env.barrier();
     }

     struct string_compare {
     bool operator ()(const dsss::string a, const dsss::string b) {
     return (dsss::string_cmp(a, b) < 0);
     }
     }; // struct string_compare

     tlx::LoserTreeCopy<false, dsss::string, string_compare> lt(env.size());

     std::size_t filled_sources = 0;
     for (std::int32_t i = 0; i < env.size(); ++i) {
     if (string_it[i] == end_it[i]) { lt.insert_start(nullptr, i, true); }
     else {
     lt.insert_start(&*string_it[i], i, false);
     ++filled_sources;
     }
     }

     lt.init();

     std::vector<dsss::string> result;
     result.reserve(local_string_set.size());
     while (filled_sources) {
     std::int32_t source = lt.min_source();
     result.push_back(*string_it[source]);
     ++string_it[source];
     if (string_it[source] != end_it[source]) {
     lt.delete_min_insert(&*string_it[source], false);
     } else {
     lt.delete_min_insert(nullptr, true);
     --filled_sources;
     }
     }
     local_string_set.update(std::move(result));

     if constexpr (debug) {
     if (env.rank() == 0) {
     std::cout << "Finished" << std::endl;
     }
     env.barrier();
     }*/
}
