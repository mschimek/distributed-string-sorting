#include "strings/stringptr.hpp"
#include "strings/stringset.hpp"
#include "strings/stringtools.hpp"

#include "sorter/local/strings/insertion_sort_unified.hpp"

#include "mpi/alltoall.hpp"
#include "mpi/allgather.hpp"

template<typename StringPtr>
static inline void merge_sort(StringPtr& local_string_ptr, dss_schimek::StringLcpContainer<unsigned char>& local_string_container, 
  dsss::mpi::environment env = dsss::mpi::environment()) {

  constexpr bool debug = true;

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

  auto nr_splitters = std::min<std::size_t>(env.size() - 1, local_n);
  auto splitter_dist = local_n / (nr_splitters + 1);
  std::vector<Char> raw_splitters;

  for (std::size_t i = 1; i <= nr_splitters; ++i) {
    const auto splitter = ss[ss.begin() + i * splitter_dist];
    std::copy_n(splitter, dsss::string_length(splitter) + 1,
      std::back_inserter(raw_splitters));
  }
  
  if constexpr (debug) {
    if (env.rank() == 0) {
      std::cout << "nr_splitters " << nr_splitters << " splitter_dist: " << splitter_dist << std::endl;
      std::cout << raw_splitters.data() << std::endl;
    }
  }

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
  

  dss_schimek::StringLcpContainer<Char> recv_containers = dsss::mpi::alltoallv(local_string_container, interval_sizes);
  env.barrier();
  for(volatile int i = 0; i < 10000000; ++i);
  if constexpr(debug) {
    if (env.rank() == 0) {
      dss_schimek::print(recv_containers.raw_strings());
      std::cout << std::endl;
      for (size_t i = 0; i < recv_containers.size(); ++i)
        std::cout << i << " " << recv_containers.lcp_array()[i] << std::endl;
    }
  }
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
/*
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
