#pragma once

#include <vector>

#include "strings/stringcontainer.hpp"
#include "merge/stringptr.hpp"
#include "merge/bingmann-lcp_losertree.hpp"
#include "mpi/environment.hpp"

template<typename AllToAllStringPolicy, size_t K, typename StringSet>
    static inline dss_schimek::StringLcpContainer<StringSet> merge(
        dss_schimek::StringLcpContainer<StringSet>&& recv_string_cont,
        const std::vector<std::pair<size_t, size_t>>& ranges,
        const size_t num_recv_elems) {


      std::vector<typename StringSet::String> sorted_string(recv_string_cont.size());
      std::vector<size_t> sorted_lcp(recv_string_cont.size());
      StringSet ss = recv_string_cont.make_string_set();

      dss_schimek::StringLcpPtrMergeAdapter<StringSet> mergeAdapter(ss, recv_string_cont.lcp_array());
      dss_schimek::LcpStringLoserTree_<K, StringSet> loser_tree(mergeAdapter, ranges.data());
      StringSet sortedSet(sorted_string.data(), sorted_string.data() + sorted_string.size());
      dss_schimek::StringLcpPtrMergeAdapter out_(sortedSet, sorted_lcp.data());

      std::vector<size_t> oldLcps;
      if (AllToAllStringPolicy::PrefixCompression) {
        loser_tree.writeElementsToStream(out_, num_recv_elems, oldLcps);
      }
      else {
        loser_tree.writeElementsToStream(out_, num_recv_elems);
      }
      dss_schimek::StringLcpContainer<StringSet> sorted_string_cont;//(std::move(recv_string_cont));

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
