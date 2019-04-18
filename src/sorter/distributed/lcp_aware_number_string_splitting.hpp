#include <algorithm>
#include <random>
#include <type_traits>

#include "strings/stringptr.hpp"
#include "strings/stringset.hpp"
#include "strings/stringtools.hpp"

#include "sorter/distributed/merging.hpp"
#include "sorter/distributed/misc.hpp"
#include "sorter/distributed/samplingStrategies.hpp"

#include "mpi/allgather.hpp"
#include "mpi/alltoall.hpp"
#include "mpi/byte_encoder.hpp"
#include "mpi/is_sorted.hpp"
#include "mpi/synchron.hpp"

#include "merge/bingmann-lcp_losertree.hpp"
#include "merge/stringtools.hpp"

#include "util/measuringTool.hpp"
#include <tlx/sort/strings/radix_sort.hpp>
#include <tlx/sort/strings/string_ptr.hpp>

namespace dss_schimek {

static constexpr bool debug = false;

//template <typename AllToAllStringPolicy, size_t K, typename StringSet>
//static inline StringLcpContainer<StringSet> merge(
//    dss_schimek::StringLcpContainer<StringSet>&& recv_string_cont,
//    const std::vector<std::pair<size_t, size_t>>& ranges,
//    const size_t num_recv_elems) {
//
//    std::vector<typename StringSet::String> sorted_string(
//        recv_string_cont.size());
//    std::vector<size_t> sorted_lcp(recv_string_cont.size());
//    StringSet ss = recv_string_cont.make_string_set();
//    // dss_schimek::mpi::execute_in_order([&] () {
//    //    dss_schimek::mpi::environment env;
//    //    env.barrier();
//    //    std::cout << "before merge: \n rank: " << env.rank() << std::endl;
//    //    ss.print();
//    //    for (size_t i = 0; i < ss.size(); ++i)
//    //    std::cout << i << " " << recv_string_cont.lcps()[i] << std::endl;
//    //    });
//    dss_schimek::StringLcpPtrMergeAdapter<StringSet> mergeAdapter(
//        ss, recv_string_cont.lcp_array());
//    dss_schimek::LcpStringLoserTree_<K, StringSet> loser_tree(
//        mergeAdapter, ranges.data());
//    StringSet sortedSet(
//        sorted_string.data(), sorted_string.data() + sorted_string.size());
//    dss_schimek::StringLcpPtrMergeAdapter out_(sortedSet, sorted_lcp.data());
//    // dss_schimek::mpi::execute_in_order([&] () {
//    //    dss_schimek::mpi::environment env;
//    //    env.barrier();
//
//    //    std::cout << "merge  rank: " << env.rank() << std::endl;
//    // loser_tree.writeElementsToStream(out_, num_recv_elems);
//    //    });
//    std::vector<size_t> oldLcps;
//    if (AllToAllStringPolicy::PrefixCompression) {
//        loser_tree.writeElementsToStream(out_, num_recv_elems, oldLcps);
//    }
//    else {
//        loser_tree.writeElementsToStream(out_, num_recv_elems);
//    }
//    StringLcpContainer<StringSet>
//        sorted_string_cont; //(std::move(recv_string_cont));
//
//    sorted_string_cont.set(std::move(recv_string_cont.raw_strings()));
//    sorted_string_cont.set(std::move(sorted_string));
//    sorted_string_cont.set(std::move(sorted_lcp));
//    sorted_string_cont.setSavedLcps(std::move(oldLcps));
//
//    return sorted_string_cont;
//}

//template <typename AllToAllStringPolicy, typename StringLcpContainer>
//static inline StringLcpContainer choose_merge(
//    StringLcpContainer&& recv_string_cont,
//    std::vector<std::pair<size_t, size_t>> ranges, size_t num_recv_elems,
//    dss_schimek::mpi::environment env = dss_schimek::mpi::environment()) {
//
//    switch (env.size()) {
//    case 1:
//        return merge<AllToAllStringPolicy, 1>(
//            std::move(recv_string_cont), ranges, num_recv_elems);
//    case 2:
//        return merge<AllToAllStringPolicy, 2>(
//            std::move(recv_string_cont), ranges, num_recv_elems);
//    case 4:
//        return merge<AllToAllStringPolicy, 4>(
//            std::move(recv_string_cont), ranges, num_recv_elems);
//    case 8:
//        return merge<AllToAllStringPolicy, 8>(
//            std::move(recv_string_cont), ranges, num_recv_elems);
//    case 16:
//        return merge<AllToAllStringPolicy, 16>(
//            std::move(recv_string_cont), ranges, num_recv_elems);
//    case 32:
//        return merge<AllToAllStringPolicy, 32>(
//            std::move(recv_string_cont), ranges, num_recv_elems);
//    case 64:
//        return merge<AllToAllStringPolicy, 64>(
//            std::move(recv_string_cont), ranges, num_recv_elems);
//    case 128:
//        return merge<AllToAllStringPolicy, 128>(
//            std::move(recv_string_cont), ranges, num_recv_elems);
//    case 264:
//        return merge<AllToAllStringPolicy, 264>(
//            std::move(recv_string_cont), ranges, num_recv_elems);
//    case 512:
//        return merge<AllToAllStringPolicy, 512>(
//            std::move(recv_string_cont), ranges, num_recv_elems);
//    default:
//        std::cout << "Error in merge: K is not 2^i for i in {0,...,9} "
//                  << std::endl;
//        std::abort();
//    }
//    return StringLcpContainer();
//}

template <typename StringPtr, typename SampleSplittersPolicy,
    typename AllToAllStringPolicy>
class DistributedMergeSort : private SampleSplittersPolicy,
                             private AllToAllStringPolicy {
public:
    dss_schimek::StringLcpContainer<typename StringPtr::StringSet> sort(
        StringPtr& local_string_ptr,
        dss_schimek::StringLcpContainer<typename StringPtr::StringSet>&&
            local_string_container,
        dss_schimek::mpi::environment env = dss_schimek::mpi::environment()) {

        // constexpr bool debug = false;
        using namespace dss_schimek;
        using namespace dss_schimek::measurement;

        using StringSet = typename StringPtr::StringSet;
        using Char = typename StringSet::Char;
        MeasuringTool& measuringTool = MeasuringTool::measuringTool();

        const StringSet& ss = local_string_ptr.active();
        const size_t lcpSummand = 25u;

        size_t charactersInSet = 0;
        for (const auto& str : ss) {
            charactersInSet += ss.get_length(str) + 1;
        }

        measuringTool.add(charactersInSet, "charactersInSet");

        // sort locally
        env.barrier();
        measuringTool.start("sort_locally");
        tlx::sort_strings_detail::radixsort_CI3(local_string_ptr, 0, 0);

        // dss_schimek::radixsort_CI3(local_string_ptr, 0, 0);
        measuringTool.stop("sort_locally");

        // There is only one PE, hence there is no need for distributed sorting
        if (env.size() == 1)
            return dss_schimek::StringLcpContainer<StringSet>(
                std::move(local_string_container));

        measuringTool.start("avg_lcp");
        const size_t globalLcpAvg = getAvgLcp(local_string_ptr) + lcpSummand;
        measuringTool.stop("avg_lcp");

        measuringTool.setPhase("splitter");
        measuringTool.start("sample_splitters");
        std::vector<Char> raw_splitters =
            SampleSplittersPolicy::sample_splitters(ss, globalLcpAvg);
        measuringTool.stop("sample_splitters");

        measuringTool.add(
            raw_splitters.size(), "allgather_splitters_bytes_sent");
        measuringTool.start("allgather_splitters");
        std::vector<Char> splitters =
            dss_schimek::mpi::allgather_strings(raw_splitters, env);
        measuringTool.stop("allgather_splitters");

        measuringTool.start("choose_splitters");
        dss_schimek::StringLcpContainer chosen_splitters_cont =
            choose_splitters(ss, splitters);
        measuringTool.stop("choose_splitters");

        const StringSet chosen_splitters_set(chosen_splitters_cont.strings(),
            chosen_splitters_cont.strings() + chosen_splitters_cont.size());

        measuringTool.start("compute_interval_sizes");
        std::vector<std::size_t> interval_sizes =
            compute_interval_binary(ss, chosen_splitters_set);
        std::vector<std::size_t> receiving_interval_sizes =
            dss_schimek::mpi::alltoall(interval_sizes);
        measuringTool.stop("compute_interval_sizes");

        measuringTool.setPhase("string_exchange");
        measuringTool.start("all_to_all_strings");
        dss_schimek::StringLcpContainer<StringSet> recv_string_cont =
            AllToAllStringPolicy::alltoallv(
                local_string_container, interval_sizes);
        measuringTool.stop("all_to_all_strings");

        measuringTool.add(
            recv_string_cont.char_size() - recv_string_cont.size(),
            "num_received_chars", false);
        measuringTool.add(
            recv_string_cont.size(),
            "num_recv_strings", false);
        measuringTool.setPhase("merging");

        size_t num_recv_elems =
            std::accumulate(receiving_interval_sizes.begin(),
                receiving_interval_sizes.end(), static_cast<size_t>(0u));

        assert(num_recv_elems == recv_string_cont.size());

        measuringTool.start("compute_ranges");
        std::vector<std::pair<size_t, size_t>> ranges =
            compute_ranges_and_set_lcp_at_start_of_range(
                recv_string_cont, receiving_interval_sizes);
        measuringTool.stop("compute_ranges");

        measuringTool.start("merge_ranges");
        auto sorted_container = choose_merge<AllToAllStringPolicy>(
            std::move(recv_string_cont), ranges, num_recv_elems);
        measuringTool.stop("merge_ranges");
        measuringTool.setPhase("none");
        return sorted_container;
    }
};

} // namespace dss_schimek
