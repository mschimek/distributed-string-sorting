#pragma once

#include "sorter/RQuick/RQuick.hpp"
#include "merge/stringtools.hpp"
#include "sorter/distributed/samplingStrategies.hpp"
#include "strings/stringcontainer.hpp"
#include <random>
//#include "merge/stringptr.hpp"
#include "merge/bingmann-lcp_losertree.hpp"

#include "mpi/allgather.hpp"
#include "mpi/environment.hpp"
#include "mpi/synchron.hpp"
#include "mpi/writeToFile.hpp"

#include <tlx/sort/strings/radix_sort.hpp>

namespace dss_schimek {

template <typename StringLcpPtr>
size_t getAvgLcp(const StringLcpPtr stringLcpPtr) {
    auto lcps = stringLcpPtr.get_lcp();
    struct LcpSumNumStrings {
        size_t lcpSum;
        size_t numStrings;
    };
    size_t localL = std::accumulate(
        lcps, lcps + stringLcpPtr.active().size(), static_cast<size_t>(0u));
    LcpSumNumStrings lcpSumNumStrings{localL, stringLcpPtr.active().size()};

    std::vector<LcpSumNumStrings> lcpSumsNumStrings =
        dss_schimek::mpi::allgather(lcpSumNumStrings);
    size_t totalL = 0;
    size_t totalNumString = 0;
    for (const auto& elem : lcpSumsNumStrings) {
        totalL += elem.lcpSum;
        totalNumString += elem.numStrings;
    }
    return totalL / totalNumString;
}

template <typename StringContainer>
std::vector<unsigned char> getSplitters(StringContainer& sortedLocalSample) {
    dss_schimek::mpi::environment env;
    uint64_t localSampleSize = sortedLocalSample.size();
    const auto allLocalSizes = dss_schimek::mpi::allgather(localSampleSize);
    const uint64_t localPrefix = std::accumulate(
        allLocalSizes.begin(), allLocalSizes.begin() + env.rank(), 0ull);
    const uint64_t totalSize = std::accumulate(
        allLocalSizes.begin() + env.rank(), allLocalSizes.end(), localPrefix);

    const size_t nr_splitters =
        std::min<std::size_t>(env.size() - 1, totalSize);
    const size_t splitter_dist = totalSize / (nr_splitters + 1);

    auto ss = sortedLocalSample.make_string_set();
    size_t splitterSize = 0u;
    for (std::size_t i = 1; i <= nr_splitters; ++i) {
        const uint64_t curIndex = i * splitter_dist;
        if (curIndex >= localPrefix &&
            curIndex < localPrefix + localSampleSize) {
            const auto str = ss[ss.begin() + curIndex - localPrefix];
            const auto length = ss.get_length(str) + 1;
            splitterSize += length;
        }
    }
    // std::cout << "rank: " << env.rank() << " localSample: " << ss.size()
    //          << " splitterSize: " << splitterSize
    //          << " nr_splitters: " << nr_splitters
    //          << " splitter_dist: " << splitter_dist
    //          << " localPrefix:" << localPrefix << " totalSize: " << totalSize
    //          << std::endl;

    std::vector<unsigned char> chosenSplitters(splitterSize);
    uint64_t curPos = 0;
    for (std::size_t i = 1; i <= nr_splitters; ++i) {
        const uint64_t curIndex = i * splitter_dist;
        if (curIndex >= localPrefix &&
            curIndex < localPrefix + localSampleSize) {
            const auto str = ss[ss.begin() + curIndex - localPrefix];
            const auto length = ss.get_length(str) + 1;
            auto chars = ss.get_chars(str, 0);
            std::copy_n(chars, length, chosenSplitters.data() + curPos);
            curPos += length;
        }
    }
    return dss_schimek::mpi::allgatherv(chosenSplitters);
}

template <typename StringLcpPtr>
std::vector<std::pair<uint64_t, uint64_t>> getDuplicateRanges(
    StringLcpPtr strptr) {
    using StartEnd = std::pair<uint64_t, uint64_t>;
    using String = typename StringLcpPtr::StringSet::String;
    std::vector<StartEnd> intervals;

    if (strptr.size() == 0) return intervals;

    mpi::environment env;

    auto ss = strptr.active();
    intervals.emplace_back(0, 0);
    uint64_t prevLength = ss.get_length(ss[ss.begin()]);
    for (size_t i = 1; i < strptr.size(); ++i) {
        const uint64_t curLcp = strptr.get_lcp(i);
        auto curString = ss[ss.begin() + i];
        const uint64_t curLength = ss.get_length(curString);
        // std::cout << "rank: " << env.rank() << " " << ss.get_chars(curString,
        // 0) <<  " " << strptr.get_lcp(i) << " length: " << curLength <<
        // std::endl;
        if (curLength != curLcp || prevLength != curLcp) {
            if (intervals.back().first + 1 != i) {
                intervals.back().second = i;
                intervals.emplace_back(i, i);
            }
            else {
                intervals.back().first = i;
            }
        }
        // std::cout << "rank: " << env.rank() << " " << intervals.back().first
        // << " " << intervals.back().first << std::endl; std::cout << "rank: "
        // << env.rank() << " " << ss.get_chars(curString, 0) <<  " " <<
        // strptr.get_lcp(i) << " length: " << curLength << std::endl;
        prevLength = curLength;
    }
    intervals.back().second = strptr.size();
    return intervals;
}

template <typename StringSet>
void sortRanges(dss_schimek::IndexStringLcpContainer<StringSet>& indexContainer,
    const std::vector<std::pair<uint64_t, uint64_t>>& ranges) {
    using IndexString = typename StringSet::String;
    for (auto [begin, end] : ranges) {
        std::sort(indexContainer.strings() + begin,
            indexContainer.strings() + end,
            [&](IndexString a, IndexString b) { return a.index < b.index; });
        ;
    }
}
template <typename StringSet>
IndexStringLcpContainer<StringSet> choose_splitters(
    IndexStringLcpContainer<StringSet>& indexContainer,
    dss_schimek::mpi::environment env = dss_schimek::mpi::environment()) {
    using Char = typename StringSet::Char;
    using String = typename StringSet::String;

    tlx::sort_strings_detail::StringLcpPtr all_splitters_strptr =
        indexContainer.make_string_lcp_ptr();
    const StringSet& all_splitters_set = all_splitters_strptr.active();

    tlx::sort_strings_detail::radixsort_CI3(all_splitters_strptr, 0, 0);
    auto ranges = getDuplicateRanges(all_splitters_strptr);
    sortRanges(indexContainer, ranges);

    const size_t nr_splitters =
        std::min<std::size_t>(env.size() - 1, all_splitters_set.size());
    const size_t splitter_dist = all_splitters_set.size() / (nr_splitters + 1);

    size_t splitterSize = 0u;
    for (std::size_t i = 1; i <= nr_splitters; ++i) {
        const auto begin = all_splitters_set.begin();
        const String splitter = all_splitters_set[begin + i * splitter_dist];
        splitterSize += all_splitters_set.get_length(splitter) + 1;
    }

    std::vector<Char> raw_chosen_splitters(splitterSize);
    std::vector<uint64_t> indices(splitterSize);
    size_t curPos = 0u;
    auto ss = all_splitters_set;

    for (std::size_t i = 1; i <= nr_splitters; ++i) {
        const auto begin = all_splitters_set.begin();
        const String splitter = all_splitters_set[begin + i * splitter_dist];
        indices[i - 1] = splitter.index;
        auto chars = ss.get_chars(splitter, 0);
        const size_t splitterLength = ss.get_length(splitter) + 1;
        std::copy(chars, chars + splitterLength,
            raw_chosen_splitters.begin() + curPos);
        curPos += splitterLength;
    }
    return IndexStringLcpContainer<StringSet>(
        std::move(raw_chosen_splitters), indices);
}
template <typename StringSet>
StringLcpContainer<StringSet> choose_splitters(const StringSet& ss,
    std::vector<typename StringSet::Char>& all_splitters,
    dss_schimek::mpi::environment env = dss_schimek::mpi::environment()) {
    using Char = typename StringSet::Char;
    using String = typename StringSet::String;

    StringLcpContainer<StringSet> all_splitters_cont(std::move(all_splitters));
    tlx::sort_strings_detail::StringLcpPtr all_splitters_strptr =
        all_splitters_cont.make_string_lcp_ptr();
    const StringSet& all_splitters_set = all_splitters_strptr.active();

    tlx::sort_strings_detail::radixsort_CI3(all_splitters_strptr, 0, 0);

    const size_t nr_splitters =
        std::min<std::size_t>(env.size() - 1, all_splitters_set.size());
    const size_t splitter_dist = all_splitters_set.size() / (nr_splitters + 1);

    size_t splitterSize = 0u;
    for (std::size_t i = 1; i <= nr_splitters; ++i) {
        const auto begin = all_splitters_set.begin();
        const String splitter = all_splitters_set[begin + i * splitter_dist];
        splitterSize += all_splitters_set.get_length(splitter) + 1;
    }

    std::vector<Char> raw_chosen_splitters(splitterSize);
    size_t curPos = 0u;

    for (std::size_t i = 1; i <= nr_splitters; ++i) {
        const auto begin = all_splitters_set.begin();
        const String splitter = all_splitters_set[begin + i * splitter_dist];
        auto chars = ss.get_chars(splitter, 0);
        const size_t splitterLength = ss.get_length(splitter) + 1;
        std::copy(chars, chars + splitterLength,
            raw_chosen_splitters.begin() + curPos);
        curPos += splitterLength;
    }
    return StringLcpContainer<StringSet>(std::move(raw_chosen_splitters));
}

template <typename StringSet>
inline std::vector<size_t> compute_interval_sizes(const StringSet& ss,
    const StringSet& splitters,
    dss_schimek::mpi::environment env = dss_schimek::mpi::environment()) {
    std::vector<size_t> interval_sizes;
    interval_sizes.reserve(splitters.size());

    size_t nr_splitters = std::min<size_t>(env.size() - 1, ss.size());
    size_t splitter_dist = ss.size() / (nr_splitters + 1);
    size_t element_pos = 0;

    for (std::size_t i = 0; i < splitters.size(); ++i) {
        element_pos = (i + 1) * splitter_dist;

        while (element_pos > 0 &&
               !dss_schimek::leq(ss.get_chars(ss[ss.begin() + element_pos], 0),
                   splitters.get_chars(splitters[splitters.begin() + i], 0))) {
            --element_pos;
        }

        while (element_pos < ss.size() &&
               dss_schimek::leq(ss.get_chars(ss[ss.begin() + element_pos], 0),
                   splitters.get_chars(splitters[splitters.begin() + i], 0))) {
            ++element_pos;
        }

        interval_sizes.emplace_back(element_pos);
    }
    interval_sizes.emplace_back(ss.size());
    for (std::size_t i = interval_sizes.size() - 1; i > 0; --i) {
        interval_sizes[i] -= interval_sizes[i - 1];
    }
    return interval_sizes;
}

template <typename StringSet>
inline static int binarySearch(
    const StringSet& ss, typename StringSet::CharIterator elem) {
    using String = typename StringSet::String;

    auto left = ss.begin();
    auto right = ss.end();

    while (left != right) {
        size_t dist = (right - left) / 2;
        String curStr = ss[left + dist];
        int res = dss_schimek::scmp(ss.get_chars(curStr, 0), elem);
        if (res < 0) {
            left = left + dist + 1;
        }
        else if (res == 0) {
            return left + dist - ss.begin();
        }
        else {
            right = left + dist;
        }
    }
    return left - ss.begin();
}

int indexStringCompare(const unsigned char* lhs, const uint64_t indexLhs,
    const unsigned char* rhs, const uint64_t indexRhs) {
    while (*lhs == *rhs && *lhs != 0) {
        ++lhs;
        ++rhs;
    } // TODO
    if (*lhs != *rhs) {
        return static_cast<int>(*lhs - *rhs);
    }
    return static_cast<int64_t>(indexLhs - indexRhs);
}

template <typename StringSet>
inline static int binarySearchIndexed(const StringSet& ss,
    const UCharLengthIndexStringSet& splitters, const uint64_t splitterIndex,
    const uint64_t localOffset) {
    using String = typename StringSet::String;

    auto left = ss.begin();
    auto right = ss.end();

    while (left != right) {
        size_t dist = (right - left) / 2;
        String curStr = ss[left + dist];
        uint64_t curIndex = left - ss.begin() + dist + localOffset;
        auto splitter = splitters[splitters.begin() + splitterIndex];
        int res = dss_schimek::indexStringCompare(ss.get_chars(curStr, 0),
            curIndex, splitters.get_chars(splitter, 0), splitter.index);
        if (res < 0) {
            left = left + dist + 1;
        }
        else if (res == 0) {
            return left + dist - ss.begin();
        }
        else {
            right = left + dist;
        }
    }
    return left - ss.begin();
}

template <typename StringSet> // TODO
inline std::vector<size_t> compute_interval_binary(const StringSet& ss,
    const StringSet& splitters,
    dss_schimek::mpi::environment env = dss_schimek::mpi::environment()) {
    using CharIt = typename StringSet::CharIterator;
    std::vector<size_t> interval_sizes;
    interval_sizes.reserve(splitters.size() + 1);

    for (std::size_t i = 0; i < splitters.size(); ++i) {
        CharIt splitter =
            splitters.get_chars(splitters[splitters.begin() + i], 0);
        size_t pos = binarySearch(ss, splitter);
        interval_sizes.emplace_back(pos);
    }
    interval_sizes.emplace_back(ss.size());
    for (std::size_t i = interval_sizes.size() - 1; i > 0; --i) {
        interval_sizes[i] -= interval_sizes[i - 1];
    }
    return interval_sizes;
}

template <typename StringSet> // TODO
inline std::vector<size_t> compute_interval_binary_index(const StringSet& ss,
    const UCharLengthIndexStringSet& splitters, const uint64_t localOffset,
    dss_schimek::mpi::environment env = dss_schimek::mpi::environment()) {
    using CharIt = typename StringSet::CharIterator;
    std::vector<size_t> interval_sizes;
    interval_sizes.reserve(splitters.size() + 1);

    for (std::size_t i = 0; i < splitters.size(); ++i) {
        size_t pos = binarySearchIndexed(ss, splitters, i, localOffset);
        interval_sizes.emplace_back(pos);
    }
    interval_sizes.emplace_back(ss.size());
    for (std::size_t i = interval_sizes.size() - 1; i > 0; --i) {
        interval_sizes[i] -= interval_sizes[i - 1];
    }
    return interval_sizes;
}

static inline void print_interval_sizes(
    const std::vector<size_t>& sent_interval_sizes,
    const std::vector<size_t>& recv_interval_sizes,
    dss_schimek::mpi::environment env = dss_schimek::mpi::environment()) {
    constexpr bool print_interval_details = true;
    if constexpr (print_interval_details) {
        for (std::uint32_t rank = 0; rank < env.size(); ++rank) {
            if (env.rank() == rank) {
                std::size_t total_size = 0;
                std::cout << "### Sending interval sizes on PE " << rank
                          << std::endl;
                for (const auto is : sent_interval_sizes) {
                    total_size += is;
                    std::cout << is << ", ";
                }
                std::cout << "Total size: " << total_size << std::endl;
            }
            env.barrier();
        }
        for (std::uint32_t rank = 0; rank < env.size(); ++rank) {
            if (env.rank() == rank) {
                std::size_t total_size = 0;
                std::cout << "### Receiving interval sizes on PE " << rank
                          << std::endl;
                for (const auto is : recv_interval_sizes) {
                    total_size += is;
                    std::cout << is << ", ";
                }
                std::cout << "Total size: " << total_size << std::endl;
            }
            env.barrier();
        }
        if (env.rank() == 0) {
            std::cout << std::endl;
        }
    }
}

template <typename StringLcpContainer>
static inline std::vector<std::pair<size_t, size_t>>
compute_ranges_and_set_lcp_at_start_of_range(
    StringLcpContainer& recv_string_cont,
    std::vector<size_t>& recv_interval_sizes,
    dss_schimek::mpi::environment env = dss_schimek::mpi::environment()) {
    std::vector<std::pair<size_t, size_t>> ranges;
    for (size_t i = 0, offset = 0; i < env.size(); ++i) {
        if (recv_interval_sizes[i] == 0) {
            ranges.emplace_back(0, 0);
            continue;
        }
        *(recv_string_cont.lcp_array() + offset) = 0;
        ranges.emplace_back(offset, recv_interval_sizes[i]);
        offset += recv_interval_sizes[i];
    }
    return ranges;
}

template <typename Sampler, typename StringPtr>
typename std::enable_if<!Sampler::isIndexed, std::vector<uint64_t>>::type
computePartition(
    StringPtr stringptr, uint64_t globalLcpAvg, uint64_t samplingFactor) {
    using StringSet = typename StringPtr::StringSet;
    using namespace dss_schimek;
    using measurement::MeasuringTool;

    auto ss = stringptr.active();
    MeasuringTool& measuringTool = MeasuringTool::measuringTool();
    measuringTool.setPhase("splitter");
    measuringTool.start("sample_splitters");
    std::vector<unsigned char> raw_splitters =
        Sampler::sample_splitters(ss, globalLcpAvg, samplingFactor);
    measuringTool.stop("sample_splitters");

    measuringTool.add(raw_splitters.size(), "allgather_splitters_bytes_sent");
    measuringTool.start("allgather_splitters");
    std::vector<unsigned char> splitters =
        dss_schimek::mpi::allgather_strings(raw_splitters);
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
    measuringTool.stop("compute_interval_sizes");
    return interval_sizes;
}

template <typename Sampler, typename StringPtr>
typename std::enable_if<Sampler::isIndexed, std::vector<uint64_t>>::type
computePartition(
    StringPtr stringptr, uint64_t globalLcpAvg, uint64_t samplingFactor) {
    using namespace dss_schimek;
    using namespace measurement;
    using IndexStringSet = UCharLengthIndexStringSet;
    MeasuringTool& measuringTool = MeasuringTool::measuringTool();

    auto ss = stringptr.active();
    measuringTool.setPhase("splitter");
    measuringTool.start("sample_splitters");
    auto sampleIndices = Sampler::sample_splitters(
        stringptr.active(), 2 * globalLcpAvg, samplingFactor);
    measuringTool.stop("sample_splitters");
    measuringTool.add(
        sampleIndices.sample.size(), "allgather_splitters_bytes_sent");
    measuringTool.start("allgather_splitters");
    auto recvSample = mpi::allgatherv(sampleIndices.sample);
    auto recvIndices = mpi::allgatherv(sampleIndices.indices);
    measuringTool.stop("allgather_splitters");

    measuringTool.start("choose_splitters");
    IndexStringLcpContainer<IndexStringSet> indexContainer(
        std::move(recvSample), recvIndices);
    indexContainer = choose_splitters(indexContainer);
    measuringTool.stop("choose_splitters");
    measuringTool.start("compute_interval_sizes");
    auto interval_sizes = compute_interval_binary_index(
        ss, indexContainer.make_string_set(), getLocalOffset(ss.size()));
    measuringTool.stop("compute_interval_sizes");
    return interval_sizes;
}

struct StringComparator {
    using String = dss_schimek::UCharLengthStringSet::String;
    bool operator()(String lhs, String rhs) {
        const unsigned char* lhsChars = lhs.string;
        const unsigned char* rhsChars = rhs.string;
        size_t counter = 0;
        // std::cout << "lhs: " << lhsChars << " rhs: " << rhsChars <<
        // std::endl;
        while (*lhsChars == *rhsChars && *lhsChars != 0) {
            ++lhsChars;
            ++rhsChars;
            counter++;
        }
        return *lhsChars < *rhsChars;
    }
};

template <typename Generator, typename Comparator>
dss_schimek::StringContainer<dss_schimek::UCharLengthStringSet> splitterSort(
    std::vector<unsigned char>&& rawStrings, Generator& generator,
    Comparator& comp) {
    dss_schimek::mpi::environment env;

    const bool isRobust = true;
    int tag = 11111;
    MPI_Comm comm = env.communicator();
    return RQuick::sort(
        generator, rawStrings, MPI_BYTE, tag, comm, comp, isRobust);
}

template <typename Sampler, typename StringPtr>
typename std::enable_if<!Sampler::isIndexed, std::vector<uint64_t>>::type
computePartition_(
    StringPtr stringptr, uint64_t globalLcpAvg, uint64_t samplingFactor) {
    using StringSet = typename StringPtr::StringSet;
    using StringContainer = dss_schimek::StringContainer<StringSet>;
    using namespace dss_schimek;
    using measurement::MeasuringTool;
    static int64_t iteration = -1;
    ++iteration;

    auto ss = stringptr.active();
    MeasuringTool& measuringTool = MeasuringTool::measuringTool();
    measuringTool.setPhase("splitter");
    measuringTool.start("sample_splitters");
    const uint64_t factor = 2;
    std::vector<unsigned char> raw_splitters =
        Sampler::sample_splitters(ss, 2 * globalLcpAvg, samplingFactor);
    measuringTool.stop("sample_splitters");

//    auto tmp1 = raw_splitters;
//    auto tmp2 = raw_splitters;
//    mpi::writeToOwnFile(
//        "TMP_Sample_iteration_" + std::to_string(iteration) + "_",
//        raw_splitters);
//    StringContainer tmp(std::move(tmp1));

    measuringTool.add(globalLcpAvg, "globalLcpAvg", false);

    dss_schimek::mpi::environment env;
    StringComparator comp;
    std::mt19937_64 generator;
    std::mt19937_64 generator2;
    int data_seed = 3469931 + env.rank();
    generator.seed(data_seed);
    generator2.seed(data_seed);
//    env.barrier();
//
//    measuringTool.start("createRBCComm");
////        RBC::Comm rbcComm;
////        RBC::Create_Comm_from_MPI(env.communicator(), &rbcComm);
//    measuringTool.stop("createRBCComm");
//    env.barrier();
//    measuringTool.start("firstRBCBarrier");
////    RBC::Barrier(rbcComm);
//    measuringTool.stop("firstRBCBarrier");
//    measuringTool.start("sort_splitterWarumup");
//    measuringTool.disable();
//    StringContainer sortedLocalSample2 =
//        splitterSort(std::move(tmp2), generator2, comp);
//    measuringTool.enable();
//    measuringTool.stop("sort_splitterWarumup");
//
    env.barrier();
    measuringTool.start("sort_splitter");
    StringContainer sortedLocalSample =
        splitterSort(std::move(raw_splitters), generator, comp);
    measuringTool.stop("sort_splitter");

    measuringTool.start("choose_splitters");
    auto rawChosenSplitters = getSplitters(sortedLocalSample);
    StringContainer chosen_splitters_cont(std::move(rawChosenSplitters));
    measuringTool.stop("choose_splitters");

    const StringSet chosen_splitters_set(chosen_splitters_cont.strings(),
        chosen_splitters_cont.strings() + chosen_splitters_cont.size());
    measuringTool.add(chosen_splitters_set.size(), "chosenSplitterSize", false);
    measuringTool.add(
        chosen_splitters_cont.char_size(), "chosenSplitterCharSize", false);

    measuringTool.start("compute_interval_sizes");
    std::vector<std::size_t> interval_sizes =
        compute_interval_binary(ss, chosen_splitters_set);
    measuringTool.stop("compute_interval_sizes");
    return interval_sizes;
}
} // namespace dss_schimek
