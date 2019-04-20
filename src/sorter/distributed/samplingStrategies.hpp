#pragma once

#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>

#include "mpi/allgather.hpp"
#include "mpi/environment.hpp"

namespace dss_schimek {

template <typename StringSet>
class SampleSplittersNumStringsPolicy {
public:
    static constexpr bool isIndexed = false;
    static std::string getName() { return "NumStrings"; }

    static std::vector<typename StringSet::Char> sample_splitters(const StringSet& ss,
        size_t maxLength, uint64_t samplingFactor,
        dss_schimek::mpi::environment env = dss_schimek::mpi::environment()) {

        using Char = typename StringSet::Char;
        using String = typename StringSet::String;

        const size_t local_num_strings = ss.size();
        const size_t nr_splitters =
            std::min<size_t>(std::log2(local_num_strings) * samplingFactor, local_num_strings);
        const double splitter_dist = static_cast<double>(local_num_strings) /
                                     static_cast<double>(nr_splitters + 1);
        std::vector<Char> raw_splitters;
        raw_splitters.reserve(nr_splitters * (maxLength + 1u));

        for (size_t i = 1; i <= nr_splitters; ++i) {
            const String splitter =
                ss[ss.begin() + static_cast<size_t>(i * splitter_dist)];
            const size_t splitterLength = ss.get_length(splitter);
            const size_t usedSplitterLength =
                splitterLength > (maxLength) ? (maxLength) : splitterLength;
            std::copy_n(ss.get_chars(splitter, 0), usedSplitterLength,
                std::back_inserter(raw_splitters));
            raw_splitters.push_back(0);
        }
        return raw_splitters;
    }
};

uint64_t getLocalOffset(uint64_t localStringSize) {
    dss_schimek::mpi::environment env;
    auto allOffsets = dss_schimek::mpi::allgather(localStringSize);
    uint64_t localOffset = 0;
    for (size_t i = 0; i < env.rank(); ++i)
        localOffset += allOffsets[i];
    return localOffset;
}

struct SampleIndices {
    std::vector<unsigned char> sample;
    std::vector<uint64_t> indices;
};

template <typename StringSet>
class SampleIndexedSplittersNumStringsPolicy {
public:
    static constexpr bool isIndexed = true;
    static std::string getName() { return "IndexedNumStrings"; }

    static SampleIndices sample_splitters(const StringSet& ss, size_t maxLength,
        uint64_t samplingFactor,
        dss_schimek::mpi::environment env = dss_schimek::mpi::environment()) {

        using Char = typename StringSet::Char;
        using String = typename StringSet::String;
        SampleIndices sampleIndices;

        uint64_t localOffset = getLocalOffset(ss.size());

        const size_t local_num_strings = ss.size();
        const size_t nr_splitters = std::min<uint64_t>(
            std::log2(local_num_strings) * samplingFactor, local_num_strings);
        const double splitter_dist = static_cast<double>(local_num_strings) /
                                     static_cast<double>(nr_splitters + 1);
        std::vector<Char>& raw_splitters = sampleIndices.sample;
        raw_splitters.reserve(nr_splitters * (maxLength + 1u));
        std::vector<uint64_t>& splitterIndices = sampleIndices.indices;
        splitterIndices.resize(nr_splitters, localOffset);

        for (size_t i = 1; i <= nr_splitters; ++i) {
            const uint64_t splitterIndex =
                static_cast<uint64_t>(i * splitter_dist);
            splitterIndices[i - 1] += splitterIndex;
            const String splitter = ss[ss.begin() + splitterIndex];
            const size_t splitterLength = ss.get_length(splitter);
            const size_t usedSplitterLength =
                splitterLength > (maxLength) ? (maxLength) : splitterLength;
            std::copy_n(ss.get_chars(splitter, 0), usedSplitterLength,
                std::back_inserter(raw_splitters));
            raw_splitters.push_back(0);
        }
        return sampleIndices;
    }
};

template <typename StringSet>
class SampleIndexedSplittersNumCharsPolicy {
public:
    static constexpr bool isIndexed = true;
    static std::string getName() { return "IndexedNumChars"; }

    static SampleIndices sample_splitters(const StringSet& ss, const size_t maxLength,
        const uint64_t samplingFactor,
        dss_schimek::mpi::environment env = dss_schimek::mpi::environment()) {

        using Char = typename StringSet::Char;
        using String = typename StringSet::String;

        SampleIndices sampleIndices;
        uint64_t localOffset = getLocalOffset(ss.size());
        const size_t num_chars =
            std::accumulate(ss.begin(), ss.end(), static_cast<size_t>(0u),
                [&ss](const size_t& sum, const String& str) {
                    return sum + ss.get_length(str);
                });

        const size_t local_num_strings = ss.size();
        const size_t nr_splitters = std::min<uint64_t>(
            std::log2(local_num_strings) * samplingFactor, local_num_strings);
        const size_t splitter_dist = num_chars / (nr_splitters + 1);

        std::vector<Char>& raw_splitters = sampleIndices.sample;
        raw_splitters.reserve(nr_splitters * (maxLength + 1u));
        std::vector<uint64_t>& splitterIndices = sampleIndices.indices;
        splitterIndices.resize(nr_splitters, localOffset);

        size_t string_index = 0;
        size_t i = 1;
        for (; i <= nr_splitters && string_index < local_num_strings; ++i) {
            size_t num_chars_seen = 0;
            while (num_chars_seen < splitter_dist &&
                   string_index < local_num_strings) {
                num_chars_seen += ss.get_length(ss[ss.begin() + string_index]);
                ++string_index;
            }
            splitterIndices[i - 1] += string_index - 1;

            const String splitter = ss[ss.begin() + string_index - 1];
            const size_t splitterLength = ss.get_length(splitter);
            const size_t usedSplitterLength =
                splitterLength > (maxLength) ? (maxLength) : splitterLength;
            std::copy_n(ss.get_chars(splitter, 0), usedSplitterLength,
                std::back_inserter(raw_splitters));
            raw_splitters.push_back(0);
        }
        if (i < nr_splitters) splitterIndices.resize(i);
        return sampleIndices;
    }
};

template <typename StringSet>
class SampleSplittersNumCharsPolicy {
public:
    static constexpr bool isIndexed = false;
    static std::string getName() { return "NumChars"; }

    static std::vector<typename StringSet::Char> sample_splitters(const StringSet& ss,
        size_t maxLength, uint64_t samplingFactor,
        dss_schimek::mpi::environment env = dss_schimek::mpi::environment()) {

        using Char = typename StringSet::Char;
        using String = typename StringSet::String;

        const size_t num_chars =
            std::accumulate(ss.begin(), ss.end(), static_cast<size_t>(0u),
                [&ss](const size_t& sum, const String& str) {
                    return sum + ss.get_length(str);
                });

        const size_t local_num_strings = ss.size();
        const size_t nr_splitters =
            std::min<size_t>(std::log2(local_num_strings) * samplingFactor, local_num_strings);
        const size_t splitter_dist = num_chars / (nr_splitters + 1);
        std::vector<Char> raw_splitters;
        raw_splitters.reserve(nr_splitters * (maxLength + 1));

        size_t string_index = 0;
        for (size_t i = 1; i <= nr_splitters && string_index < local_num_strings; ++i) {
            size_t num_chars_seen = 0;
            while (num_chars_seen < splitter_dist &&
                   string_index < local_num_strings) {
                num_chars_seen += ss.get_length(ss[ss.begin() + string_index]);
                ++string_index;
            }

            const String splitter = ss[ss.begin() + string_index - 1];
            const size_t splitterLength = ss.get_length(splitter);
            const size_t usedSplitterLength =
                splitterLength > (maxLength) ? (maxLength) : splitterLength;
            std::copy_n(ss.get_chars(splitter, 0), usedSplitterLength,
                std::back_inserter(raw_splitters));
            raw_splitters.push_back(0);
        }
        return raw_splitters;
    }
};
} // namespace dss_schimek
