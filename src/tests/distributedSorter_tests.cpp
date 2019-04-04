#include "sorter/distributed/lcp_aware_number_string_splitting.hpp"

void test() {
    using namespace dss_schimek;
    using measurement::MeasuringTool;
    using StringSet = UCharStringSet;
    using MPIAllToAllRoutine = mpi::AllToAllvSmall;
    using ByteEncoder = EmptyLcpByteEncoderMemCpy;
    using SampleSplittersPolicy = SampleSplittersNumStringsPolicy<StringSet>;
    using AllToAllPolicy = dss_schimek::mpi::AllToAllStringImpl<true, StringSet,
        MPIAllToAllRoutine, ByteEncoder>;
    using StringLcpPtr =
        typename tlx::sort_strings_detail::StringLcpPtr<StringSet, size_t>;

    dss_schimek::mpi::environment env;
    MeasuringTool& measuringTool = MeasuringTool::measuringTool();
    measuringTool.setVerbose(false);

    std::vector<unsigned char> rawStrings = {'a', 'a', 0};
    StringLcpContainer<StringSet> container(std::move(rawStrings));
    auto stringPtr = container.make_string_lcp_ptr();

    DistributedMergeSort<StringLcpPtr, SampleSplittersPolicy, AllToAllPolicy>
        sorter;
    StringLcpContainer<StringSet> sorted_string_cont =
        sorter.sort(stringPtr, std::move(container));
    std::cout << "rank: " << env.rank() << " "
              << sorted_string_cont.make_string_set().size() << std::endl;

    std::cout << sorted_string_cont.savedLcps().size() << std::endl;
    if (AllToAllPolicy::PrefixCompression && env.size() > 1)
        sorted_string_cont.extendPrefix(sorted_string_cont.make_string_set(),
            sorted_string_cont.savedLcps());
}

int main() {
    dss_schimek::mpi::environment env;
    test();
    env.finalize();
}
