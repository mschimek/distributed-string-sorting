#include "sorter/distributed/prefix_doubling.hpp"
#include "mpi/readInput.hpp"
#include "mpi/synchron.hpp"
#include "util/measuringTool.hpp"
#include "util/random_string_generator.hpp"
#include <map>

#include <tlx/cmdline_parser.hpp>

struct GeneratedStringsArgs {
    size_t numOfStrings = 0;
    size_t stringLength = 0;
    size_t minStringLength = 0;
    size_t maxStringLength = 0;
    double dToNRatio = 0.5;
    std::string path = "";
};

template <typename StringGenerator, typename StringSet>
StringGenerator getGeneratedStringContainer(const GeneratedStringsArgs& args) {
    if constexpr (std::is_same_v<StringGenerator,
                      dss_schimek::DNRatioGenerator<StringSet>>) {
        return StringGenerator(
            args.numOfStrings, args.stringLength, args.dToNRatio);
    }
    else if constexpr (std::is_same_v<StringGenerator,
                      dss_schimek::FileDistributer<StringSet>>) {
        return dss_schimek::FileDistributer<StringSet>(args.path);
    }
    else if constexpr (std::is_same_v<StringGenerator,
                      dss_schimek::SuffixGenerator<StringSet>>) {
        return dss_schimek::SuffixGenerator<StringSet>(args.path);
    }
    else {
        return StringGenerator(
            args.numOfStrings, args.minStringLength, args.maxStringLength);
    }
}

template <typename StringSet, typename StringGenerator,
    typename SampleSplittersPolicy, typename MPIAllToAllRoutine,
    typename ByteEncoder, typename GolombEncoding, bool compressLcps>
void execute_sorter(size_t numOfStrings, const bool check,
    const bool exhaustiveCheck, size_t iteration, const bool strongScaling,
    GeneratedStringsArgs genStringArgs,
    dss_schimek::mpi::environment env = dss_schimek::mpi::environment()) {
    using StringLcpPtr =
        typename tlx::sort_strings_detail::StringLcpPtr<StringSet, size_t>;
    using namespace dss_schimek;
    using dss_schimek::measurement::MeasuringTool;

    std::string prefix =
        std::string("RESULT") +
        " numberProcessors=" + std::to_string(env.size()) +
        " samplePolicy=" + SampleSplittersPolicy::getName() +
        " StringGenerator=" + StringGenerator::getName() +
        " dToNRatio=" + std::to_string(genStringArgs.dToNRatio) +
        " stringLength=" + std::to_string(genStringArgs.stringLength) +
        " MPIAllToAllRoutine=" + MPIAllToAllRoutine::getName() +
        " ByteEncoder=" + ByteEncoder::getName() +
        " GolombEncoding=" + GolombEncoding::getName() +
        " StringSet=" + StringSet::getName() +
        " iteration=" + std::to_string(iteration) +
        " size=" + std::to_string(numOfStrings) +
        " strongScaling=" + std::to_string(strongScaling);

    MeasuringTool& measuringTool = MeasuringTool::measuringTool();
    measuringTool.setPrefix(prefix);

    CheckerWithCompleteExchange<StringLcpPtr> checker;

    if (!strongScaling) genStringArgs.numOfStrings *= env.size();
    if (env.rank() == 0)
      std::cout << "string generation started " << std::endl;
    env.barrier();

    // std::cout << "rank: " << env.rank() <<  " generate strings" << std::endl;
    StringGenerator generatedContainer =
        getGeneratedStringContainer<StringGenerator, StringSet>(genStringArgs);
    if (check || exhaustiveCheck)
        checker.storeLocalInput(generatedContainer.raw_strings());
    // std::cout << "rank: " << env.rank() <<  " generate strings completed" <<
    // std::endl;
    StringLcpPtr rand_string_ptr = generatedContainer.make_string_lcp_ptr();
    // dss_schimek::mpi::execute_in_order([&]() {
    //    env.barrier();
    //    std::cout << "rank: " << env.rank() << std::endl;
    //    multikey_quicksort(rand_container.make_string_lcp_ptr(), 0, 0);
    //    rand_container.make_string_set().print();
    //    });
    const size_t numGeneratedChars = generatedContainer.char_size();
    env.barrier();
    if (env.rank() == 0)
      std::cout << "string generation completed " << std::endl;
    const size_t numGeneratedStrings = generatedContainer.size();

    measuringTool.start("sorting_overall");
    using AllToAllPolicy =
        dss_schimek::mpi::AllToAllStringImplPrefixDoubling<compressLcps,
            StringLcpPtr, MPIAllToAllRoutine>;

    DistributedPrefixDoublingSort<StringLcpPtr, SampleSplittersPolicy,
        AllToAllPolicy, GolombEncoding>
        prefixDoublingSorter;
    std::vector<StringIndexPEIndex> permutation = prefixDoublingSorter.sort(
        std::move(generatedContainer), rand_string_ptr);

    measuringTool.stop("sorting_overall");

    if (check || exhaustiveCheck) {
        if (env.size() > 1) {

            StringLcpContainer<StringSet> originalInput(
                checker.getLocalInput());
            auto originalInputStrPtr = originalInput.make_string_lcp_ptr();
            tlx::sort_strings_detail::radixsort_CI3(originalInputStrPtr, 0, 0);

            auto CompleteStringsCont =
                dss_schimek::mpi::getStrings(permutation.begin(),
                    permutation.end(), originalInput.make_string_set());

            auto completeStringSet = CompleteStringsCont.make_string_set();
            reorder(completeStringSet, permutation.begin(), permutation.end());
            const StringLcpPtr sorted_strptr =
                CompleteStringsCont.make_string_lcp_ptr();
            const bool is_complete_and_sorted =
                dss_schimek::is_complete_and_sorted(sorted_strptr,
                    numGeneratedChars, CompleteStringsCont.char_size(),
                    numGeneratedStrings, CompleteStringsCont.size());
            //
            if (!is_complete_and_sorted) {
                std::cout << "not sorted" << std::endl;
                std::abort();
            }

            if (exhaustiveCheck) {
                const bool isSorted = checker.check(sorted_strptr, false);
                if (!isSorted) {
                    std::cout << "not complete sorted" << std::endl;
                    std::abort();
                }
            }
        }
    }
    env.barrier();
    std::stringstream buffer;
    measuringTool.writeToStream(buffer);

    if (env.rank() == 0) {
        std::cout << buffer.str() << std::endl;
    }
    measuringTool.reset();
}

namespace PolicyEnums {
enum class GolombEncoding {
    noGolombEncoding = 0,
    sequentialGolombEncoding = 1,
    pipelinedGolombEncoding = 2
};
GolombEncoding getGolombEncoding(size_t i) {
    switch (i) {
    case 0:
        return GolombEncoding::noGolombEncoding;
    case 1:
        return GolombEncoding::sequentialGolombEncoding;
    case 2:
        return GolombEncoding::pipelinedGolombEncoding;
    default:
        std::abort();
    }
}
enum class StringSet { UCharLengthStringSet = 0, UCharStringSet = 1 };
StringSet getStringSet(size_t i) {
    switch (i) {
    case 0:
        return StringSet::UCharLengthStringSet;
    case 1:
        return StringSet::UCharStringSet;
    default:
        std::abort();
    }
}
enum class StringGenerator {
    skewedRandomStringLcpContainer = 0,
    DNRatioGenerator = 1,
    File = 2,
    Suffix = 3
};
StringGenerator getStringGenerator(size_t i) {
    switch (i) {
    case 0:
        return StringGenerator::skewedRandomStringLcpContainer;
    case 1:
        return StringGenerator::DNRatioGenerator;
    case 2:
        return StringGenerator::File;
    case 3:
        return StringGenerator::Suffix;
    default:
        std::abort();
    }
}
enum class SampleString { numStrings = 0, numChars = 1 };
SampleString getSampleString(size_t i) {
    switch (i) {
    case 0:
        return SampleString::numStrings;
    case 1:
        return SampleString::numChars;
    default:
        std::abort();
    }
}
enum class MPIRoutineAllToAll { small = 0, directMessages = 1, combined = 2 };
MPIRoutineAllToAll getMPIRoutineAllToAll(size_t i) {
    switch (i) {
    case 0:
        return MPIRoutineAllToAll::small;
    case 1:
        return MPIRoutineAllToAll::directMessages;
    case 2:
        return MPIRoutineAllToAll::combined;
    default:
        std::abort();
    }
}
enum class ByteEncoder {
    emptyByteEncoderCopy = 0,
    emptyByteEncoderMemCpy = 1,
    sequentialDelayedByteEncoder = 2,
    sequentialByteEncoder = 3,
    interleavedByteEncoder = 4,
    emptyLcpByteEncoderMemCpy = 5
};
ByteEncoder getByteEncoder(size_t i) {
    switch (i) {
    case 0:
        return ByteEncoder::emptyByteEncoderCopy;
    case 1:
        return ByteEncoder::emptyByteEncoderMemCpy;
    case 2:
        return ByteEncoder::sequentialDelayedByteEncoder;
    case 3:
        return ByteEncoder::sequentialByteEncoder;
    case 4:
        return ByteEncoder::interleavedByteEncoder;
    case 5:
        return ByteEncoder::emptyLcpByteEncoderMemCpy;
    default:
        std::cout << "Enum ByteEncoder not defined" << std::endl;
        std::abort();
    }
}
} // namespace PolicyEnums

namespace PolicyEnums {
struct CombinationKey {
    CombinationKey(StringSet stringSet, GolombEncoding golombEncoding,
        StringGenerator stringGenerator, SampleString sampleStringPolicy,
        MPIRoutineAllToAll mpiAllToAllRoutine, ByteEncoder byteEncoder,
        bool compressLcps)
        : stringSet_(stringSet), golombEncoding_(golombEncoding),
          stringGenerator_(stringGenerator),
          sampleStringPolicy_(sampleStringPolicy),
          mpiRoutineAllToAll_(mpiAllToAllRoutine), byteEncoder_(byteEncoder),
          compressLcps_(compressLcps) {}

    StringSet stringSet_;
    GolombEncoding golombEncoding_;
    StringGenerator stringGenerator_;
    SampleString sampleStringPolicy_;
    MPIRoutineAllToAll mpiRoutineAllToAll_;
    ByteEncoder byteEncoder_;
    bool compressLcps_;

    bool operator==(const CombinationKey& other) {
        return stringSet_ == other.stringSet_ &&
               sampleStringPolicy_ == other.sampleStringPolicy_ &&
               other.mpiRoutineAllToAll_ == mpiRoutineAllToAll_ &&
               other.byteEncoder_ == byteEncoder_;
    }
};
} // namespace PolicyEnums

using namespace dss_schimek;

// template<typename StringSet>
// void secondArg_(const PolicyEnums::CombinationKey& key) {
//  //using SampleSplittersPolicy = SampleSplittersNumStringsPolicy<StringSet>;
//  //using StringGenerator = SkewedRandomStringLcpContainer<StringSet>;
//  //using MPIAllToAllRoutine = dss_schimek::mpi::AllToAllvSmall;
//  //using ByteEncoder = dss_schimek::EmptyByteEncoderCopy;
//}

struct SorterArgs {
    size_t size;
    bool check;
    bool exhaustiveCheck;
    size_t iteration;
    bool strongScaling;
    GeneratedStringsArgs generatorArgs;
};

template <typename StringSet, typename StringGenerator, typename SampleString,
    typename MPIRoutineAllToAll, typename ByteEncoder, typename GolombEncoding,
    bool compressLcps>
void eighthArg(const PolicyEnums::CombinationKey& key, const SorterArgs& args) {
    // execute_sorter<StringSet,
    //               StringGenerator,
    //               SampleString,
    //               MPIRoutineAllToAll,
    //               ByteEncoder,
    //               Timer>(args.size, args.checkInput, args.iteration,
    //               args.strongScaling, args.generatorArgs);
    execute_sorter<StringSet, StringGenerator, SampleString, MPIRoutineAllToAll,
        ByteEncoder, GolombEncoding, compressLcps>(args.size, args.check,
        args.exhaustiveCheck, args.iteration, args.strongScaling,
        args.generatorArgs);
}

template <typename StringSet, typename StringGenerator, typename SampleString,
    typename MPIRoutineAllToAll, typename ByteEncoder, typename GolombEncoding>
void seventhArg(
    const PolicyEnums::CombinationKey& key, const SorterArgs& args) {
    switch (key.compressLcps_) {
    case true:
        eighthArg<StringSet, StringGenerator, SampleString, MPIRoutineAllToAll,
            ByteEncoder, GolombEncoding, true>(key, args);
        break;

    case false:
        eighthArg<StringSet, StringGenerator, SampleString, MPIRoutineAllToAll,
            ByteEncoder, GolombEncoding, false>(key, args);
        break;
    }
}
template <typename StringSet, typename StringGenerator, typename SampleString,
    typename MPIRoutineAllToAll, typename ByteEncoder>
void sixthArg(const PolicyEnums::CombinationKey& key, const SorterArgs& args) {
    switch (key.golombEncoding_) {
    case PolicyEnums::GolombEncoding::noGolombEncoding: {
        using GolombEncoding = dss_schimek::AllToAllHashesNaive;
        seventhArg<StringSet, StringGenerator, SampleString, MPIRoutineAllToAll,
            ByteEncoder, GolombEncoding>(key, args);
        break;
    }
    case PolicyEnums::GolombEncoding::sequentialGolombEncoding: {
        using GolombEncoding = dss_schimek::AllToAllHashesGolomb;
        seventhArg<StringSet, StringGenerator, SampleString, MPIRoutineAllToAll,
            ByteEncoder, GolombEncoding>(key, args);
        break;
    }
    case PolicyEnums::GolombEncoding::pipelinedGolombEncoding: {
        using GolombEncoding = dss_schimek::AllToAllHashValuesPipeline;
        seventhArg<StringSet, StringGenerator, SampleString, MPIRoutineAllToAll,
            ByteEncoder, GolombEncoding>(key, args);
        break;
    }
    }
}

template <typename StringSet, typename StringGenerator, typename SampleString,
    typename MPIRoutineAllToAll>
void fifthArg(const PolicyEnums::CombinationKey& key, const SorterArgs& args) {
    switch (key.byteEncoder_) {
    case PolicyEnums::ByteEncoder::emptyByteEncoderCopy: {
        using ByteEncoder = dss_schimek::EmptyByteEncoderCopy;
        sixthArg<StringSet, StringGenerator, SampleString, MPIRoutineAllToAll,
            ByteEncoder>(key, args);
        break;
    }
    case PolicyEnums::ByteEncoder::emptyByteEncoderMemCpy: {
        using ByteEncoder = dss_schimek::EmptyByteEncoderMemCpy;
        sixthArg<StringSet, StringGenerator, SampleString, MPIRoutineAllToAll,
            ByteEncoder>(key, args);
        break;
    }
    case PolicyEnums::ByteEncoder::sequentialDelayedByteEncoder: {
        using ByteEncoder = dss_schimek::SequentialDelayedByteEncoder;
        sixthArg<StringSet, StringGenerator, SampleString, MPIRoutineAllToAll,
            ByteEncoder>(key, args);
        break;
    }
    case PolicyEnums::ByteEncoder::sequentialByteEncoder: {
        using ByteEncoder = dss_schimek::SequentialByteEncoder;
        sixthArg<StringSet, StringGenerator, SampleString, MPIRoutineAllToAll,
            ByteEncoder>(key, args);
        break;
    }
    case PolicyEnums::ByteEncoder::interleavedByteEncoder: {
        using ByteEncoder = dss_schimek::InterleavedByteEncoder;
        sixthArg<StringSet, StringGenerator, SampleString, MPIRoutineAllToAll,
            ByteEncoder>(key, args);
        break;
    }
    case PolicyEnums::ByteEncoder::emptyLcpByteEncoderMemCpy: {
        using ByteEncoder = dss_schimek::EmptyLcpByteEncoderMemCpy;
        sixthArg<StringSet, StringGenerator, SampleString, MPIRoutineAllToAll,
            ByteEncoder>(key, args);
        break;
    }
    };
}
template <typename StringSet, typename StringGenerator, typename SampleString>
void fourthArg(const PolicyEnums::CombinationKey& key, const SorterArgs& args) {
    switch (key.mpiRoutineAllToAll_) {
    case PolicyEnums::MPIRoutineAllToAll::small: {
        using MPIRoutineAllToAll = dss_schimek::mpi::AllToAllvSmall;
        fifthArg<StringSet, StringGenerator, SampleString, MPIRoutineAllToAll>(
            key, args);
        break;
    }
    case PolicyEnums::MPIRoutineAllToAll::directMessages: {
        using MPIRoutineAllToAll = dss_schimek::mpi::AllToAllvDirectMessages;
        fifthArg<StringSet, StringGenerator, SampleString, MPIRoutineAllToAll>(
            key, args);
        break;
    }
    case PolicyEnums::MPIRoutineAllToAll::combined: {
        using MPIRoutineAllToAll = dss_schimek::mpi::AllToAllvCombined<
            dss_schimek::mpi::AllToAllvSmall>;
        fifthArg<StringSet, StringGenerator, SampleString, MPIRoutineAllToAll>(
            key, args);
        break;
    }
    }
}
template <typename StringSet, typename StringGenerator>
void thirdArg(const PolicyEnums::CombinationKey& key, const SorterArgs& args) {
    switch (key.sampleStringPolicy_) {
    case PolicyEnums::SampleString::numStrings: {
        using SampleString =
            dss_schimek::SampleSplittersNumStringsPolicy<StringSet>;
        fourthArg<StringSet, StringGenerator, SampleString>(key, args);
        break;
    }
    case PolicyEnums::SampleString::numChars: {
        using SampleString =
            dss_schimek::SampleSplittersNumCharsPolicy<StringSet>;
        fourthArg<StringSet, StringGenerator, SampleString>(key, args);
        break;
    }
    };
}
template <typename StringSet>
void secondArg(const PolicyEnums::CombinationKey& key, const SorterArgs& args) {
    switch (key.stringGenerator_) {
    case PolicyEnums::StringGenerator::skewedRandomStringLcpContainer: {
        using StringGenerator =
            dss_schimek::SkewedRandomStringLcpContainer<StringSet>;
        thirdArg<StringSet, StringGenerator>(key, args);
        break;
    }
    case PolicyEnums::StringGenerator::DNRatioGenerator: {
        using StringGenerator = dss_schimek::DNRatioGenerator<StringSet>;
        thirdArg<StringSet, StringGenerator>(key, args);
        break;
    }
    case PolicyEnums::StringGenerator::File: {
        using StringGenerator = dss_schimek::FileDistributer<StringSet>;
        thirdArg<StringSet, StringGenerator>(key, args);
        break;
    }
    case PolicyEnums::StringGenerator::Suffix: {
        using StringGenerator = dss_schimek::SuffixGenerator<StringSet>;
        thirdArg<StringSet, StringGenerator>(key, args);
        break;
    }
    };
}

void firstArg(const PolicyEnums::CombinationKey& key, const SorterArgs& args) {
    switch (key.stringSet_) {
    case PolicyEnums::StringSet::UCharLengthStringSet:
        secondArg<UCharLengthStringSet>(key, args);
        break;
    case PolicyEnums::StringSet::UCharStringSet:
        // secondArg<UCharStringSet>(key, args);
        break;
    };
}

int main(std::int32_t argc, char const* argv[]) {
    using namespace dss_schimek;

    dss_schimek::mpi::environment env;
    env.barrier();

    bool check = false;
    bool exhaustiveCheck = false;
    bool compressLcps = false;
    unsigned int generator = 0;
    bool strongScaling = false;
    unsigned int sampleStringsPolicy =
        static_cast<int>(PolicyEnums::SampleString::numStrings);
    unsigned int golombPolicy =
        static_cast<int>(PolicyEnums::GolombEncoding::noGolombEncoding);
    unsigned int byteEncoder =
        static_cast<int>(PolicyEnums::ByteEncoder::emptyByteEncoderMemCpy);
    unsigned int mpiRoutineAllToAll =
        static_cast<int>(PolicyEnums::MPIRoutineAllToAll::small);
    unsigned int numberOfStrings = 100000;
    unsigned int numberOfIterations = 5;
    unsigned int stringLength = 50;
    double dToNRatio = 0.5;
    std::string path = "";

    tlx::CmdlineParser cp;
    cp.set_description("a distributed sorter");
    cp.set_author("Matthias Schimek");
    cp.add_double('r', "dToNRatio", dToNRatio, "D/N ratio");
    cp.add_string('y', "path", path, " path to file");
    cp.add_unsigned(
        's', "size", numberOfStrings, " number of strings to be generated");
    cp.add_unsigned('x', "golombEncodingPolicy", golombPolicy,
        " strategy for encoding of hashes sent to filter: noGolombEncoding = "
        "0, sequentialGolombEncoding = 1, pipelinedGolombEncoding = 2");
    cp.add_unsigned('p', "sampleStringsPolicy", sampleStringsPolicy,
        "0 = NumStrings, 1 = NumChars");
    cp.add_unsigned('b', "byteEncoder", byteEncoder,
        "emptyByteEncoderCopy = 0, emptyByteEncoderMemCpy = 1, "
        "sequentialDelayedByteEncoder = 2, sequentialByteEncoder = 3, "
        "interleavedByteEncoder = 4, emptyLcpByteEncoderMemCpy = 5");
    cp.add_unsigned('m', "MPIRoutineAllToAll", mpiRoutineAllToAll,
        "small = 0, directMessages = 1, combined = 2");
    cp.add_unsigned('i', "numberOfIterations", numberOfIterations, "");
    cp.add_flag('c', "check", check, " ");
    cp.add_flag('c', "exhaustiveCheck", exhaustiveCheck, " ");
    cp.add_unsigned('k', "generator", generator, " 0 = skewed, 1 = DNGen , 2 = file, 3 = suffix");
    cp.add_flag('x', "strongScaling", strongScaling, " ");
    cp.add_flag('v', "compressLcps", compressLcps,
        " compress Lcp values in alltoall-exchange");
    cp.add_unsigned('a', "stringLength", stringLength, " string Length ");

    if (!cp.process(argc, argv)) {
        return -1;
    }

    PolicyEnums::CombinationKey key(
        PolicyEnums::StringSet::UCharLengthStringSet,
        PolicyEnums::getGolombEncoding(golombPolicy),
        PolicyEnums::getStringGenerator(generator),
        PolicyEnums::getSampleString(sampleStringsPolicy),
        PolicyEnums::getMPIRoutineAllToAll(mpiRoutineAllToAll),
        PolicyEnums::getByteEncoder(byteEncoder), compressLcps);
    GeneratedStringsArgs generatorArgs;
    generatorArgs.numOfStrings = numberOfStrings;
    generatorArgs.stringLength = stringLength;
    generatorArgs.minStringLength = stringLength;
    generatorArgs.maxStringLength = stringLength + 10;
    generatorArgs.dToNRatio = dToNRatio;
    generatorArgs.path = path;
    for (size_t i = 0; i < numberOfIterations; ++i) {
        SorterArgs args = {numberOfStrings, check, exhaustiveCheck, i,
            strongScaling, generatorArgs};
        firstArg(key, args);
    }
    env.finalize();
}
