#include "sorter/distributed/lcp_aware_number_string_splitting.hpp"

#include "util/random_string_generator.hpp"
#include "util/timer.hpp"
#include "mpi/synchron.hpp"
#include <map>

#include <tlx/cmdline_parser.hpp>

namespace dss_schimek::execution {



template <typename StringSet, typename StringGenerator, typename SampleSplittersPolicy, typename MPIAllToAllRoutine, typename ByteEncoder, typename Timer>
  void execute_sorter(const size_t numOfStrings, const bool checkInput,
      dsss::mpi::environment env = dsss::mpi::environment()) { 
    using StringLcpPtr = typename dss_schimek::StringLcpPtr<StringSet>;
    using namespace dss_schimek;
    static size_t iteration = 0;
    ++iteration;

    std::string prefix = std::string("RESULT") +
      " numberProcessors=" + std::to_string(env.size()) +
      " samplePolicy=" + SampleSplittersPolicy::getName() +
      " StringGenerator=" + StringGenerator::getName() +
      " MPIAllToAllRoutine=" + MPIAllToAllRoutine::getName() + 
      " ByteEncoder= " + ByteEncoder::getName() + 
      " Timer= " + Timer::getName() + 
      " StringSet=" + StringSet::getName() + 
      " iteration=" + std::to_string(iteration) +
      " size=" + std::to_string(numOfStrings);

    Timer timer(prefix);

    StringGenerator rand_container(numOfStrings);
    StringLcpPtr rand_string_ptr = 
      rand_container.make_string_lcp_ptr();

    const size_t numGeneratedChars = rand_container.char_size();

    timer.start("sorting_overall");
    using AllToAllPolicy = dss_schimek::mpi::AllToAllStringImpl<StringSet, dss_schimek::mpi::AllToAllvSmall, dss_schimek::EmptyByteEncoder, Timer>;
    DistributedMergeSort<StringLcpPtr, SampleSplittersPolicy, AllToAllPolicy> sorter;
    StringLcpContainer<StringSet> sorted_string_cont = 
      sorter.sort(rand_string_ptr, std::move(rand_container), timer);

    timer.end("sorting_overall");

    const StringLcpPtr sorted_strptr = sorted_string_cont.make_string_lcp_ptr();
    const bool is_complete_and_sorted = dss_schimek::is_complete_and_sorted(sorted_strptr,
        numGeneratedChars,
        sorted_string_cont.char_size(),
        numOfStrings,
        sorted_string_cont.size()); 

    std::stringstream buffer;
    timer.writeToStream(buffer);
    if (env.rank() == 0) {
      std::cout << buffer.str() << std::endl;
    }
  }
}
enum StringSet {UCharLengthStringSetE = 0, UCharStringSetE = 1};
enum SampleStringPolicys { NumStrings = 0, NumChars = 1};
enum MPIAllToAllRoutine { small = 0, directMessages = 1, combined = 2};
enum ByteEncoder { emptyByteEncoder = 0, sequentialDelayedByeEncoder = 1, sequentialByteEncoder = 2, interleavedByteEncoder = 3 };


template <typename T>
struct TypeContainer{
  using Type = T;
};


struct CombinationKey  {
  CombinationKey(StringSet stringSet, SampleStringPolicys sampleStringPolicy, MPIAllToAllRoutine mpiAllToAllRoutine, ByteEncoder byteEncoder) :
    stringSet_(stringSet), sampleStringPolicy_(sampleStringPolicy), mpiAllToAllRoutine_(mpiAllToAllRoutine), byteEncoder_(byteEncoder) {}
  StringSet stringSet_;
  SampleStringPolicys sampleStringPolicy_;
  MPIAllToAllRoutine mpiAllToAllRoutine_;
  ByteEncoder byteEncoder_;

  bool operator==(const CombinationKey& other) {
    return stringSet_ == other.stringSet_ && sampleStringPolicy_ == other.sampleStringPolicy_ 
      && other.mpiAllToAllRoutine_ == mpiAllToAllRoutine_ && other.byteEncoder_ == byteEncoder_;
  }
};


using namespace dss_schimek;
using namespace dss_schimek::execution;

template<typename StringSet>
void secondArg(const CombinationKey& key) {
  using SampleSplittersPolicy = SampleSplittersNumStringsPolicy<StringSet>;
  using StringGenerator = SkewedRandomStringLcpContainer<StringSet>;
  using MPIAllToAllRoutine = dss_schimek::mpi::AllToAllvSmall;
  using ByteEncoder = dss_schimek::EmptyByteEncoder;
  execute_sorter<StringSet, StringGenerator, SampleSplittersPolicy, MPIAllToAllRoutine, ByteEncoder, Timer>(1000000, false);
}

void firstArg(const CombinationKey& key) {
  switch (key.stringSet_) {
    case 0 : secondArg<UCharLengthStringSet>(key); break;
    case 1 : secondArg<UCharStringSet>(key); break;
  };
}
int main(std::int32_t argc, char const *argv[]) {
  using namespace dss_schimek;
  using StringSet = UCharLengthStringSet;

  dsss::mpi::environment env;

  bool check = false;
  bool skewedInput = false;
  unsigned int sampleStringsPolicy = NumStrings;
  unsigned int numberOfStrings = 10;
  unsigned int numberOfIterations = 5;

  tlx::CmdlineParser cp;
  cp.set_description("a distributed sorter");
  cp.set_author("Matthias Schimek");
  cp.add_unsigned('s', "size", numberOfStrings, " number of strings to be generated");
  cp.add_unsigned('p', "sampleStringsPolicy", sampleStringsPolicy, "0 = NumStrings, 1 = NumChars");
  cp.add_unsigned('i', "numberOfIterations", numberOfIterations, "");
  cp.add_flag('c', "checkSortedness", check, " ");
  cp.add_flag('k', "skewed", skewedInput, " ");
  
  if (!cp.process(argc, argv)) {
    return -1;
  }

  CombinationKey key =  {UCharLengthStringSetE, NumStrings, small, emptyByteEncoder};
  firstArg(key);
/*
  switch (sampleStringsPolicy) {
    case NumStrings : {
                        using SampleSplittersPolicy = SampleSplittersNumStringsPolicy<StringSet>;
                        if (skewedInput) {
                          using StringGenerator = SkewedRandomStringLcpContainer<StringSet>;
                          for (size_t i = 0; i < numberOfIterations; ++i)
                          execute_sorter<StringSet, StringGenerator, SampleSplittersPolicy>
                            (numberOfStrings, check);
                        } else {
                          using StringGenerator = RandomStringLcpContainer<StringSet>;
                          for (size_t i = 0; i < numberOfIterations; ++i)
                          execute_sorter<StringSet, StringGenerator, SampleSplittersPolicy>
                            (numberOfStrings, check);
                        }
                      } 
                      break;
    case NumChars : {
                      using SampleSplittersPolicy = SampleSplittersNumCharsPolicy<StringSet>;
                      if (skewedInput) {
                        using StringGenerator = SkewedRandomStringLcpContainer<StringSet>;
                          for (size_t i = 0; i < numberOfIterations; ++i)
                        execute_sorter<StringSet, StringGenerator, SampleSplittersPolicy>
                          (numberOfStrings, check);
                      } else {
                        using StringGenerator = RandomStringLcpContainer<StringSet>;
                          for (size_t i = 0; i < numberOfIterations; ++i)
                        execute_sorter<StringSet, StringGenerator, SampleSplittersPolicy>
                          (numberOfStrings, check);
                      }
                      break;
                    }
  };*/
  env.finalize();
}
