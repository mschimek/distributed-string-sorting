#include "sorter/distributed/lcp_aware_number_string_splitting.hpp"

#include "util/random_string_generator.hpp"
#include "util/timer.hpp"
#include "mpi/synchron.hpp"
#include <map>

#include <tlx/cmdline_parser.hpp>




template <typename StringSet, typename StringGenerator, 
         typename SampleSplittersPolicy, 
         typename MPIAllToAllRoutine, 
         typename ByteEncoder, 
         typename Timer>
           void execute_sorter(const size_t numOfStrings, const bool checkInput, size_t iteration,
               dsss::mpi::environment env = dsss::mpi::environment()) { 
             using StringLcpPtr = typename dss_schimek::StringLcpPtr<StringSet>;
             using namespace dss_schimek;

             std::string prefix = std::string("RESULT") +
               " numberProcessors=" + std::to_string(env.size()) +
               " samplePolicy=" + SampleSplittersPolicy::getName() +
               " StringGenerator=" + StringGenerator::getName() +
               " MPIAllToAllRoutine=" + MPIAllToAllRoutine::getName() + 
               " ByteEncoder=" + ByteEncoder::getName() + 
               " Timer=" + Timer::getName() + 
               " StringSet=" + StringSet::getName() + 
               " iteration=" + std::to_string(iteration) +
               " size=" + std::to_string(numOfStrings);

             dss_schimek::Timer timer(prefix);

             StringGenerator rand_container(numOfStrings);
             StringLcpPtr rand_string_ptr = 
               rand_container.make_string_lcp_ptr();

             const size_t numGeneratedChars = rand_container.char_size();

             timer.start("sorting_overall");
             using AllToAllPolicy = dss_schimek::mpi::AllToAllStringImpl<StringSet, MPIAllToAllRoutine, ByteEncoder, Timer>;
             DistributedMergeSort<StringLcpPtr, SampleSplittersPolicy, AllToAllPolicy, Timer> sorter;
             StringLcpContainer<StringSet> sorted_string_cont = 
               sorter.sort(rand_string_ptr, std::move(rand_container), timer);

             timer.end("sorting_overall");

             const StringLcpPtr sorted_strptr = sorted_string_cont.make_string_lcp_ptr();
             const bool is_complete_and_sorted = dss_schimek::is_complete_and_sorted(sorted_strptr,
                 numGeneratedChars,
                 sorted_string_cont.char_size(),
                 numOfStrings,
                 sorted_string_cont.size()); 
             
             if (!is_complete_and_sorted) {
               std::cout << "not sorted" << std::endl;
               std::abort(); 
             }
             std::stringstream buffer;
             timer.writeToStream(buffer);
             if (env.rank() == 0) {
               std::cout << buffer.str() << std::endl;
             }
           }

namespace PolicyEnums {
  enum class StringSet {UCharLengthStringSet = 0, UCharStringSet = 1};
  StringSet getStringSet(size_t i) {
    switch (i) {
      case 0: return StringSet::UCharLengthStringSet;
      case 1: return StringSet::UCharStringSet;
      default: std::abort();
    }
  }
  enum class StringGenerator {skewedRandomStringLcpContainer = 0};
  StringGenerator getStringGenerator(size_t i) {
    switch(i) {
      case 0 : return StringGenerator::skewedRandomStringLcpContainer;
      default: std::abort();
    }
  }
  enum class SampleString { numStrings = 0, numChars = 1};
  SampleString getSampleString(size_t i) {
    switch(i) {
      case 0 : return SampleString::numStrings;
      case 1 : return SampleString::numChars;
      default: std::abort();
    }
  }
  enum class MPIRoutineAllToAll { small = 0, directMessages = 1, combined = 2};
  MPIRoutineAllToAll getMPIRoutineAllToAll(size_t i) {
    switch(i) {
      case 0 : return MPIRoutineAllToAll::small;
      case 1 : return MPIRoutineAllToAll::directMessages;
      case 2 : return MPIRoutineAllToAll::combined;
      default : std::abort();
    }
  }
  enum class ByteEncoder { emptyByteEncoderCopy = 0, emptyByteEncoderMemCpy = 1, emptyByteEncoderMemCpyIntervalwise = 2, sequentialDelayedByteEncoder = 3, sequentialByteEncoder = 4, interleavedByteEncoder = 5 };
  ByteEncoder getByteEncoder(size_t i) {
    switch(i) {
      case 0 : return ByteEncoder::emptyByteEncoderCopy;
      case 1 : return ByteEncoder::emptyByteEncoderMemCpy;
      case 2 : return ByteEncoder::emptyByteEncoderMemCpyIntervalwise;
      case 3 : return ByteEncoder::sequentialDelayedByteEncoder;
      case 4 : return ByteEncoder::sequentialByteEncoder;
      case 5 : return ByteEncoder::interleavedByteEncoder;
      default : std::abort();
    }
  }
}

namespace PolicyEnums {
  struct CombinationKey  {
    CombinationKey(StringSet stringSet, 
        StringGenerator stringGenerator,
        SampleString sampleStringPolicy, 
        MPIRoutineAllToAll mpiAllToAllRoutine, 
        ByteEncoder byteEncoder) :
      stringSet_(stringSet), 
      stringGenerator_(stringGenerator), 
      sampleStringPolicy_(sampleStringPolicy), 
      mpiRoutineAllToAll_(mpiAllToAllRoutine), 
      byteEncoder_(byteEncoder) {}

    StringSet stringSet_;
    StringGenerator stringGenerator_;
    SampleString sampleStringPolicy_;
    MPIRoutineAllToAll mpiRoutineAllToAll_;
    ByteEncoder byteEncoder_;

    bool operator==(const CombinationKey& other) {
      return stringSet_ == other.stringSet_ && sampleStringPolicy_ == other.sampleStringPolicy_ 
        && other.mpiRoutineAllToAll_ == mpiRoutineAllToAll_ && other.byteEncoder_ == byteEncoder_;
    }
  };
}

using namespace dss_schimek;

template<typename StringSet>
void secondArg_(const PolicyEnums::CombinationKey& key) {
  using SampleSplittersPolicy = SampleSplittersNumStringsPolicy<StringSet>;
  using StringGenerator = SkewedRandomStringLcpContainer<StringSet>;
  using MPIAllToAllRoutine = dss_schimek::mpi::AllToAllvSmall;
  using ByteEncoder = dss_schimek::EmptyByteEncoderCopy;
}

struct SorterArgs {
  size_t size;
  bool checkInput;
  size_t iteration;
};

template<typename StringSet, typename StringGenerator, typename SampleString,
  typename MPIRoutineAllToAll, typename ByteEncoder>
   void sixthArg(const PolicyEnums::CombinationKey& key, const SorterArgs& args) {
   execute_sorter<StringSet,
                  StringGenerator,
                  SampleString,
                  MPIRoutineAllToAll,
                  ByteEncoder,
                  Timer>(args.size, args.checkInput, args.iteration);
   execute_sorter<StringSet,
                  StringGenerator,
                  SampleString,
                  MPIRoutineAllToAll,
                  ByteEncoder,
                  EmptyTimer>(args.size, args.checkInput, args.iteration);
   }

template<typename StringSet, typename StringGenerator, typename SampleString,
  typename MPIRoutineAllToAll>
  void fifthArg(const PolicyEnums::CombinationKey& key, const SorterArgs& args) {
    switch(key.byteEncoder_) {
      case PolicyEnums::ByteEncoder::emptyByteEncoderCopy : 
        {
          using ByteEncoder = dss_schimek::EmptyByteEncoderCopy;
          sixthArg<StringSet, StringGenerator, SampleString, MPIRoutineAllToAll, ByteEncoder>(key, args);
          break;
        }
      case PolicyEnums::ByteEncoder::emptyByteEncoderMemCpy : 
        {
          using ByteEncoder = dss_schimek::EmptyByteEncoderMemCpy;
          sixthArg<StringSet, StringGenerator, SampleString, MPIRoutineAllToAll, ByteEncoder>(key, args);
          break;
        }
      case PolicyEnums::ByteEncoder::emptyByteEncoderMemCpyIntervalwise : 
        {
          using ByteEncoder = dss_schimek::EmptyByteEncoderMemCpy;
          sixthArg<StringSet, StringGenerator, SampleString, MPIRoutineAllToAll, ByteEncoder>(key, args);
          break;
        }
      case PolicyEnums::ByteEncoder::sequentialDelayedByteEncoder : 
        {
          using ByteEncoder = dss_schimek::SequentialDelayedByteEncoder;
          sixthArg<StringSet, StringGenerator, SampleString, MPIRoutineAllToAll, ByteEncoder>(key, args);
          break;
        }
      case PolicyEnums::ByteEncoder::sequentialByteEncoder : 
        {
          using ByteEncoder = dss_schimek::SequentialByteEncoder;
          sixthArg<StringSet, StringGenerator, SampleString, MPIRoutineAllToAll, ByteEncoder>(key, args);
          break;
        }
      case PolicyEnums::ByteEncoder::interleavedByteEncoder : 
        {
          using ByteEncoder = dss_schimek::InterleavedByteEncoder;
          sixthArg<StringSet, StringGenerator, SampleString, MPIRoutineAllToAll, ByteEncoder>(key, args);
          break;
        }
    };
  }
template<typename StringSet, typename StringGenerator, typename SampleString>
void fourthArg(const PolicyEnums::CombinationKey& key, const SorterArgs& args) {
  switch(key.mpiRoutineAllToAll_) {
    case PolicyEnums::MPIRoutineAllToAll::small : 
      {
        using MPIRoutineAllToAll = dss_schimek::mpi::AllToAllvSmall;
        fifthArg<StringSet, StringGenerator, SampleString, MPIRoutineAllToAll>(key, args);        
        break;
      }
    case PolicyEnums::MPIRoutineAllToAll::directMessages : 
      {
        using MPIRoutineAllToAll = dss_schimek::mpi::AllToAllvDirectMessages;
        fifthArg<StringSet, StringGenerator, SampleString, MPIRoutineAllToAll>(key, args);        
        break;
      }
    case PolicyEnums::MPIRoutineAllToAll::combined : 
      {
        using MPIRoutineAllToAll = 
          dss_schimek::mpi::AllToAllvCombined<dss_schimek::mpi::AllToAllvSmall>;
        fifthArg<StringSet, StringGenerator, SampleString, MPIRoutineAllToAll>(key, args);        
        break;
      }
  }  
}
template<typename StringSet, typename StringGenerator>
void thirdArg(const PolicyEnums::CombinationKey& key, const SorterArgs& args) {
  switch(key.sampleStringPolicy_) {
    case PolicyEnums::SampleString::numStrings:
      {
        using SampleString = dss_schimek::SampleSplittersNumStringsPolicy<StringSet>;
        fourthArg<StringSet, StringGenerator, SampleString>(key, args); break;
      }
    case PolicyEnums::SampleString::numChars: 
      {
        using SampleString = dss_schimek::SampleSplittersNumCharsPolicy<StringSet>;
        fourthArg<StringSet, StringGenerator, SampleString>(key, args); break;
      }
  };
}
template<typename StringSet>
void secondArg(const PolicyEnums::CombinationKey& key, const SorterArgs& args) {
  switch(key.stringGenerator_) {
    case PolicyEnums::StringGenerator::skewedRandomStringLcpContainer : 
      using StringGenerator = dss_schimek::SkewedRandomStringLcpContainer<StringSet>; 
      thirdArg<StringSet, StringGenerator>(key, args); break;
  };
}

void firstArg(const PolicyEnums::CombinationKey& key, const SorterArgs& args) {
  switch (key.stringSet_) {
    case PolicyEnums::StringSet::UCharLengthStringSet : 
      secondArg<UCharLengthStringSet>(key, args); break;
    case PolicyEnums::StringSet::UCharStringSet : 
      //secondArg<UCharStringSet>(key, args); 
      break;
  };
}
 
int main(std::int32_t argc, char const *argv[]) {
  using namespace dss_schimek;
  using StringSet = UCharLengthStringSet;

  dsss::mpi::environment env;

  bool check = true;
  bool skewedInput = false;
  unsigned int sampleStringsPolicy = static_cast<int>(PolicyEnums::SampleString::numStrings);
  unsigned int byteEncoder = static_cast<int>(PolicyEnums::ByteEncoder::emptyByteEncoderCopy);
  unsigned int mpiRoutineAllToAll = static_cast<int>(PolicyEnums::MPIRoutineAllToAll::small);
  unsigned int numberOfStrings = 100000;
  unsigned int numberOfIterations = 5;

  tlx::CmdlineParser cp;
  cp.set_description("a distributed sorter");
  cp.set_author("Matthias Schimek");
  cp.add_unsigned('s', "size", numberOfStrings, " number of strings to be generated");
  cp.add_unsigned('p', "sampleStringsPolicy", sampleStringsPolicy, "0 = NumStrings, 1 = NumChars");
  cp.add_unsigned('b', "byteEncoder", byteEncoder, "emptyByteEncoder = 0, sequentialDelayedByteEncoder = 1, sequentialByteEncoder = 2, interleavedByteEncoder = 3");
  cp.add_unsigned('m', "MPIRoutineAllToAll", mpiRoutineAllToAll, "small = 0, directMessages = 1, combined = 2");
  cp.add_unsigned('i', "numberOfIterations", numberOfIterations, "");
  cp.add_flag('c', "checkSortedness", check, " ");
  cp.add_flag('k', "skewed", skewedInput, " ");

  if (!cp.process(argc, argv)) {
    return -1;
  }
  
  PolicyEnums::CombinationKey key(
      PolicyEnums::StringSet::UCharLengthStringSet, 
      PolicyEnums::StringGenerator::skewedRandomStringLcpContainer,
      PolicyEnums::getSampleString(sampleStringsPolicy),
      PolicyEnums::getMPIRoutineAllToAll(mpiRoutineAllToAll),
      PolicyEnums::getByteEncoder(byteEncoder));
  
  for (size_t i = 0; i < numberOfIterations; ++i) {
    SorterArgs args =  {numberOfStrings, check, i};
    firstArg(key, args);
  }
  env.finalize();
}
