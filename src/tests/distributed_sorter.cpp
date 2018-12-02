#include "sorter/distributed/lcp_aware_number_string_splitting.hpp"

#include "util/random_string_generator.hpp"
#include "util/timer.hpp"
#include "mpi/synchron.hpp"

#include <tlx/cmdline_parser.hpp>

template <typename StringSet, typename StringGenerator, typename SampleSplittersPolicy>
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
                       " StringSet=" + StringSet::getName() + 
                       " iteration=" + std::to_string(iteration) +
                       " size=" + std::to_string(numOfStrings);
                       
  Timer timer(prefix);

  StringGenerator rand_container(numOfStrings);
  StringLcpPtr rand_string_ptr = 
    rand_container.make_string_lcp_ptr();

  const size_t numGeneratedChars = rand_container.char_size();

  timer.start("sorting_overall");

  DistributedMergeSort<StringLcpPtr, SampleSplittersPolicy> sorter;
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

int main(std::int32_t argc, char const *argv[]) {
  using namespace dss_schimek;
  using StringSet = UCharLengthStringSet;

  dsss::mpi::environment env;
  enum SampleStringPolicy { NumStrings = 0, NumChars = 1};

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
  };
  env.finalize();
}
