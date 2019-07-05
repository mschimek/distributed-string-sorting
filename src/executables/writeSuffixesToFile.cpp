#include "mpi/warmup.hpp"
#include "sorter/distributed/lcp_aware_number_string_splitting.hpp"

#include "mpi/readInput.hpp"
#include "mpi/synchron.hpp"
#include "util/measuringTool.hpp"
#include "util/random_string_generator.hpp"
#include <map>

#include <tlx/cmdline_parser.hpp>


int main(std::int32_t argc, char const* argv[]) {
    using namespace dss_schimek;
    using StringSet = dss_schimek::UCharLengthStringSet;

    dss_schimek::mpi::environment env;
    env.barrier();

    std::string path = "";

    tlx::CmdlineParser cp;
    cp.set_description("a distributed sorter");
    cp.set_author("Matthias Schimek");
    cp.add_string('p', "path", path, " path to file");

    if (!cp.process(argc, argv)) {
        return -1;
    }

    auto cont = dss_schimek::SuffixGenerator<StringSet>(path);
    auto ss = cont.make_string_set();
    for (auto& str : ss)
      std::cout << ss.get_chars(str, 0) << std::endl;
    
    env.finalize();
}
