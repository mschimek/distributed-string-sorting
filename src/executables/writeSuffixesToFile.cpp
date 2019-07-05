#include "mpi/warmup.hpp"
#include "sorter/distributed/lcp_aware_number_string_splitting.hpp"

#include "mpi/readInput.hpp"
#include "mpi/synchron.hpp"
#include "util/measuringTool.hpp"
#include "util/random_string_generator.hpp"
#include <map>

#include <tlx/cmdline_parser.hpp>

namespace makeSuffixes {
std::vector<unsigned char> readFile(const std::string& path) {
        using dss_schimek::RawStringsLines;
        using namespace dss_schimek;
        using namespace dss_schimek::mpi;

        RawStringsLines data;
        const size_t fileSize = getFileSize(path);
        std::ifstream in(path);
        std::vector<unsigned char>& rawStrings = data.rawStrings;
        rawStrings.reserve(1.5 * fileSize);

        std::string line;
        data.lines = 0u;
        while (std::getline(in, line)) {
            ++data.lines;
            for (unsigned char curChar : line)
                rawStrings.push_back(curChar);
        }
        rawStrings.push_back(0);
        in.close();
        return rawStrings;
    }

}
int main(std::int32_t argc, char const* argv[]) {
    using namespace dss_schimek;

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

    const auto completeString = makeSuffixes::readFile(path);
    for (uint64_t i = 0; i + 1< completeString.size(); ++i) {
      std::cout << completeString.data() + i << std::endl;
    }
    
    env.finalize();
}
