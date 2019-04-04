#include "mpi/readInput.hpp"
#include "mpi/allgather.hpp"

#include "strings/stringcontainer.hpp"
#include "strings/stringset.hpp"

#include <tlx/die/core.hpp>

int main() {
    using namespace dss_schimek;
    using StringSet = UCharLengthStringSet;
    dss_schimek::mpi::environment env;
    const std::string path = "../../../testData/smallTest.dat";

    std::vector<unsigned char> rawStrings = readFileAndDistribute(path);
    StringLcpContainer<StringSet> cont(std::move(rawStrings));
    StringSet ss = cont.make_string_set();
    size_t recvLines = ss.size();
    auto allRecvLines = dss_schimek::mpi::allgather(recvLines);
    const size_t sumRecvLines = std::accumulate(allRecvLines.begin(), allRecvLines.end(), 0ull);

    tlx_die_unless(sumRecvLines == 9);
    env.finalize();
}
