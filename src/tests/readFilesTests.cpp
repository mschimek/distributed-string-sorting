#include "mpi/allgather.hpp"
#include "mpi/readInput.hpp"
#include "mpi/synchron.hpp"

#include "strings/stringcontainer.hpp"
#include "strings/stringset.hpp"

#include <tlx/die/core.hpp>
void testReadFileInParallel(const std::string& path) {
    using namespace dss_schimek;
    using StringSet = UCharLengthStringSet;
    mpi::environment env;

    std::vector<unsigned char> rawStrings = readFileInParallel(path);
    StringLcpContainer<StringSet> cont(std::move(rawStrings));
    StringSet ss = cont.make_string_set();
    size_t recvLines = ss.size();
    auto allRecvLines = dss_schimek::mpi::allgather(recvLines);
    const size_t sumRecvLines =
        std::accumulate(allRecvLines.begin(), allRecvLines.end(), 0ull);
    mpi::execute_in_order([&]() {
        std::cout << "rank: " << env.rank() << std::endl;
        ss.print();
    });

    tlx_die_unless(sumRecvLines == 9);
}

void testReadSequentially(const std::string& path) {

    using namespace dss_schimek;
    using StringSet = UCharLengthStringSet;

    std::vector<unsigned char> rawStrings = readFileAndDistribute(path);
    StringLcpContainer<StringSet> cont(std::move(rawStrings));
    StringSet ss = cont.make_string_set();
    size_t recvLines = ss.size();
    auto allRecvLines = dss_schimek::mpi::allgather(recvLines);
    const size_t sumRecvLines =
        std::accumulate(allRecvLines.begin(), allRecvLines.end(), 0ull);

    tlx_die_unless(sumRecvLines == 9);
}

int main() {
    using namespace dss_schimek;
    dss_schimek::mpi::environment env;
    const std::string path = "../../../testData/smallTest.dat";

    testReadFileInParallel(path);

    env.finalize();
}
