#include "mpi/environment.hpp"
#include "mpi/gather.hpp"
#include <tlx/die/core.hpp>

namespace dss_schimek {
namespace tests {
void gatherTest() {
    dsss::mpi::environment env;
    size_t root = 0;
    size_t data = env.rank();
    auto res = dss_schimek::mpi::gather(data, root);
    if (env.rank() != root) {
        tlx_die_unless(res.size() == 0);
    }
    else {
        std::vector<size_t> expectedData(env.size());
        std::iota(expectedData.begin(), expectedData.end(), 0u);

        tlx_die_unless(res.size() == env.size());
        tlx_die_unless(res == expectedData);
    }
}

void gathervTest() {
    dsss::mpi::environment env;
    size_t root = 0;
    size_t numberElems = 10000000ull;
    std::vector<size_t> data(numberElems, env.rank());
    auto res = dss_schimek::mpi::gatherv(data, root);
    if (env.rank() != root) {
        tlx_die_unless(res.size() == 0);
    }
    else {
        const size_t totalSum =
            std::accumulate(res.begin(), res.end(), static_cast<size_t>(0u));
        const size_t expectedSum =
            numberElems *
            (((env.size() - 1) * (env.size() - 1) + env.size() - 1)) / 2;
        tlx_die_unless(expectedSum == totalSum);
    }
}
} // namespace tests
} // namespace dss_schimek

int main() {
    using namespace dss_schimek;
    using namespace tests;
    using measurement::MeasuringTool;
    MeasuringTool& measuringTool = MeasuringTool::measuringTool();
    measuringTool.disable();

    dsss::mpi::environment env;
    std::cout << "start tests " << std::endl;
    gatherTest();
    gathervTest();
    std::cout << "completed tests successfully " << std::endl;
    env.finalize();
}

