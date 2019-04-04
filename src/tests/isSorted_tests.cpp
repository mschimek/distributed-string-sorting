#include <algorithm>

#include "strings/stringcontainer.hpp"
#include "strings/stringset.hpp"

#include "mpi/environment.hpp"
#include "mpi/is_sorted.hpp"
#include <tlx/die/core.hpp>

namespace dss_schimek {
namespace tests {
void test_AllPEHaveSortedData() {
    using namespace dss_schimek;
    using StringSet = UCharLengthStringSet;
    dss_schimek::mpi::environment env;
    unsigned char firstChar = env.rank() + 1;
    std::vector<unsigned char> myData{firstChar, 'a', 'a', 0};
    StringLcpContainer<StringSet> container(std::move(myData));
    auto stringPtr = container.make_string_lcp_ptr();
    const bool isSorted = is_sorted(stringPtr);

    tlx_die_unless(isSorted);
}

void test_AllButOnePEHaveSortedData_1() {
    using namespace dss_schimek;
    using StringSet = UCharLengthStringSet;
    dss_schimek::mpi::environment env;
    unsigned char firstChar = env.rank() + 65;
    std::vector<unsigned char> myData{
        firstChar, 'a', 'a', 0, firstChar, 'b', 'b', 0};
    if (env.rank() == 0) myData.clear();
    StringLcpContainer<StringSet> container(std::move(myData));
    auto stringPtr = container.make_string_lcp_ptr();
    const bool isSorted = is_sorted(stringPtr);

    tlx_die_unless(isSorted);
}

void test_AllButOnePEHaveSortedData_2() {
    using namespace dss_schimek;
    using StringSet = UCharLengthStringSet;
    dss_schimek::mpi::environment env;
    unsigned char firstChar = env.rank() + 65;
    std::vector<unsigned char> myData{
        firstChar, 'a', 'a', 0, firstChar, 'b', 'b', 0};
    if (env.rank() == env.size() / 2) myData.clear();
    StringLcpContainer<StringSet> container(std::move(myData));
    auto stringPtr = container.make_string_lcp_ptr();
    const bool isSorted = is_sorted(stringPtr);

    tlx_die_unless(isSorted);
}

void test_AllButOnePEHaveSortedData_3() {
    using namespace dss_schimek;
    using StringSet = UCharLengthStringSet;
    dss_schimek::mpi::environment env;
    unsigned char firstChar = env.rank() + 65;
    std::vector<unsigned char> myData{
        firstChar, 'a', 'a', 0, firstChar, 'b', 'b', 0};
    if (env.rank() + 1 == env.size()) myData.clear();
    StringLcpContainer<StringSet> container(std::move(myData));
    auto stringPtr = container.make_string_lcp_ptr();
    const bool isSorted = is_sorted(stringPtr);

    tlx_die_unless(isSorted);
}

void test_AllButnPEHaveSortedData(size_t n) {
    using namespace dss_schimek;
    using StringSet = UCharLengthStringSet;
    dss_schimek::mpi::environment env;

    std::vector<size_t> peIndices(env.size());
    std::iota(peIndices.begin(), peIndices.end(), 0);
    std::random_shuffle(peIndices.begin(), peIndices.end());
    std::vector<size_t> emptyPEs(peIndices.begin(), peIndices.begin() + n);

    tlx_die_unless(emptyPEs.size() == n);

    unsigned char firstChar = env.rank() + 65;
    std::vector<unsigned char> myData{
        firstChar, 'a', 'a', 0, firstChar, 'b', 'b', 0};

    if (std::find(emptyPEs.begin(), emptyPEs.end(), env.rank()) !=
        peIndices.end())
        myData.clear();

    StringLcpContainer<StringSet> container(std::move(myData));
    auto stringPtr = container.make_string_lcp_ptr();
    const bool isSorted = is_sorted(stringPtr);

    tlx_die_unless(isSorted);
}

void test_onePELocallyUnsorted() {
    using namespace dss_schimek;
    using StringSet = UCharLengthStringSet;
    dss_schimek::mpi::environment env;

    unsigned char firstChar = env.rank() + 65;
    std::vector<unsigned char> myData{
        firstChar, 'a', 'a', 0, firstChar, 'b', 'b', 0};

    if (env.rank() == 0)
        myData = {firstChar, 'b', 'a', 0, firstChar, 'a', 'a', 0};

    StringLcpContainer<StringSet> container(std::move(myData));
    auto stringPtr = container.make_string_lcp_ptr();
    const bool isSorted = is_sorted(stringPtr);

    tlx_die_if(isSorted);
}

void test_nPEGloballyUnsorted(size_t n) {
    using namespace dss_schimek;
    using StringSet = UCharLengthStringSet;
    dss_schimek::mpi::environment env;

    std::vector<size_t> peIndices(env.size());
    std::iota(peIndices.begin(), peIndices.end(), 0);
    std::random_shuffle(peIndices.begin(), peIndices.end());
    std::vector<size_t> emptyPEs(peIndices.begin(), peIndices.begin() + n);

    tlx_die_unless(emptyPEs.size() == n);

    unsigned char firstChar = env.rank() + 65;
    if (std::find(emptyPEs.begin(), emptyPEs.end(), env.rank()) !=
        emptyPEs.end())
        firstChar = env.size() - env.rank() + 65;

    std::vector<unsigned char> myData{
        firstChar, 'a', 'a', 0, firstChar, 'b', 'b', 0};

    StringLcpContainer<StringSet> container(std::move(myData));
    auto stringPtr = container.make_string_lcp_ptr();
    const bool isSorted = is_sorted(stringPtr);

    tlx_die_if(isSorted);
}

} // namespace tests
} // namespace dss_schimek

int main() {
    using namespace dss_schimek;
    using namespace tests;
    dss_schimek::mpi::environment env;

    if (env.rank() == 0)
    std::cout << "start isSorted tests" << std::endl;

    test_AllPEHaveSortedData();
    test_AllButOnePEHaveSortedData_1();
    test_AllButOnePEHaveSortedData_2();
    test_AllButOnePEHaveSortedData_3();
    for (size_t i = 0; i < env.size(); ++i)
        test_AllButnPEHaveSortedData(i);
    test_onePELocallyUnsorted();

    // i needs to start with 1, otherwise all PEs have globally sorted data
    for (size_t i = 1; i < env.size(); ++i)
        test_nPEGloballyUnsorted(i);
    if (env.rank() == 0)
    std::cout << "isSorted tests completed successfully" << std::endl;
    env.finalize();
}
