
#include "mpi/environment.hpp"
#include "mpi/gather.hpp"
#include "util/random_string_generator.hpp"
#include <tlx/die/core.hpp>
#include <tlx/sort/strings/radix_sort.hpp>

namespace dss_schimek {
namespace tests {
std::string createFile() {
    using namespace std::chrono;
    const size_t systemTime =
        duration_cast<milliseconds>(system_clock::now().time_since_epoch())
            .count();
    std::string fileName = "input_" + std::to_string(systemTime);
    std::cout << "Filename: " << fileName << std::endl;
    std::ofstream outfile(fileName);
    outfile << "This is a very short sample text about nothing. Its only "
               "purpose is to create some characters for suffixes"
            << std::endl;
    outfile.close();
    return fileName;
}
void testSuffixes() {
    using namespace dss_schimek;
    using StringSet = UCharLengthStringSet;
    std::string fileName = "./" + createFile();
    SuffixGenerator<StringSet> suffixGenerator(fileName);
    StringSet ss = suffixGenerator.make_string_set();
    tlx::sort_strings_detail::radixsort_CI3(
        suffixGenerator.make_string_lcp_ptr(), 0, 0);
    // ss.print();
    std::ofstream outfile(fileName + "_sorted");
    for (size_t i = 0; i < ss.size(); ++i) {
        auto str = ss[ss.begin() + i];
        const char* chars = reinterpret_cast<const char*>(ss.get_chars(str, 0));
        std::string curString(chars);
        outfile << curString << std::endl;
    }
    outfile.close();
}
} // namespace tests
} // namespace dss_schimek

int main() {
    using namespace dss_schimek::tests;
    testSuffixes();
}
