#include "encoding/integer_compression.hpp"
#include <bitset>
#include <iostream>
#include <tlx/die/core.hpp>
#include <vector>

namespace dss_schimek {
namespace tests {
std::vector<uint64_t> numbersToEncode() {
    const uint64_t oneByte = 1ull;
    const uint64_t twoBytes = 1ull << 8;
    const uint64_t threeBytes = 1ull << 16;
    const uint64_t fourBytes = 1ull << 24;
    const uint64_t fiveBytes = 1ull << 32;
    const uint64_t sixBytes = 1ull << 40;
    const uint64_t sevenBytes = 1ull << 48;
    const uint64_t eightBytes = 1ull << 56;
    return std::vector<uint64_t>{oneByte, twoBytes, threeBytes, fourBytes,
        fiveBytes, sixBytes, sevenBytes, eightBytes};
}
void integerCompressionTest() {
    using namespace dss_schimek;

    auto numbers = numbersToEncode();
    auto numbersRef = numbers;
    Writer writer(reinterpret_cast<unsigned char*>(numbers.data()));
    for (const uint64_t number : numbers) {
        writer.PutVarint(number);
    }

    std::vector<uint64_t> output(numbers.size());
    unsigned char* ch = reinterpret_cast<unsigned char*>(numbers.data());
    Reader reader(ch, output.begin(), output.end());
    reader.decode();

    tlx_die_unless(output == numbersRef);
}
} // namespace tests
} // namespace dss_schimek

int main() {
    using namespace dss_schimek::tests;
    std::cout << "start test" << std::endl;
    integerCompressionTest();
    std::cout << " completed all tests successfully" << std::endl;
}
