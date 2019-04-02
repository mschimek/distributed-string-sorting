#include "encoding/integer_compression.hpp"
#include "mpi/byte_encoder.hpp"
#include <algorithm>
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
    std::vector<uint8_t> compressedIntegers(numbers.size() * sizeof(uint64_t));
    Writer writer(compressedIntegers.begin());
    for (const uint64_t number : numbers) {
        writer.PutVarint(number);
    }

    std::vector<uint64_t> output(numbers.size());
    Reader reader(compressedIntegers.begin(), output.begin(), output.end());
    reader.decode();

    tlx_die_unless(output == numbersRef);
}

void integerCompressionTestRanges() {
    using namespace dss_schimek;

    auto numbers = numbersToEncode();
    std::vector<size_t> numbersNumbers(numbers);
    std::copy_n(
        numbers.begin(), numbers.size(), std::back_inserter(numbersNumbers));
    auto numbersNumbersRef = numbersNumbers;
    std::vector<size_t> counts(2, numbers.size());
    const auto& [compressedIntegers, compressedCounts] =
        IntegerCompression::writeRanges(
            counts.begin(), counts.end(), numbersNumbers.begin());

    const auto decompressedIntegers = IntegerCompression::readRanges(counts.begin(), counts.end(), compressedIntegers.begin());

    // ++++++++ test +++++++++
    const std::vector<uint64_t> expectedNumUsedBytes{37, 37};
    tlx_die_unless(expectedNumUsedBytes == compressedCounts);
    tlx_die_unless(numbersNumbersRef == decompressedIntegers);


}
} // namespace tests
} // namespace dss_schimek

int main() {
    using namespace dss_schimek::tests;
    std::cout << "start test" << std::endl;
    integerCompressionTest();
    integerCompressionTestRanges();
    std::cout << " completed all tests successfully" << std::endl;
}
