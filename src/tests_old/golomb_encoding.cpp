#include <bitset>
#include "encoding/golomb_encoding.hpp"
#include <tlx/define/likely.hpp>
int main() {
  std::vector<size_t> values {6, 7, 7, 8, 9, 16, 222};
  values = getDeltaEncoding(values); 
  std::cout << "result before decoding, values's size: " << values.size() << std::endl;
  for (const size_t value : values)
    printBits(value);
  std::vector<size_t> decodedValues = getDeltaDecoding(values);
  for (const size_t value: decodedValues)
    std::cout << value << std::endl;
  decodedValues = getDecoding(values);
  for (const size_t value: decodedValues)
    std::cout << value << std::endl;

}
