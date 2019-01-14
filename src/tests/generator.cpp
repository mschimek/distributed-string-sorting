#include "util/random_string_generator.hpp"

int main() {
  using namespace dss_schimek;

  DNRationGeneratorStringLcpContainer<UCharLengthStringSet> generator(19);
  std::vector<unsigned char> initialChar{'A', 'A', 0};

  //std::vector<unsigned char> curChar(initialChar);
  //for (size_t i = 0; i < 1000; ++i) {
  //  std::cout << "output: " << curChar.data();
  //  std::cout << std::endl;
  //  curChar = generator.nextChar(curChar, 65, 90);
  //}
 
  const size_t numStrings = 50000000;
  auto [rawStrings, stringLength] = generator.rawStrings(numStrings, 40, 0.90);
  size_t offset = 0;
  for (size_t j = 0; j < numStrings; ++j) {
    if (j % 50000 == 0) {
      std::cout << j << " " << rawStrings.data() + offset;
      std::cout << std::endl;
    }
    offset += stringLength + 1;
  }
}

