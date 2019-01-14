#include "util/random_string_generator.hpp"
#include "mpi/environment.hpp"

int main() {
  using namespace dss_schimek;
  dsss::mpi::environment env;

  DNRationGeneratorStringLcpContainer<UCharLengthStringSet> generator(19);
  std::vector<unsigned char> initialChar{'A', 'A', 0};

  //std::vector<unsigned char> curChar(initialChar);
  //for (size_t i = 0; i < 1000; ++i) {
  //  std::cout << "output: " << curChar.data();
  //  std::cout << std::endl;
  //  curChar = generator.nextChar(curChar, 65, 90);
  //}
 
  const size_t numStrings = 5;
  size_t genNumStrings = 0;
  std::vector<unsigned char> rawStrings;
  size_t stringLength = 0;
  std::tie(rawStrings, genNumStrings, stringLength) = generator.getRawStrings(numStrings, 40, 0.90);
  size_t offset = 0;
  dss_schimek::mpi::execute_in_order([&] () {
      std::cout << "rank: " << env.rank() << std::endl;
      for (size_t j = 0; j < genNumStrings; ++j) {
        std::cout << j << " " << rawStrings.data() + offset;
        std::cout << std::endl;
        offset += stringLength + 1;
        }
        });
  DNRatioGenerator<UCharLengthStringSet> gen(200,60, 0.2);
  UCharLengthStringSet ss = gen.make_string_set();
  dss_schimek::mpi::execute_in_order([&] () {
      std::cout << "rank: " << env.rank() << std::endl;
      ss.print();
      env.rank();
      });
  env.finalize();
}

