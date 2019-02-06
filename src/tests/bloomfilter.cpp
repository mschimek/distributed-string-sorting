#include "sorter/distributed/bloomfilter.hpp"
#include <iostream>

int main() {
  using namespace dss_schimek;
  std::vector<HashPEIndex> hashPEIndices {{5,1}, {7,1}, {3,2}, {5,2}, {1,3}, {2, 3}}; 
  std::vector<size_t> intervalSizes{2, 2, 2};
  auto results = FindDuplicates::findDuplicates(hashPEIndices, intervalSizes);
  //for (size_t i = 0; i < results.size(); ++i) {
  //  std::cout << "rank: " << i << std::endl;
  //  for (size_t j = 0; j < results[i].size(); ++j) 
  //    std::cout << "\t\t" <<  j << " " << results[i][j] << std::endl;
  }
}
