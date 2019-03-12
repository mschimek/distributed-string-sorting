#include <vector>
#include <iostream>
#include <iostream>

//#include "mpi/environment.hpp"

int main() {
  //dsss::mpi::environment env;
  const size_t gigabyte = 1000000000;
  const size_t megabyte = 1000000;
  std::vector<unsigned char*> vec(2019);
  std::cout << "start allocation: " << std::endl;
  for (auto& elem : vec)
    elem = new unsigned char[ 20u * gigabyte];
  std::cout << "allocation done " << std::endl;
  size_t i = 0;
  std::cout << "get char" << std::endl;
  std::cin >> i; 
  std::cout << "start writing" << std::endl;
  for (auto& elem : vec) {
    for (size_t i = 0; i < megabyte; ++i) {
      elem[i] = 20;
    }
  }
  std::cout << "writing done " << std::endl;
  for (volatile size_t i = 0; i < 1000000000; ++i) {

  }
  //env.finalize();
}
