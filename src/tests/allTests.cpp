#include <iostream>
#include "tests/valueTracker_test.hpp"

int main() {
  dsss::mpi::environment env;
  std::cout << "start tests" << std::endl;
  using namespace tlx;
  valueTracker_test_1();
  valueTracker_test_2();
  std::cout << "end tests" << std::endl;
  env.finalize();
}

