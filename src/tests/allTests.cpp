#include <iostream>
#include <sstream>
#include "tests/measuringTool_tests.hpp"
#include "util/measuringTool.hpp"

int main() {
  dsss::mpi::environment env;
  std::cout << "start tests" << std::endl;
  measuringTool_basicTest1();
  std::cout << "end tests" << std::endl;
  env.finalize();
}

