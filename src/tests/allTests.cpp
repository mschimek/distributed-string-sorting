#include <iostream>
#include <sstream>
#include "tests/measuringTool_tests.hpp"
#include "tests/DNGenerator_tests.hpp"
#include "util/measuringTool.hpp"

int main() {
  using namespace dss_schimek::tests;
  dsss::mpi::environment env;
  std::cout << "start tests" << std::endl;
  measuringTool_basicTest1();
  DNRatioGenerator_test();
  std::cout << "end tests" << std::endl;
  env.finalize();
  return 0;
}


