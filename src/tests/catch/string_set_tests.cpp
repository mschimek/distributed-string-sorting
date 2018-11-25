#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
#include <catch/catch.hpp>

#include <vector>
#include "strings/stringcontainer.hpp"


TEST_CASE( "check method in string sets work", "check method" ) {
  using namespace dss_schimek;
  std::vector<unsigned char> ordered_strings = {
    'A', 'A', 0,
    'B', 'A', 0,
    'B', 'A', 'B', 0,
    'C', 'C', 'C', 0,
    'C', 'C', 'D', 0
  };
  std::vector<unsigned char> unordered_strings_1 = {
    'A', 'A', 0,
    'B', 'A', 0,
    'B', 'A', 'B', 0,
    'C', 'C', 'C', 0,
    'C', 'C', 'A', 0
  };
  std::vector<unsigned char> unordered_strings_2 = {
    'A', 'A', 0,
    'B', 'B', 0,
    'B', 'A', 'B', 0,
    'C', 'C', 'C', 0,
    'C', 'C', 'D', 0
  };
  StringLcpContainer<unsigned char> ordered_cont(std::move(ordered_strings));
  StringLcpContainer<unsigned char> unordered_cont_1(std::move(unordered_strings_1));
  StringLcpContainer<unsigned char> unordered_cont_2(std::move(unordered_strings_2));
  
  REQUIRE(ordered_cont.make_string_set<UCharStringSet>().check_order() == 1);
  REQUIRE(unordered_cont_1.make_string_set<UCharStringSet>().check_order() == 0);
  REQUIRE(unordered_cont_2.make_string_set<UCharStringSet>().check_order() == 0);
}

