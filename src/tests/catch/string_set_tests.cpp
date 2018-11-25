#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
#include <catch/catch.hpp>

#include <vector>
#include "strings/stringset.hpp"
#include "strings/stringcontainer.hpp"

TEST_CASE( "check length property") {
  using namespace dss_schimek;
  using StringSet = UCharLengthStringSet;
  std::vector<typename StringSet::Char> raw_strings = {
      'A', 'A', 'B', 'C', 'A', 0,
      'B', 'A', 0,
      'B', 'A', 'B', 0,
      0,
      'C', 'C', 'C', 0,
      'C', 'C', 'D', 0
    };
  StringLcpContainer<StringSet> container(std::move(raw_strings));
  const StringSet ss = container.make_string_set();
  auto i = ss.begin();
  REQUIRE(ss[i++].length == 5);
  REQUIRE(ss[i++].length == 2);
  REQUIRE(ss[i++].length == 3);
  REQUIRE(ss[i++].length == 0);
  REQUIRE(ss[i++].length == 3);
  REQUIRE(ss[i++].length == 3);

}

TEST_CASE( "check method in string sets", "check method") {
  using namespace dss_schimek;
  SECTION("UChar") {
    using StringSet = UCharStringSet;
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
    StringLcpContainer<StringSet> ordered_cont(std::move(ordered_strings));
    StringLcpContainer<StringSet> unordered_cont_1(std::move(unordered_strings_1));
    StringLcpContainer<StringSet> unordered_cont_2(std::move(unordered_strings_2));

    REQUIRE(ordered_cont.make_string_set().check_order() == 1);
    REQUIRE(unordered_cont_1.make_string_set().check_order() == 0);
    REQUIRE(unordered_cont_2.make_string_set().check_order() == 0);

  }
  SECTION("UCharLength") {
    using StringSet = UCharLengthStringSet;
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
    StringLcpContainer<StringSet> ordered_cont(std::move(ordered_strings));
    StringLcpContainer<StringSet> unordered_cont_1(std::move(unordered_strings_1));
    StringLcpContainer<StringSet> unordered_cont_2(std::move(unordered_strings_2));

    REQUIRE(ordered_cont.make_string_set().check_order() == 1);
    REQUIRE(unordered_cont_1.make_string_set().check_order() == 0);
    REQUIRE(unordered_cont_2.make_string_set().check_order() == 0);
  }
}


