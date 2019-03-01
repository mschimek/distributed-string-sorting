#pragma once

#include <vector>
#include <algorithm>
#include <numeric>

#include "mpi/environment.hpp"

namespace dss_schimek {

  template <typename StringSet>
    class  SampleSplittersNumStringsPolicy
    {
      public:
        static std::string getName() {
          return "NumStrings";
        }
      protected:
        std::vector<typename StringSet::Char> sample_splitters(const StringSet& ss,
            dsss::mpi::environment env = dsss::mpi::environment()) {

          using Char = typename StringSet::Char;
          using String = typename StringSet::String;

          const size_t local_num_strings = ss.size();
          const size_t nr_splitters = std::min<size_t>(env.size() - 1, local_num_strings);
          const size_t splitter_dist = local_num_strings / (nr_splitters + 1);
          std::vector<Char> raw_splitters;

          for (size_t i = 1; i <= nr_splitters; ++i) {
            const String splitter = ss[ss.begin() + i * splitter_dist];
            const size_t splitterLength = ss.get_length(splitter);
            const size_t usedSplitterLength = splitterLength > 100 ? 100 : splitterLength;
            std::copy_n(ss.get_chars(splitter, 0), usedSplitterLength,
                std::back_inserter(raw_splitters));
            raw_splitters.push_back(0);
          }
          return raw_splitters; 
        }
    };

  template <typename StringSet>
    class  SampleSplittersNumCharsPolicy
    {
      public: 
        static std::string getName() {
          return "NumChars";
        }
      protected:
        std::vector<typename StringSet::Char> sample_splitters(const StringSet& ss,
            dsss::mpi::environment env = dsss::mpi::environment()) {

          using Char = typename StringSet::Char;
          using String = typename StringSet::String;

          const size_t num_chars = std::accumulate(ss.begin(), ss.end(), 0, 
              [&ss](const size_t& sum, const String& str) {
              return sum + ss.get_length(str);
              });

          const size_t local_num_strings = ss.size();
          const size_t nr_splitters = std::min<size_t>(env.size() - 1, local_num_strings);
          const size_t splitter_dist = num_chars / (nr_splitters + 1);
          std::vector<Char> raw_splitters;

          size_t string_index = 0;
          for (size_t i = 1; i <= nr_splitters; ++i) {
            size_t num_chars_seen = 0;
            while (num_chars_seen < splitter_dist) {
              num_chars_seen += ss.get_length(ss[ss.begin() + string_index]);
              ++string_index;
            }

            const String splitter = ss[ss.begin() + string_index - 1];
            std::copy_n(ss.get_chars(splitter, 0), ss.get_length(splitter) + 1,
                std::back_inserter(raw_splitters));
          }
          return raw_splitters; 
        }
    };
}
