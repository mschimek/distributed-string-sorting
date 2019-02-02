#pragma once 

#include "strings/stringptr.hpp"
#include "strings/stringset.hpp"
#include "mpi/shift.hpp" 
#include "mpi/allreduce.hpp"
#include "mpi/environment.hpp"

namespace dss_schimek {

  template<typename StringPtr>
  bool is_sorted(const StringPtr& strptr, dsss::mpi::environment env = dsss::mpi::environment())
  {
    using StringSet = typename StringPtr::StringSet;
    using String = typename StringSet::String;
    const StringSet& ss = strptr.active();
    bool is_locally_sorted = ss.check_order();

    if (env.size() == 1)
      return is_locally_sorted;
    
    size_t has_strings = ss.size() > 0;
    size_t number_PE_with_data = dsss::mpi::allreduce_sum(has_strings);

    if (number_PE_with_data <= 1)
      return is_locally_sorted;
   
    int own_min_number = has_strings ? env.rank() : env.size();
    int own_max_number = has_strings ? env.rank() : -1;
    int min_PE_with_data = dsss::mpi::allreduce_min(own_min_number);
    int max_PE_with_data = dsss::mpi::allreduce_max(own_max_number);
    const bool is_left_shift = true;
    std::vector<unsigned char> greater_string;
    std::vector<unsigned char> smaller_string; 
    unsigned char* front;
    unsigned char* back;
    if constexpr (std::is_same<StringSet, UCharLengthStringSet>::value) {
      front = has_strings ? (*ss.begin()).string : ss.empty_string().string;
      back = has_strings ? (*(ss.end() - 1)).string : ss.empty_string().string;
      greater_string = dss_schimek::mpi::shift_string<is_left_shift>(
          front, env
          ); 
      smaller_string = dss_schimek::mpi::shift_string<!is_left_shift>(
          back, env
          ); 
    } else {
      front = has_strings ? *ss.begin() : ss.empty_string();
      back = has_strings ? *(ss.end() - 1) : ss.empty_string();
      greater_string = dss_schimek::mpi::shift_string<is_left_shift>(
          front, env
          ); 
      smaller_string = dss_schimek::mpi::shift_string<!is_left_shift>(
          back, env
          );
    }


    bool is_overall_sorted = is_locally_sorted;
    if (!has_strings)
      return dsss::mpi::allreduce_and(is_overall_sorted, env);

    if (env.rank() != min_PE_with_data)
      is_overall_sorted &= dss_schimek::scmp(smaller_string.data(), front) <= 0;
    if (env.rank() != max_PE_with_data)
      is_overall_sorted &= dss_schimek::scmp(back, greater_string.data()) <= 0;

    return dsss::mpi::allreduce_and(is_overall_sorted, env);
  }

  template<typename StringPtr>
  bool is_complete_and_sorted(const StringPtr& strptr,
      size_t initial_local_num_chars,
      size_t current_local_num_chars,
      size_t initial_local_num_strings,
      size_t current_local_num_strings,
      dsss::mpi::environment env = dsss::mpi::environment())
  {
      if (env.size() == 0)
        return is_sorted(strptr);

      const size_t initial_total_num_chars = dsss::mpi::allreduce_sum(initial_local_num_chars);
      const size_t initial_total_num_strings = dsss::mpi::allreduce_sum(initial_local_num_strings);

      const size_t current_total_num_chars = dsss::mpi::allreduce_sum(current_local_num_chars);
      const size_t current_total_num_strings = dsss::mpi::allreduce_sum(current_local_num_strings);

      if (initial_total_num_chars != current_total_num_chars) {
        std::cout << "initial total num chars: " << initial_total_num_chars <<
           " current_total_num_chars: " << current_total_num_chars << std::endl;
        std::cout << "We've lost some chars" << std::endl;
        return false;
      }
      if (initial_total_num_strings != current_total_num_strings) {
        std::cout << "initial total num strings: " << initial_total_num_strings << " current total num strings: " << current_total_num_strings << std::endl;
        std::cout << "We've lost some strings" << std::endl;
        return false;
      }
      return is_sorted(strptr);
  }
}
