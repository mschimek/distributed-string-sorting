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
    const StringSet& ss = strptr.active();
    using String = UCharStringSet::Char*;
    std::cout << "rank: " << env.rank() << " init" << std::endl;
    bool is_locally_sorted = ss.check_order();

    std::cout << "rank: " << env.rank() << " is_locally_sorted: " << is_locally_sorted << std::endl;
    if (env.size() == 1)
      return is_locally_sorted;
    
    size_t has_strings = ss.size() > 0;
    size_t number_PE_with_data = dsss::mpi::allreduce_sum(has_strings);
    std::cout << "rank: " << env.rank() << " number_PE_with_data: " << number_PE_with_data << std::endl;

    if (number_PE_with_data <= 1)
      return is_locally_sorted;
   
    int own_min_number = has_strings ? env.rank() : env.size();
    int own_max_number = has_strings ? env.rank() : -1;
    int min_PE_with_data = dsss::mpi::allreduce_min(own_min_number);
    int max_PE_with_data = dsss::mpi::allreduce_max(own_max_number);
    std::cout << "rank: " << env.rank() << " min_PE_with_data: " << min_PE_with_data << std::endl;
    std::cout << "rank: " << env.rank() << " max_PE_with_data: " << max_PE_with_data << std::endl;
    const bool is_left_shift = true;
    const String front = has_strings ? *ss.begin() : ss.empty_string();
    const String back = has_strings ? *(ss.end() - 1) : ss.empty_string();
    std::vector<unsigned char> greater_string = dss_schimek::mpi::shift_string<is_left_shift>(
        front, env
        ); 
    std::vector<unsigned char> smaller_string = dss_schimek::mpi::shift_string<!is_left_shift>(
        back, env
        ); 


    std::cout << "rank: " << env.rank() << " greater_string: " << greater_string.data() << " back: " << back <<  
      " smaller_string: " << smaller_string.data()  << " front: " << front << std::endl;

    std::cout << "rank " << env.rank() << " smaller vs front " << dss_schimek::scmp(smaller_string.data(), front) << std::endl;
    bool is_overall_sorted = is_locally_sorted;
    if (!has_strings)
      return dsss::mpi::allreduce_and(is_overall_sorted, env);

    if (env.rank() != min_PE_with_data)
      is_overall_sorted &= dss_schimek::scmp(smaller_string.data(), front) <= 0;
    if (env.rank() != max_PE_with_data)
      is_overall_sorted &= dss_schimek::scmp(back, greater_string.data()) <= 0;

    return dsss::mpi::allreduce_and(is_overall_sorted, env);
  }
}
