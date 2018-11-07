#include "../external/dsss/string_sorting/distributed/merge_sort.hpp" 
#include "../external/dsss/util/string_set.hpp"
#include "util/random_string_generator.hpp"

#include "sequential/bingmann-mkqs.hpp"
#include "sequential/bingmann-radix_sort.hpp"

#include "string_sorting/distributed/merge_sort.hpp"

#include "utils/draw_from_distribution.hpp"

void sample_sort_kurpicz(dsss::string_set& local_string_set) 
{            
  dsss::sample_sort::sample_sort<bingmann::bingmann_msd_CE0>(local_string_set);
}

dsss::random_string_set create_random_string_set(const size_t size, 
  const size_t min_length, const size_t max_length)
{
  return dsss::random_string_set(size, min_length, max_length);
}


int main()
{
  dsss::mpi::environment env;
  constexpr size_t size = 200;
  constexpr size_t min_string_length = 10;
  constexpr size_t max_string_length = 20;
 
  dsss::string_set strings_to_be_sorted = create_random_string_set(size, 
                                                                   min_string_length, 
                                                                   max_string_length);
  std::vector<double> weights = dss::get_weights(std::vector<int>{1,6});
  std::vector<int> chosenElements = dss::draw_elements_from_distribution(50, weights);
  for(const auto& element : chosenElements)
    std::cout << element << " ";
  std::cout << std::endl;
  env.barrier();
  sample_sort_kurpicz(strings_to_be_sorted);
  env.finalize();
}
