#pragma once

#include <algorithm>
#include <chrono>
#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <iomanip>
#include <numeric>
#include <tuple>

#include "mpi/allreduce.hpp"

namespace dss_schimek {
  class Timer {
    using Clock = std::chrono::high_resolution_clock;
    using PointInTime = std::chrono::time_point<Clock>;

    public:
      Timer(const std::string& prefix) : prefix(prefix) {};
      void start(const std::string& description) {
        if (descriptionToStart.find(description) != descriptionToStart.end())
          std::abort();
        const PointInTime start = Clock::now();
        descriptionToStart.emplace(description, start);
      }

      void end(const std::string& description) {
        auto itToPairInMap = descriptionToStart.find(description);
        if (itToPairInMap == descriptionToStart.end())
          std::abort();

        PointInTime startPoint;
        std::tie(std::ignore, startPoint) = *itToPairInMap;

        const PointInTime endPoint = Clock::now();
        size_t elapsed_time =
          std::chrono::duration_cast<std::chrono::nanoseconds>(endPoint - startPoint).count();
        descriptionToTime.emplace(description, elapsed_time);
      } 

      void print_sum(const std::vector<std::string>& descriptions) const {

        size_t overall_time = 0;
        std::vector<std::string> containedDescriptions;
        for (const std::string& description : descriptions) {
          auto candidate = descriptionToTime.find(description);
          if (candidate != descriptionToTime.end()) {
            containedDescriptions.emplace_back(description);
            overall_time += (*candidate).second;
          }
        }

        print_n(50, "-");
        size_t max_description_length = get_max_string_length(containedDescriptions);
        for (const std::string& description : containedDescriptions) {
            std::cout << prefix 
                      << std::setw(max_description_length + 2)
                      << description
                      << std::endl;
        }
        print_n(50, "-");
        std::cout << prefix
                  << std::setw(max_description_length + 2)
                  << "result: "
                  << std::setw(12)
                  << overall_time 
                  << std::endl;
      }
      
      void print_n(const size_t n, const std::string& string) const{
        for (size_t i = 0; i < n; ++i)
          std::cout << string;
        std::cout << std::endl;
      }

      void print() const {
        std::vector<std::string> descriptions;
        for (auto [description, not_used] : descriptionToTime)
          descriptions.push_back(description);
        size_t max_description_length = get_max_string_length(descriptions);

        
        print_n(50, "-");
        for (auto [description, interval_length] : descriptionToTime) {
          std::cout << prefix
                    << std::setw(max_description_length + 2) 
                    << description 
                    << std::setw(12) 
                    << interval_length
                    << std::endl;
        }
        print_n(50, "-");
      }

      void print_synchronized(dsss::mpi::environment env = dsss::mpi::environment()) {
         std::vector<std::string> descriptions;
        for (auto [description, not_used] : descriptionToSynchronizedTime)
          descriptions.push_back(description);
        size_t max_description_length = get_max_string_length(descriptions);

        
        print_n(50, "-");
        std::cout << "synchronized" << std::endl;
        print_n(50, "-");
        for (auto [description, interval_length] : descriptionToSynchronizedTime) {
          std::cout << prefix
                    << std::setw(max_description_length + 2) 
                    << description 
                    << std::setw(12) 
                    << interval_length
                    << std::endl;
        }
        print_n(50, "-");

      }
      void synchronize(dsss::mpi::environment env = dsss::mpi::environment()) {
        descriptionToSynchronizedTime.clear();
        for (auto [description, interval_length] : descriptionToTime) {
          size_t global_max_interval_length = dsss::mpi::allreduce_max(interval_length);
          descriptionToSynchronizedTime.emplace(description, global_max_interval_length);
        }
      }


      //size_t synchronize(const std::string& description, 
      //    dsss::mpi::environment env = dsss::mpi::environment()) {

      //  auto itToPair = descriptionToTime.find(description);
      //  if (descriptionToTime.find(description) == descriptionToTime.end())  
      //    std::abort();

      //  size_t interval_length = itToPair->second;
      //  return dsss::mpi::allreduce_max(interval_length);

      //}
    private:
      std::string prefix;
      std::map<std::string, PointInTime> descriptionToStart;
      std::map<std::string, size_t> descriptionToTime;
      std::map<std::string, size_t> descriptionToSynchronizedTime;

      size_t get_max_string_length(const std::vector<std::string>& strings) const {
        size_t maxStringLength = 0;
        for (const std::string& str : strings)
          maxStringLength = std::max(str.size(), maxStringLength);
        return maxStringLength;
      }
  };
}
