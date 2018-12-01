#pragma once

#include <algorithm>
#include <chrono>
#include <vector>
#include <map>
#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <numeric>
#include <tuple>

#include "mpi/allreduce.hpp"

namespace dss_schimek {
  class Timer {
    using Clock = std::chrono::high_resolution_clock;
    using PointInTime = std::chrono::time_point<Clock>;
    using TimeUnit = std::chrono::nanoseconds;
    using TimeIntervalDataType = size_t;

    public:
      Timer(const std::string& prefix) : prefix(prefix) {};
      //void start(const std::string& description) {
      //  if (descriptionToStart.find(description) != descriptionToStart.end())
      //    std::abort();
      //  const PointInTime start = Clock::now();
      //  descriptionToStart.emplace(description, start);
      //}
      void start(const std::string& description,
                          dsss::mpi::environment env = dsss::mpi::environment()) {
        if (descriptionToStart.find(description) != descriptionToStart.end())
          std::abort();
        env.barrier();
        const PointInTime start = Clock::now();
        descriptionToStart.emplace(description, start);
      }

      void end(const std::string& description,
                        dsss::mpi::environment env = dsss::mpi::environment()) {
        auto itToPairInMap = descriptionToStart.find(description);
        if (itToPairInMap == descriptionToStart.end())
          std::abort();

        PointInTime startPoint;
        std::tie(std::ignore, startPoint) = *itToPairInMap;

        const PointInTime endPoint = Clock::now();
        size_t elapsedActiveTime =
          std::chrono::duration_cast<std::chrono::nanoseconds>(endPoint - startPoint).count();
        descriptionToActiveTime.emplace(description, elapsedActiveTime);

        env.barrier();

        TimeIntervalDataType elapsedTotalTime =
          std::chrono::duration_cast<TimeUnit>(endPoint - startPoint).count();
        descriptionToTotalTime.emplace(description, elapsedTotalTime);
      }

      TimeIntervalDataType getLoss(const std::string& description) {
        auto itDescriptionActiveTime = descriptionToActiveTime.find(description);
        auto itDescriptionTotalTime = descriptionToTotalTime.find(description);
        if (itDescriptionActiveTime == descriptionToActiveTime.end() || 
            itDescriptionTotalTime == descriptionToTotalTime.end())
          std::abort();

        return  (*itDescriptionTotalTime).second - (*itDescriptionActiveTime).second;
      }

      TimeIntervalDataType avgLoss(const std::string& description,
          dsss::mpi::environment env = dsss::mpi::environment()) {

        TimeIntervalDataType localLoss = getLoss(description);
        TimeIntervalDataType sum = dsss::mpi::allreduce_sum(localLoss);
        return sum / env.size();
      }

      TimeIntervalDataType maxLoss(const std::string& description,
          dsss::mpi::environment env = dsss::mpi::environment()) {

        TimeIntervalDataType localLoss = getLoss(description);
        return dsss::mpi::allreduce_max(localLoss);
      }

      TimeIntervalDataType minLoss(const std::string& description,
          dsss::mpi::environment env = dsss::mpi::environment()) {

        TimeIntervalDataType localLoss = getLoss(description);
        return dsss::mpi::allreduce_min(localLoss);
      }

      TimeIntervalDataType avgTime(const std::string& description,
          dsss::mpi::environment env = dsss::mpi::environment()) {
        auto itDescriptionTime = descriptionToActiveTime.find(description);
        if (itDescriptionTime == descriptionToActiveTime.end())
          std::abort();
        TimeIntervalDataType sum = dsss::mpi::allreduce_sum((*itDescriptionTime).second);
        return sum / env.size();
      }
      
      TimeIntervalDataType maxTime(const std::string& description,
          dsss::mpi::environment env = dsss::mpi::environment()) {
        auto itDescriptionTime = descriptionToActiveTime.find(description);
        if (itDescriptionTime == descriptionToActiveTime.end())
          std::abort();
        return dsss::mpi::allreduce_max((*itDescriptionTime).second);
      }
      
      TimeIntervalDataType minTime(const std::string& description,
          dsss::mpi::environment env = dsss::mpi::environment()) {
        auto itDescriptionTime = descriptionToActiveTime.find(description);
        if (itDescriptionTime == descriptionToActiveTime.end())
          std::abort();
        return dsss::mpi::allreduce_min((*itDescriptionTime).second);
      }

      //void end(const std::string& description) {
      //  auto itToPairInMap = descriptionToStart.find(description);
      //  if (itToPairInMap == descriptionToStart.end())
      //    std::abort();

      //  PointInTime startPoint;
      //  std::tie(std::ignore, startPoint) = *itToPairInMap;

      //  const PointInTime endPoint = Clock::now();
      //  TimeIntervalDataType elapsed_time =
      //    std::chrono::duration_cast<TimeUnit>(endPoint - startPoint).count();
      //  descriptionToTime.emplace(description, elapsed_time);
      //} 

      void writeToStream(std::stringstream& buffer, 
          const std::string& description,
          const std::string& type, 
          const TimeIntervalDataType time) {

        buffer <<  prefix
          << std::setw(alignmentLong)  << ("operation=" + description)
          << std::setw(alignmentSmall) << ("type=" + type)
          << std::setw(alignmentSmall) << ("time=" + std::to_string(time))
          << std::endl;;
      }

      void writeToStream(std::stringstream& buffer, const std::string& description, 
          dsss::mpi::environment env = dsss::mpi::environment()) {

        writeToStream(buffer, description, "avgTime", avgTime(description));
        writeToStream(buffer, description, "maxTime", maxTime(description));
        writeToStream(buffer, description, "minTime", minTime(description));
        writeToStream(buffer, description, "avgLoss", avgLoss(description));
        writeToStream(buffer, description, "maxLoss", maxLoss(description));
        writeToStream(buffer, description, "minLoss", minLoss(description));
      }

      void writeToStream(std::stringstream& buffer,
          dsss::mpi::environment env = dsss::mpi::environment()) {

        std::vector<std::string> descriptions;
        for (auto [description, not_used] : descriptionToStart)
          descriptions.push_back(description);
        for (const std::string& description : descriptions)
          writeToStream(buffer, description);
      }

      //void print_sum(const std::vector<std::string>& descriptions) const {

      //  size_t overall_time = 0;
      //  std::vector<std::string> containedDescriptions;
      //  for (const std::string& description : descriptions) {
      //    auto candidate = descriptionToTime.find(description);
      //    if (candidate != descriptionToTime.end()) {
      //      containedDescriptions.emplace_back(description);
      //      overall_time += (*candidate).second;
      //    }
      //  }

      //  print_n(50, "-");
      //  size_t max_description_length = get_max_string_length(containedDescriptions);
      //  for (const std::string& description : containedDescriptions) {
      //      std::cout << prefix 
      //                << std::setw(max_description_length + 2)
      //                << description
      //                << std::endl;
      //  }
      //  print_n(50, "-");
      //  std::cout << prefix
      //            << std::setw(max_description_length + 2)
      //            << "result: "
      //            << std::setw(12)
      //            << overall_time 
      //            << std::endl;
      //}
      //
      //void print_n(const size_t n, const std::string& string) const{
      //  for (size_t i = 0; i < n; ++i)
      //    std::cout << string;
      //  std::cout << std::endl;
      //}

      //void print() const {
      //  std::vector<std::string> descriptions;
      //  for (auto [description, not_used] : descriptionToTime)
      //    descriptions.push_back(description);
      //  size_t max_description_length = get_max_string_length(descriptions);

      //  
      //  print_n(50, "#");
      //  for (auto [description, interval_length] : descriptionToTime) {
      //    std::cout << "RESULT "
      //              << prefix 
      //              << std::setw(max_description_length + 2) 
      //              << "name=" << description  
      //              << std::setw(12) 
      //              << "time="
      //              << interval_length
      //              << std::endl;
      //  }
      //  print_n(50, "#");
      //}

      //void print_synchronized(dsss::mpi::environment env = dsss::mpi::environment()) {
      //   std::vector<std::string> descriptions;
      //  for (auto [description, not_used] : descriptionToSynchronizedTime)
      //    descriptions.push_back(description);
      //  size_t max_description_length = get_max_string_length(descriptions);

      //  
      //  print_n(50, "-");
      //  std::cout << "synchronized" << std::endl;
      //  print_n(50, "-");
      //  for (auto [description, interval_length] : descriptionToSynchronizedTime) {
      //    std::cout << prefix
      //              << std::setw(max_description_length + 2) 
      //              << description 
      //              << std::setw(12) 
      //              << interval_length
      //              << std::endl;
      //  }
      //  print_n(50, "-");

      //}
      //void synchronize(dsss::mpi::environment env = dsss::mpi::environment()) {
      //  descriptionToSynchronizedTime.clear();
      //  for (auto [description, interval_length] : descriptionToTime) {
      //    size_t global_max_interval_length = dsss::mpi::allreduce_max(interval_length);
      //    descriptionToSynchronizedTime.emplace(description, global_max_interval_length);
      //  }
      //}


    private:
      const size_t alignmentLong = 40;
      const size_t alignmentSmall = 15;
      std::string prefix;
      std::map<std::string, PointInTime> descriptionToStart;
      std::map<std::string, TimeIntervalDataType> descriptionToTime;
      std::map<std::string, TimeIntervalDataType> descriptionToActiveTime;
      std::map<std::string, TimeIntervalDataType> descriptionToTotalTime;
      std::map<std::string, TimeIntervalDataType> descriptionToSynchronizedTime;

      size_t get_max_string_length(const std::vector<std::string>& strings) const {
        size_t maxStringLength = 0;
        for (const std::string& str : strings)
          maxStringLength = std::max(str.size(), maxStringLength);
        return maxStringLength;
      }
  };
}
