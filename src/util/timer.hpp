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
#include "mpi/allgather.hpp"

namespace dss_schimek {
  class EmptyTimer { // do not measure, used in alltoall-variants 
    public:
      static std::string getName() {
        return "EmptyTimer";
      }
      void start(const std::string& description,
          dsss::mpi::environment env = dsss::mpi::environment()) {}
      void end(const std::string& description,
          dsss::mpi::environment env = dsss::mpi::environment()) {}

  };
  class Timer {
    using Clock = std::chrono::high_resolution_clock;
    using PointInTime = std::chrono::time_point<Clock>;
    using TimeUnit = std::chrono::nanoseconds;
    using TimeIntervalDataType = size_t;

    public:
    static std::string getName() {
      return "Timer";
    }
    Timer(const std::string& prefix) : prefix(prefix) {};
    
    void add(const std::string& description, size_t value) {
      if (descriptionToValue.find(description) != descriptionToValue.end())
      {
        std::cout << description << " has already been added to timer" << std::endl;
        std::abort();
      }
      descriptionToValue.emplace(description, value);
    }

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

      const PointInTime endPointAfterBarrier = Clock::now();
      TimeIntervalDataType elapsedTotalTime =
        std::chrono::duration_cast<TimeUnit>(endPointAfterBarrier - startPoint).count();
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

    void writeToStream(std::stringstream& buffer, 
        const std::string& description,
        const std::string& type, 
        const TimeIntervalDataType time) {

      buffer <<  prefix
        << " "/*std::setw(alignmentLong) */<< ("operation=" + description)
        << " "/*std::setw(alignmentSmall)*/<< ("type=" + type)
        << " "/*std::setw(alignmentLong) */<< ("value=" + std::to_string(time))
        << std::endl;
    }
    
    void collectAndWriteToStream(std::stringstream& buffer, 
        const std::string& description,
        dsss::mpi::environment env = dsss::mpi::environment()) {

      size_t localValue = (*descriptionToValue.find(description)).second;
      std::vector<size_t> values = dsss::mpi::allgather(localValue, env);
      for (const size_t value : values) {
        buffer << prefix
          << " "/*std::setw(alignmentLong) */<< ("operation=" + description)
          << " "/*std::setw(alignmentSmall)*/<< ("type=number")
          << " "/*std::setw(alignmentLong) */<< ("value=" + std::to_string(value))
          << std::endl;
      }
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

      // time related values
      std::vector<std::string> descriptions;
      for (auto [description, not_used] : descriptionToStart)
        descriptions.push_back(description);
      for (const std::string& description : descriptions)
        writeToStream(buffer, description);
     
      // write remaining values 
      for (auto [description, value] : descriptionToValue)
        collectAndWriteToStream(buffer, description);
    }

    private:
    std::string prefix;
    std::map<std::string, size_t> descriptionToValue;
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
