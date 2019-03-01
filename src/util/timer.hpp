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
      void add(const std::string& , size_t) {}

  };

  
  class Timer {
    using Clock = std::chrono::high_resolution_clock;
    using PointInTime = std::chrono::time_point<Clock>;
    using TimeUnit = std::chrono::nanoseconds;
    using TimeIntervalDataType = size_t;
    struct DescriptionIteration {
      std::string description;
      size_t iteration;
      DescriptionIteration(const std::string& description, size_t iteration) 
        : description(description), iteration(iteration) {}

      bool operator< (const DescriptionIteration& rhs) const {
        if (description == rhs.description)
          return iteration < rhs.iteration;
        return description < rhs.description;
      }

      friend std::ostream& operator<<(std::ostream& stream, const DescriptionIteration& elem) {
        return stream << "[" << elem.description << ", " << elem.iteration << "]"; 
      }
    };

    private: 
    
    

    // make methods more generic to work with iteration and without iterations
    
    //template <typename KeyToValue, typename Key>
    //  void add(KeyToValue& keyToValue,const Key& key, size_t value) {
    //    if (keyToValue.find(key) != keyToValue.end())
    //      std::abort();
    //    keyToValue.emplace(key, value);
    //  }

    //DescriptionIteration addString(DescriptionIteration descriptionIteration, const std::string& str) {
    //  descriptionIteration.first += str;
    //  return descriptionIteration;
    //}
    //
    //std::string addString(std::string description, const std::string& str) {
    //  description += str;
    //  return description;
    //}

//template <typename KeyStartingPointMap, typename KeyDurationMap, typename Key>
//      void localStart(KeyStartingPointMap& keyToStartingPoint, KeyStartingPointMap& barrierKeyToStartingPoint, KeyDurationMap& keyToActiveTime, const Key& key,
//          dsss::mpi::environment env = dsss::mpi::environment())
//      {
//        if (!measurementEnabled)
//          return;
//        if (keyToStartingPoint.find(key) != keyToStartingPoint.end())
//          std::abort();
//        const PointInTime startBarrier = Clock::now();
//        //env.barrier();
//        const PointInTime start = Clock::now();
//        Key barrierKey = addString(key, "_Barrier");
//        barrierKeyToStartingPoint.emplace(barrierKey, start);
//        size_t elapsedActiveTimeInBarrier = 
//          std::chrono::duration_cast<std::chrono::nanoseconds>(start - startBarrier).count();
//        keyToActiveTime.emplace(barrierKey, elapsedActiveTimeInBarrier);
//
//        keyToStartingPoint.emplace(key, start);
//      }

    Timer() = default;

    Timer(const std::string& prefix) : prefix(prefix) {};
    void enableMeasurement() {
      measurementEnabled = true;
    }
    
    public:
    Timer(const Timer& timer) = delete;
    static std::string getName() {
      return "Timer";
    }

    static Timer& timer() {
      static Timer timer;
      return timer;
    }

    void setPrefix(const std::string& prefix_) {
      prefix = prefix_;
    }

    void setInternIteration(const size_t internIteration) {
      curInternIteration = internIteration;
    }

    void reset() {
     curInternIteration = 0;
     measurementEnabled = true;
     prefix = "";
     descriptionIterationToValue.clear();
     descriptionIterationToStart.clear();
     descriptionIterationToBarrierTime.clear();
     descriptionIterationToActiveTime.clear();
     descriptionIterationToTotalTime.clear();
    }
    
    void disableMeasurement() {
      measurementEnabled = false;
    }

    void add(const std::string& description, size_t iteration, size_t value) {
      DescriptionIteration key(description, iteration);
      if (descriptionIterationToValue.find(key) != descriptionIterationToValue.end()) {
        std::cout << key << std::endl;
        std::cout << "already added" << std::endl;
        std::abort();

      }
      descriptionIterationToValue.emplace(key, value);
    }

    void add(const std::string& description, size_t value) {
      add(description, curInternIteration, value); 
    }
    
    void start(const std::string& description, size_t iteration) {
      if (!measurementEnabled)
        return;
      DescriptionIteration key(description, iteration);
      if (descriptionIterationToStart.find(key) != descriptionIterationToStart.end()) {
        std::cout << key << std::endl;
        std::cout << "already added" << std::endl;
        std::abort();
      }
      const PointInTime start = Clock::now();
      descriptionIterationToStart.emplace(key, start);
    }

    void start(const std::string& description) {
      start(description, curInternIteration);
    }

    void localStart(const std::string& description, size_t iteration) {
      start(description, iteration);
    }

    void localStart(const std::string& description) {
      start(description, curInternIteration);
    }
    
    void end(const std::string& description, size_t iteration) {
        if (!measurementEnabled)
          return;
        DescriptionIteration key(description, iteration);
        auto itToPairInMap = descriptionIterationToStart.find(key);
        if (itToPairInMap == descriptionIterationToStart.end()) {
          std::cout << key << std::endl;
          std::cout << "no corresponding start" << std::endl;
          std::abort();
        }

        PointInTime startPoint;
        std::tie(std::ignore, startPoint) = *itToPairInMap;

        const PointInTime endPoint = Clock::now();
        size_t elapsedActiveTime =
          std::chrono::duration_cast<std::chrono::nanoseconds>(endPoint - startPoint).count();
        descriptionIterationToActiveTime.emplace(key, elapsedActiveTime);
        env.barrier();
        const PointInTime endPointAfterBarrier = Clock::now();
        TimeIntervalDataType elapsedTotalTime =
          std::chrono::duration_cast<TimeUnit>(endPointAfterBarrier - startPoint).count();
        descriptionIterationToTotalTime.emplace(key, elapsedTotalTime);
      }

    void end(const std::string& description) {
      end(description, curInternIteration); 
    }

    TimeIntervalDataType getLoss(const std::string& description, size_t iteration) {
        DescriptionIteration key(description, iteration);
        auto itKeyActiveTime = descriptionIterationToActiveTime.find(key);
        auto itKeyTotalTime = descriptionIterationToTotalTime.find(key);
        if (itKeyActiveTime == descriptionIterationToActiveTime.end() || 
            itKeyTotalTime == descriptionIterationToTotalTime.end()) {
          std::cout << key << std::endl;
          std::cout << "not added yet" << std::endl;
          std::abort();
        }
        return  (*itKeyTotalTime).second - (*itKeyActiveTime).second;
    }
   
    TimeIntervalDataType avgLoss(const std::string& description, size_t iteration, 
        dsss::mpi::environment env = dsss::mpi::environment()) {
      TimeIntervalDataType localLoss = getLoss(description, iteration);
      TimeIntervalDataType sum = dsss::mpi::allreduce_sum(localLoss);
      return sum / env.size();
    }

    TimeIntervalDataType maxLoss(const std::string& description, size_t iteration) {
      TimeIntervalDataType localLoss = getLoss(description, iteration);
      return dsss::mpi::allreduce_max(localLoss);
    }

    TimeIntervalDataType minLoss(const std::string& description, size_t iteration) {
      TimeIntervalDataType localLoss = getLoss(description, iteration);
      return dsss::mpi::allreduce_min(localLoss);
    }

    TimeIntervalDataType avgTime(const std::string& description, size_t iteration) {
      DescriptionIteration key(description, iteration); 
      auto itKeyTime = descriptionIterationToActiveTime.find(key);
      if (itKeyTime == descriptionIterationToActiveTime.end()) {
          std::cout << key << std::endl;
          std::cout << "not added yet" << std::endl;
          std::abort();
      }
      TimeIntervalDataType sum = dsss::mpi::allreduce_sum((*itKeyTime).second);
      return sum / env.size();
    }

    TimeIntervalDataType maxTime(const std::string& description, size_t iteration) {
      DescriptionIteration key(description, iteration); 
      auto itKeyTime = descriptionIterationToActiveTime.find(key);
      if (itKeyTime == descriptionIterationToActiveTime.end()) {
          std::cout << key << std::endl;
          std::cout << "not added yet" << std::endl;
          std::abort();
      }
      return dsss::mpi::allreduce_max((*itKeyTime).second);
    }

    TimeIntervalDataType minTime(const std::string& description, size_t iteration) {
      DescriptionIteration key(description, iteration); 
      auto itKeyTime = descriptionIterationToActiveTime.find(key);
      if (itKeyTime == descriptionIterationToActiveTime.end()) {
          std::cout << key << std::endl;
          std::cout << "not added yet" << std::endl;
          std::abort();
      }
      return dsss::mpi::allreduce_min((*itKeyTime).second);
    }

    //void writeToStream(std::stringstream& buffer, 
    //    const std::string& description,
    //    const std::string& type, 
    //    const TimeIntervalDataType time) {

    //  buffer <<  prefix
    //    << " "/*std::setw(alignmentLong) */<< ("operation=" + description)
    //    << " "                             << ("internIteration=0")
    //    << " "/*std::setw(alignmentSmall)*/<< ("type=" + type)
    //    << " "/*std::setw(alignmentLong) */<< ("value=" + std::to_string(time))
    //    << std::endl;
    //}
    
    void writeToStream(std::stringstream& buffer, 
        const DescriptionIteration& descriptionIteration,
        const std::string& type, 
        const TimeIntervalDataType time) {

      buffer <<  prefix
        << " " << ("operation=" + descriptionIteration.description)
        << " " << ("internIteration=" + std::to_string(descriptionIteration.iteration))
        << " " << ("type=" + type)
        << " " << ("value=" + std::to_string(time))
        << std::endl;
    }

    //void collectAndWriteToStream(std::stringstream& buffer, 
    //    const std::string& description,
    //    dsss::mpi::environment env = dsss::mpi::environment()) {

    //  size_t localValue = (*descriptionToValue.find(description)).second;
    //  std::vector<size_t> values = dsss::mpi::allgather(localValue, env);
    //  for (const size_t value : values) {
    //    buffer << prefix
    //      << " "/*std::setw(alignmentLong) */<< ("operation=" + description)
    //      << " "                             << ("internIteration=0")
    //      << " "/*std::setw(alignmentSmall)*/<< ("type=number")
    //      << " "/*std::setw(alignmentLong) */<< ("value=" + std::to_string(value))
    //      << std::endl;
    //  }
    //}

    void collectAndWriteToStream(std::stringstream& buffer, 
        const DescriptionIteration key,
        dsss::mpi::environment env = dsss::mpi::environment()) {

      size_t localValue = (*descriptionIterationToValue.find(key)).second;
      std::vector<size_t> values = dsss::mpi::allgather(localValue, env);
      for (const size_t value : values) {
        buffer << prefix
          << " " << ("operation=" + key.description)
          << " " << ("internIteration=" + std::to_string(key.iteration))
          << " " << ("type=number")
          << " " << ("value=" + std::to_string(value))
          << std::endl;
      }
    }

    //void writeToStream(std::stringstream& buffer, const std::string& key, 
    //    dsss::mpi::environment env = dsss::mpi::environment()) {

    //  writeToStream(buffer, key, "avgTime", avgTime(key));
    //  writeToStream(buffer, key, "maxTime", maxTime(key));
    //  writeToStream(buffer, key, "minTime", minTime(key));
    //  writeToStream(buffer, key, "avgLoss", avgLoss(key));
    //  writeToStream(buffer, key, "maxLoss", maxLoss(key));
    //  writeToStream(buffer, key, "minLoss", minLoss(key));
    //}
    
    void writeToStream(std::stringstream& buffer, const DescriptionIteration& key, 
        dsss::mpi::environment env = dsss::mpi::environment()) {

      writeToStream(buffer, key, "avgTime", avgTime(key.description, key.iteration));
      writeToStream(buffer, key, "maxTime", maxTime(key.description, key.iteration));
      writeToStream(buffer, key, "minTime", minTime(key.description, key.iteration));
      writeToStream(buffer, key, "avgLoss", avgLoss(key.description, key.iteration));
      writeToStream(buffer, key, "maxLoss", maxLoss(key.description, key.iteration));
      writeToStream(buffer, key, "minLoss", minLoss(key.description, key.iteration));
    }
    
    //void writeToStreamNoLoss(std::stringstream& buffer, const std::string& key, 
    //    dsss::mpi::environment env = dsss::mpi::environment()) {

    //  writeToStream(buffer, key, "avgTime", avgTime(key));
    //  writeToStream(buffer, key, "maxTime", maxTime(key));
    //  writeToStream(buffer, key, "minTime", minTime(key));
    //}
    
    //void writeToStreamNoLoss(std::stringstream& buffer, const std::pair<std::string, size_t>& key, 
    //    dsss::mpi::environment env = dsss::mpi::environment()) {

    //  writeToStream(buffer, key, "avgTime", avgTime(key.first, key.second));
    //  writeToStream(buffer, key, "maxTime", maxTime(key.first, key.second));
    //  writeToStream(buffer, key, "minTime", minTime(key.first, key.second));
    //}

    void writeToStream(std::stringstream& buffer,
        dsss::mpi::environment env = dsss::mpi::environment()) {

      // time related values
      //std::vector<std::string> descriptions;
      std::vector<DescriptionIteration> descriptionIterations;
      //for (auto [description, not_used] : descriptionToStart)
      //  writeToStream(buffer, description); 
      for (auto [descriptionIteration, not_used] : descriptionIterationToStart)
        writeToStream(buffer, descriptionIteration);
      //for (auto [description, not_used] : barrierDescriptionToStart)
      //  writeToStreamNoLoss(buffer, description); 
      //for (auto [descriptionIteration, not_used] : barrierDescriptionIterationToStart)
      //  writeToStreamNoLoss(buffer, descriptionIteration);

      //for (const std::string& description : descriptions)
      //  writeToStream(buffer, description);
      //for (const auto& descriptionIteration : descriptionIterations)



      // write remaining values 
      //for (auto [description, value] : descriptionToValue)
      //  collectAndWriteToStream(buffer, description);
      for (auto [descriptionIteration, value] : descriptionIterationToValue)
        collectAndWriteToStream(buffer, descriptionIteration);
    }

    private:
    bool measurementEnabled = true;
    size_t curInternIteration = 0;
    dsss::mpi::environment env;

    std::string prefix;
    std::map<DescriptionIteration, size_t> descriptionIterationToValue;
    std::map<DescriptionIteration, PointInTime> descriptionIterationToStart;
    std::map<DescriptionIteration, TimeIntervalDataType> descriptionIterationToBarrierTime;
    std::map<DescriptionIteration, TimeIntervalDataType> descriptionIterationToActiveTime;
    std::map<DescriptionIteration, TimeIntervalDataType> descriptionIterationToTotalTime;
  };
}
