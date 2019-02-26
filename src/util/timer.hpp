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
    using DescriptionIteration = std::pair<std::string, size_t>;

    private: 

    // make methods more generic to work with iteration and without iterations
    
    template <typename KeyToValue, typename Key>
      void add(KeyToValue& keyToValue,const Key& key, size_t value) {
        if (keyToValue.find(key) != keyToValue.end())
          std::abort();
        keyToValue.emplace(key, value);
      }

    DescriptionIteration addString(DescriptionIteration descriptionIteration, const std::string& str) {
      descriptionIteration.first += str;
      return descriptionIteration;
    }
    
    std::string addString(std::string description, const std::string& str) {
      description += str;
      return description;
    }

    template <typename KeyStartingPointMap, typename KeyDurationMap, typename Key>
      void start(KeyStartingPointMap& keyToStartingPoint, KeyStartingPointMap& barrierKeyToStartingPoint, KeyDurationMap& keyToActiveTime, const Key& key,
          dsss::mpi::environment env = dsss::mpi::environment())
      {
        if (!measurementEnabled)
          return;
        if (keyToStartingPoint.find(key) != keyToStartingPoint.end())
          std::abort();
        const PointInTime startBarrier = Clock::now();
        env.barrier();
        const PointInTime start = Clock::now();
        Key barrierKey = addString(key, "_Barrier");
        barrierKeyToStartingPoint.emplace(barrierKey, start);
        size_t elapsedActiveTimeInBarrier = 
          std::chrono::duration_cast<std::chrono::nanoseconds>(start - startBarrier).count();
        keyToActiveTime.emplace(barrierKey, elapsedActiveTimeInBarrier);

        keyToStartingPoint.emplace(key, start);
      }

    template<typename KeyStartingPointMap, typename KeyDurationMap, typename Key>
      void end(KeyStartingPointMap& keyToStartingPoint, KeyDurationMap& keyToActiveTime,
          KeyDurationMap& keyToTotalTime, const Key& key,
          dsss::mpi::environment env = dsss::mpi::environment()) {
        if (!measurementEnabled)
          return;
        auto itToPairInMap = keyToStartingPoint.find(key);
        if (itToPairInMap == keyToStartingPoint.end())
          std::abort();

        PointInTime startPoint;
        std::tie(std::ignore, startPoint) = *itToPairInMap;

        const PointInTime endPoint = Clock::now();
        size_t elapsedActiveTime =
          std::chrono::duration_cast<std::chrono::nanoseconds>(endPoint - startPoint).count();
        keyToActiveTime.emplace(key, elapsedActiveTime);

        //env.barrier();

        const PointInTime endPointAfterBarrier = Clock::now();
        TimeIntervalDataType elapsedTotalTime =
          std::chrono::duration_cast<TimeUnit>(endPointAfterBarrier - startPoint).count();
        keyToTotalTime.emplace(key, elapsedTotalTime);
      }

    template <typename KeyToDurationMap, typename Key>
      TimeIntervalDataType getLoss(const KeyToDurationMap& keyToActiveTime, 
          const KeyToDurationMap& keyToTotalTime, 
          const Key& key) {

        auto itKeyActiveTime = keyToActiveTime.find(key);
        auto itKeyTotalTime = keyToTotalTime.find(key);
        if (itKeyActiveTime == keyToActiveTime.end() || 
            itKeyTotalTime == keyToTotalTime.end())
          std::abort();
        return  (*itKeyTotalTime).second - (*itKeyActiveTime).second;
      }

    template <typename KeyToDurationTime, typename Key>
      TimeIntervalDataType avgTime(const KeyToDurationTime& keyToActiveTime, const Key& key,
          dsss::mpi::environment env = dsss::mpi::environment()) {
        auto itKeyTime = keyToActiveTime.find(key);
        if (itKeyTime == keyToActiveTime.end())
          std::abort();
        TimeIntervalDataType sum = dsss::mpi::allreduce_sum((*itKeyTime).second);
        return sum / env.size();
      }

    template <typename KeyToDurationTime, typename Key>
      TimeIntervalDataType maxTime(KeyToDurationTime& keyToActiveTime, const Key& key,
          dsss::mpi::environment env = dsss::mpi::environment()) {
        auto itKeyTime = keyToActiveTime.find(key);
        if (itKeyTime == keyToActiveTime.end())
          std::abort();
        return dsss::mpi::allreduce_max((*itKeyTime).second);
      }

    template <typename KeyToDurationTime, typename Key>
      TimeIntervalDataType minTime(KeyToDurationTime& keyToActiveTime, const Key& key,
          dsss::mpi::environment env = dsss::mpi::environment()) { auto itKeyTime = keyToActiveTime.find(key);
        if (itKeyTime == keyToActiveTime.end())
          std::abort();
        return dsss::mpi::allreduce_min((*itKeyTime).second);
      }

    public:
    static std::string getName() {
      return "Timer";
    }
    Timer(const std::string& prefix) : prefix(prefix) {};
    void enableMeasurement() {
      measurementEnabled = true;
    }
    
    void disableMeasurement() {
      measurementEnabled = false;
    }

    void add(const std::string& description, size_t value) {
      add(descriptionToValue, description, value);
    }

    void add(const std::string& description, size_t iteration, size_t value) {
      add(descriptionIterationToValue, make_pair(description, iteration), value);
    }

    void start(const std::string& description, size_t iteration) {
      start(descriptionIterationToStart, barrierDescriptionIterationToStart, descriptionIterationToActiveTime, make_pair(description, iteration)); 
    }

    void start(const std::string& description) {
      start(descriptionToStart, barrierDescriptionToStart, descriptionToActiveTime, description); 
    }

    void end(const std::string& description) {
      end(descriptionToStart, descriptionToActiveTime, descriptionToTotalTime, description);
    }

    void end(const std::string& description, size_t iteration) {
      end(descriptionIterationToStart, descriptionIterationToActiveTime, descriptionIterationToTotalTime, make_pair(description, iteration));
    }

    TimeIntervalDataType getLoss(const std::string& description) {
      return getLoss(descriptionToActiveTime, descriptionToTotalTime, description);
    }
    
    TimeIntervalDataType getLoss(const std::string& description, size_t iteration) {
      return getLoss(descriptionIterationToActiveTime, descriptionIterationToTotalTime, make_pair(description, iteration));
    }

    TimeIntervalDataType avgLoss(const std::string& description,
        dsss::mpi::environment env = dsss::mpi::environment()) {
      TimeIntervalDataType localLoss = getLoss(description);
      TimeIntervalDataType sum = dsss::mpi::allreduce_sum(localLoss);
      return sum / env.size();
    }
    
    TimeIntervalDataType avgLoss(const std::string& description, size_t iteration, 
        dsss::mpi::environment env = dsss::mpi::environment()) {
      TimeIntervalDataType localLoss = getLoss(description, iteration);
      TimeIntervalDataType sum = dsss::mpi::allreduce_sum(localLoss);
      return sum / env.size();
    }

    TimeIntervalDataType maxLoss(const std::string& description) {
      TimeIntervalDataType localLoss = getLoss(description);
      return dsss::mpi::allreduce_max(localLoss);
    }
    
    TimeIntervalDataType maxLoss(const std::string& description, size_t iteration) {
      TimeIntervalDataType localLoss = getLoss(description, iteration);
      return dsss::mpi::allreduce_max(localLoss);
    }

    TimeIntervalDataType minLoss(const std::string& description) {
      TimeIntervalDataType localLoss = getLoss(description);
      return dsss::mpi::allreduce_min(localLoss);
    }
    
    TimeIntervalDataType minLoss(const std::string& description, size_t iteration) {
      TimeIntervalDataType localLoss = getLoss(description, iteration);
      return dsss::mpi::allreduce_min(localLoss);
    }

    TimeIntervalDataType avgTime(const std::string& description) {
      return avgTime(descriptionToActiveTime, description);
    }
    
    TimeIntervalDataType avgTime(const std::string& description, size_t iteration) {
      return avgTime(descriptionIterationToActiveTime, make_pair(description, iteration));
    }

    TimeIntervalDataType maxTime(const std::string& description) {
      return maxTime(descriptionToActiveTime, description);   
    }
    
    TimeIntervalDataType maxTime(const std::string& description, size_t iteration) {
      return maxTime(descriptionIterationToActiveTime, make_pair(description, iteration));   
    }

    TimeIntervalDataType minTime(const std::string& description) {
      return minTime(descriptionToActiveTime, description);
    }
    
    TimeIntervalDataType minTime(const std::string& description, size_t iteration) {
      return minTime(descriptionIterationToActiveTime, make_pair(description, iteration));
    }

    void writeToStream(std::stringstream& buffer, 
        const std::string& description,
        const std::string& type, 
        const TimeIntervalDataType time) {

      buffer <<  prefix
        << " "/*std::setw(alignmentLong) */<< ("operation=" + description)
        << " "                             << ("internIteration=0")
        << " "/*std::setw(alignmentSmall)*/<< ("type=" + type)
        << " "/*std::setw(alignmentLong) */<< ("value=" + std::to_string(time))
        << std::endl;
    }
    
    void writeToStream(std::stringstream& buffer, 
        const std::pair<std::string, size_t>& descriptionIteration,
        const std::string& type, 
        const TimeIntervalDataType time) {

      buffer <<  prefix
        << " "/*std::setw(alignmentLong) */<< ("operation=" + descriptionIteration.first)
        << " "                             << ("internIteration=" + std::to_string(descriptionIteration.second))
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
          << " "                             << ("internIteration=0")
          << " "/*std::setw(alignmentSmall)*/<< ("type=number")
          << " "/*std::setw(alignmentLong) */<< ("value=" + std::to_string(value))
          << std::endl;
      }
    }

    void collectAndWriteToStream(std::stringstream& buffer, 
        const std::pair<std::string, size_t>& descriptionIteration,
        dsss::mpi::environment env = dsss::mpi::environment()) {

      size_t localValue = (*descriptionIterationToValue.find(descriptionIteration)).second;
      std::vector<size_t> values = dsss::mpi::allgather(localValue, env);
      for (const size_t value : values) {
        buffer << prefix
          << " "/*std::setw(alignmentLong) */<< ("operation=" + descriptionIteration.first)
          << " "                             << ("internIteration=" + std::to_string(descriptionIteration.second))
          << " "/*std::setw(alignmentSmall)*/<< ("type=number")
          << " "/*std::setw(alignmentLong) */<< ("value=" + std::to_string(value))
          << std::endl;
      }
    }

    void writeToStream(std::stringstream& buffer, const std::string& key, 
        dsss::mpi::environment env = dsss::mpi::environment()) {

      writeToStream(buffer, key, "avgTime", avgTime(key));
      writeToStream(buffer, key, "maxTime", maxTime(key));
      writeToStream(buffer, key, "minTime", minTime(key));
      writeToStream(buffer, key, "avgLoss", avgLoss(key));
      writeToStream(buffer, key, "maxLoss", maxLoss(key));
      writeToStream(buffer, key, "minLoss", minLoss(key));
    }
    
    void writeToStream(std::stringstream& buffer, const std::pair<std::string, size_t>& key, 
        dsss::mpi::environment env = dsss::mpi::environment()) {

      writeToStream(buffer, key, "avgTime", avgTime(key.first, key.second));
      writeToStream(buffer, key, "maxTime", maxTime(key.first, key.second));
      writeToStream(buffer, key, "minTime", minTime(key.first, key.second));
      writeToStream(buffer, key, "avgLoss", avgLoss(key.first, key.second));
      writeToStream(buffer, key, "maxLoss", maxLoss(key.first, key.second));
      writeToStream(buffer, key, "minLoss", minLoss(key.first, key.second));
    }
    
    void writeToStreamNoLoss(std::stringstream& buffer, const std::string& key, 
        dsss::mpi::environment env = dsss::mpi::environment()) {

      writeToStream(buffer, key, "avgTime", avgTime(key));
      //writeToStream(buffer, key, "maxTime", maxTime(key));
      //writeToStream(buffer, key, "minTime", minTime(key));
    }
    
    void writeToStreamNoLoss(std::stringstream& buffer, const std::pair<std::string, size_t>& key, 
        dsss::mpi::environment env = dsss::mpi::environment()) {

      writeToStream(buffer, key, "avgTime", avgTime(key.first, key.second));
      //writeToStream(buffer, key, "maxTime", maxTime(key.first, key.second));
      //writeToStream(buffer, key, "minTime", minTime(key.first, key.second));
    }

    void writeToStream(std::stringstream& buffer,
        dsss::mpi::environment env = dsss::mpi::environment()) {

      // time related values
      std::vector<std::string> descriptions;
      std::vector<DescriptionIteration> descriptionIterations;
      for (auto [description, not_used] : descriptionToStart)
        writeToStream(buffer, description); 
      for (auto [descriptionIteration, not_used] : descriptionIterationToStart)
        writeToStream(buffer, descriptionIteration);
      for (auto [description, not_used] : barrierDescriptionToStart)
        writeToStreamNoLoss(buffer, description); 
      for (auto [descriptionIteration, not_used] : barrierDescriptionIterationToStart)
        writeToStreamNoLoss(buffer, descriptionIteration);

      //for (const std::string& description : descriptions)
      //  writeToStream(buffer, description);
      //for (const auto& descriptionIteration : descriptionIterations)



      // write remaining values 
      for (auto [description, value] : descriptionToValue)
        collectAndWriteToStream(buffer, description);
      for (auto [descriptionIteration, value] : descriptionIterationToValue)
        collectAndWriteToStream(buffer, descriptionIteration);

    }

    private:
    bool measurementEnabled = true;

    std::string prefix;
    std::map<std::string, size_t> descriptionToValue;
    std::map<DescriptionIteration, size_t> descriptionIterationToValue;

    std::map<std::string, PointInTime> descriptionToStart;
    std::map<DescriptionIteration, PointInTime> descriptionIterationToStart;

    std::map<std::string, PointInTime> barrierDescriptionToStart;
    std::map<DescriptionIteration, PointInTime> barrierDescriptionIterationToStart;

    std::map<std::string, TimeIntervalDataType> descriptionToTime;
    std::map<DescriptionIteration, TimeIntervalDataType> descriptionIterationToTime;

    std::map<std::string, TimeIntervalDataType> descriptionToBarrierTime;
    std::map<DescriptionIteration, TimeIntervalDataType> descriptionIterationToBarrierTime;

    std::map<std::string, TimeIntervalDataType> descriptionToActiveTime;
    std::map<DescriptionIteration, TimeIntervalDataType> descriptionIterationToActiveTime;

    std::map<std::string, TimeIntervalDataType> descriptionToTotalTime;
    std::map<DescriptionIteration, TimeIntervalDataType> descriptionIterationToTotalTime;

    std::map<std::string, TimeIntervalDataType> descriptionToSynchronizedTime;
    std::map<DescriptionIteration, TimeIntervalDataType> descriptionIterationToSynchronizedTime;

    size_t get_max_string_length(const std::vector<std::string>& strings) const {
      size_t maxStringLength = 0;
      for (const std::string& str : strings)
        maxStringLength = std::max(str.size(), maxStringLength);
      return maxStringLength;
    }
  };
}
