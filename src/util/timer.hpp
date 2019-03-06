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
//j  class EmptyTimer { // do not measure, used in alltoall-variants 
//j    public:
//j      static std::string getName() {
//j        return "EmptyTimer";
//j      }
//j      void start(const std::string& description,
//j          dsss::mpi::environment env = dsss::mpi::environment()) {}
//j      void end(const std::string& description,
//j          dsss::mpi::environment env = dsss::mpi::environment()) {}
//j      void add(const std::string& , size_t) {}
//j
//j  };
//j
//j  
//j  class Timer_ {
//j    using Clock = std::chrono::high_resolution_clock;
//j    using PointInTime = std::chrono::time_point<Clock>;
//j    using TimeUnit = std::chrono::nanoseconds;
//j    using TimeIntervalDataType = size_t;
//j    struct DescriptionIteration {
//j      std::string description;
//j      size_t iteration;
//j      DescriptionIteration(const std::string& description, size_t iteration) 
//j        : description(description), iteration(iteration) {}
//j
//j      bool operator< (const DescriptionIteration& rhs) const {
//j        if (description == rhs.description)
//j          return iteration < rhs.iteration;
//j        return description < rhs.description;
//j      }
//j
//j      friend std::ostream& operator<<(std::ostream& stream, const DescriptionIteration& elem) {
//j        return stream << "[" << elem.description << ", " << elem.iteration << "]"; 
//j      }
//j    };
//j
//j    private: 
//j    
//j    
//j
//j    // make methods more generic to work with iteration and without iterations
//j    
//j    //template <typename KeyToValue, typename Key>
//j    //  void add(KeyToValue& keyToValue,const Key& key, size_t value) {
//j    //    if (keyToValue.find(key) != keyToValue.end())
//j    //      std::abort();
//j    //    keyToValue.emplace(key, value);
//j    //  }
//j
//j    //DescriptionIteration addString(DescriptionIteration descriptionIteration, const std::string& str) {
//j    //  descriptionIteration.first += str;
//j    //  return descriptionIteration;
//j    //}
//j    //
//j    //std::string addString(std::string description, const std::string& str) {
//j    //  description += str;
//j    //  return description;
//j    //}
//j
//j//template <typename KeyStartingPointMap, typename KeyDurationMap, typename Key>
//j//      void localStart(KeyStartingPointMap& keyToStartingPoint, KeyStartingPointMap& barrierKeyToStartingPoint, KeyDurationMap& keyToActiveTime, const Key& key,
//j//          dsss::mpi::environment env = dsss::mpi::environment())
//j//      {
//j//        if (!measurementEnabled)
//j//          return;
//j//        if (keyToStartingPoint.find(key) != keyToStartingPoint.end())
//j//          std::abort();
//j//        const PointInTime startBarrier = Clock::now();
//j//        //env.barrier();
//j//        const PointInTime start = Clock::now();
//j//        Key barrierKey = addString(key, "_Barrier");
//j//        barrierKeyToStartingPoint.emplace(barrierKey, start);
//j//        size_t elapsedActiveTimeInBarrier = 
//j//          std::chrono::duration_cast<std::chrono::nanoseconds>(start - startBarrier).count();
//j//        keyToActiveTime.emplace(barrierKey, elapsedActiveTimeInBarrier);
//j//
//j//        keyToStartingPoint.emplace(key, start);
//j//      }
//j
//j    Timer() = default;
//j
//j    Timer(const std::string& prefix) : prefix(prefix) {};
//j    void enableMeasurement() {
//j      measurementEnabled = true;
//j    }
//j    
//j    public:
//j    Timer(const Timer& timer) = delete;
//j    static std::string getName() {
//j      return "Timer";
//j    }
//j
//j    static Timer& timer() {
//j      static Timer timer;
//j      return timer;
//j    }
//j
//j    void setPrefix(const std::string& prefix_) {
//j      prefix = prefix_;
//j    }
//j
//j    void setInternIteration(const size_t internIteration) {
//j      curInternIteration = internIteration;
//j    }
//j
//j    void reset() {
//j     curInternIteration = 0;
//j     measurementEnabled = true;
//j     prefix = "";
//j     descriptionIterationToValue.clear();
//j     descriptionIterationToStart.clear();
//j     descriptionIterationToBarrierTime.clear();
//j     descriptionIterationToActiveTime.clear();
//j     descriptionIterationToTotalTime.clear();
//j    }
//j    
//j    void disableMeasurement() {
//j      measurementEnabled = false;
//j    }
//j
//j    void add(const std::string& description, size_t iteration, size_t value) {
//j      DescriptionIteration key(description, iteration);
//j      if (descriptionIterationToValue.find(key) != descriptionIterationToValue.end()) {
//j        std::cout << key << std::endl;
//j        std::cout << "already added" << std::endl;
//j        std::abort();
//j
//j      }
//j      descriptionIterationToValue.emplace(key, value);
//j    }
//j
//j    void add(const std::string& description, size_t value) {
//j      add(description, curInternIteration, value); 
//j    }
//j    
//j    void start(const std::string& description, size_t iteration) {
//j      if (!measurementEnabled)
//j        return;
//j      DescriptionIteration key(description, iteration);
//j      if (descriptionIterationToStart.find(key) != descriptionIterationToStart.end()) {
//j        std::cout << key << std::endl;
//j        std::cout << "already added" << std::endl;
//j        std::abort();
//j      }
//j      const PointInTime start = Clock::now();
//j      descriptionIterationToStart.emplace(key, start);
//j    }
//j
//j    void start(const std::string& description) {
//j      start(description, curInternIteration);
//j    }
//j
//j    void localStart(const std::string& description, size_t iteration) {
//j      start(description, iteration);
//j    }
//j
//j    void localStart(const std::string& description) {
//j      start(description, curInternIteration);
//j    }
//j    
//j    void end(const std::string& description, size_t iteration) {
//j        if (!measurementEnabled)
//j          return;
//j        DescriptionIteration key(description, iteration);
//j        auto itToPairInMap = descriptionIterationToStart.find(key);
//j        if (itToPairInMap == descriptionIterationToStart.end()) {
//j          std::cout << key << std::endl;
//j          std::cout << "no corresponding start" << std::endl;
//j          std::abort();
//j        }
//j
//j        PointInTime startPoint;
//j        std::tie(std::ignore, startPoint) = *itToPairInMap;
//j
//j        const PointInTime endPoint = Clock::now();
//j        size_t elapsedActiveTime =
//j          std::chrono::duration_cast<std::chrono::nanoseconds>(endPoint - startPoint).count();
//j        descriptionIterationToActiveTime.emplace(key, elapsedActiveTime);
//j        env.barrier();
//j        const PointInTime endPointAfterBarrier = Clock::now();
//j        TimeIntervalDataType elapsedTotalTime =
//j          std::chrono::duration_cast<TimeUnit>(endPointAfterBarrier - startPoint).count();
//j        descriptionIterationToTotalTime.emplace(key, elapsedTotalTime);
//j      }
//j
//j    void end(const std::string& description) {
//j      end(description, curInternIteration); 
//j    }
//j
//j    TimeIntervalDataType getLoss(const std::string& description, size_t iteration) {
//j        DescriptionIteration key(description, iteration);
//j        auto itKeyActiveTime = descriptionIterationToActiveTime.find(key);
//j        auto itKeyTotalTime = descriptionIterationToTotalTime.find(key);
//j        if (itKeyActiveTime == descriptionIterationToActiveTime.end() || 
//j            itKeyTotalTime == descriptionIterationToTotalTime.end()) {
//j          std::cout << key << std::endl;
//j          std::cout << "not added yet" << std::endl;
//j          std::abort();
//j        }
//j        return  (*itKeyTotalTime).second - (*itKeyActiveTime).second;
//j    }
//j   
//j    TimeIntervalDataType avgLoss(const std::string& description, size_t iteration, 
//j        dsss::mpi::environment env = dsss::mpi::environment()) {
//j      TimeIntervalDataType localLoss = getLoss(description, iteration);
//j      TimeIntervalDataType sum = dsss::mpi::allreduce_sum(localLoss);
//j      return sum / env.size();
//j    }
//j
//j    TimeIntervalDataType maxLoss(const std::string& description, size_t iteration) {
//j      TimeIntervalDataType localLoss = getLoss(description, iteration);
//j      return dsss::mpi::allreduce_max(localLoss);
//j    }
//j
//j    TimeIntervalDataType minLoss(const std::string& description, size_t iteration) {
//j      TimeIntervalDataType localLoss = getLoss(description, iteration);
//j      return dsss::mpi::allreduce_min(localLoss);
//j    }
//j
//j    TimeIntervalDataType avgTime(const std::string& description, size_t iteration) {
//j      DescriptionIteration key(description, iteration); 
//j      auto itKeyTime = descriptionIterationToActiveTime.find(key);
//j      if (itKeyTime == descriptionIterationToActiveTime.end()) {
//j          std::cout << key << std::endl;
//j          std::cout << "not added yet" << std::endl;
//j          std::abort();
//j      }
//j      TimeIntervalDataType sum = dsss::mpi::allreduce_sum((*itKeyTime).second);
//j      return sum / env.size();
//j    }
//j
//j    TimeIntervalDataType maxTime(const std::string& description, size_t iteration) {
//j      DescriptionIteration key(description, iteration); 
//j      auto itKeyTime = descriptionIterationToActiveTime.find(key);
//j      if (itKeyTime == descriptionIterationToActiveTime.end()) {
//j          std::cout << key << std::endl;
//j          std::cout << "not added yet" << std::endl;
//j          std::abort();
//j      }
//j      return dsss::mpi::allreduce_max((*itKeyTime).second);
//j    }
//j
//j    TimeIntervalDataType minTime(const std::string& description, size_t iteration) {
//j      DescriptionIteration key(description, iteration); 
//j      auto itKeyTime = descriptionIterationToActiveTime.find(key);
//j      if (itKeyTime == descriptionIterationToActiveTime.end()) {
//j          std::cout << key << std::endl;
//j          std::cout << "not added yet" << std::endl;
//j          std::abort();
//j      }
//j      return dsss::mpi::allreduce_min((*itKeyTime).second);
//j    }
//j
//j    //void writeToStream(std::stringstream& buffer, 
//j    //    const std::string& description,
//j    //    const std::string& type, 
//j    //    const TimeIntervalDataType time) {
//j
//j    //  buffer <<  prefix
//j    //    << " "/*std::setw(alignmentLong) */<< ("operation=" + description)
//j    //    << " "                             << ("internIteration=0")
//j    //    << " "/*std::setw(alignmentSmall)*/<< ("type=" + type)
//j    //    << " "/*std::setw(alignmentLong) */<< ("value=" + std::to_string(time))
//j    //    << std::endl;
//j    //}
//j    
//j    void writeToStream(std::stringstream& buffer, 
//j        const DescriptionIteration& descriptionIteration,
//j        const std::string& type, 
//j        const TimeIntervalDataType time) {
//j
//j      buffer <<  prefix
//j        << " " << ("operation=" + descriptionIteration.description)
//j        << " " << ("internIteration=" + std::to_string(descriptionIteration.iteration))
//j        << " " << ("type=" + type)
//j        << " " << ("value=" + std::to_string(time))
//j        << std::endl;
//j    }
//j
//j    //void collectAndWriteToStream(std::stringstream& buffer, 
//j    //    const std::string& description,
//j    //    dsss::mpi::environment env = dsss::mpi::environment()) {
//j
//j    //  size_t localValue = (*descriptionToValue.find(description)).second;
//j    //  std::vector<size_t> values = dsss::mpi::allgather(localValue, env);
//j    //  for (const size_t value : values) {
//j    //    buffer << prefix
//j    //      << " "/*std::setw(alignmentLong) */<< ("operation=" + description)
//j    //      << " "                             << ("internIteration=0")
//j    //      << " "/*std::setw(alignmentSmall)*/<< ("type=number")
//j    //      << " "/*std::setw(alignmentLong) */<< ("value=" + std::to_string(value))
//j    //      << std::endl;
//j    //  }
//j    //}
//j
//j    void collectAndWriteToStream(std::stringstream& buffer, 
//j        const DescriptionIteration key,
//j        dsss::mpi::environment env = dsss::mpi::environment()) {
//j
//j      size_t localValue = (*descriptionIterationToValue.find(key)).second;
//j      std::vector<size_t> values = dsss::mpi::allgather(localValue, env);
//j      for (const size_t value : values) {
//j        buffer << prefix
//j          << " " << ("operation=" + key.description)
//j          << " " << ("internIteration=" + std::to_string(key.iteration))
//j          << " " << ("type=number")
//j          << " " << ("value=" + std::to_string(value))
//j          << std::endl;
//j      }
//j    }
//j
//j    //void writeToStream(std::stringstream& buffer, const std::string& key, 
//j    //    dsss::mpi::environment env = dsss::mpi::environment()) {
//j
//j    //  writeToStream(buffer, key, "avgTime", avgTime(key));
//j    //  writeToStream(buffer, key, "maxTime", maxTime(key));
//j    //  writeToStream(buffer, key, "minTime", minTime(key));
//j    //  writeToStream(buffer, key, "avgLoss", avgLoss(key));
//j    //  writeToStream(buffer, key, "maxLoss", maxLoss(key));
//j    //  writeToStream(buffer, key, "minLoss", minLoss(key));
//j    //}
//j    
//j    void writeToStream(std::stringstream& buffer, const DescriptionIteration& key, 
//j        dsss::mpi::environment env = dsss::mpi::environment()) {
//j
//j      writeToStream(buffer, key, "avgTime", avgTime(key.description, key.iteration));
//j      writeToStream(buffer, key, "maxTime", maxTime(key.description, key.iteration));
//j      writeToStream(buffer, key, "minTime", minTime(key.description, key.iteration));
//j      writeToStream(buffer, key, "avgLoss", avgLoss(key.description, key.iteration));
//j      writeToStream(buffer, key, "maxLoss", maxLoss(key.description, key.iteration));
//j      writeToStream(buffer, key, "minLoss", minLoss(key.description, key.iteration));
//j    }
//j    
//j    //void writeToStreamNoLoss(std::stringstream& buffer, const std::string& key, 
//j    //    dsss::mpi::environment env = dsss::mpi::environment()) {
//j
//j    //  writeToStream(buffer, key, "avgTime", avgTime(key));
//j    //  writeToStream(buffer, key, "maxTime", maxTime(key));
//j    //  writeToStream(buffer, key, "minTime", minTime(key));
//j    //}
//j    
//j    //void writeToStreamNoLoss(std::stringstream& buffer, const std::pair<std::string, size_t>& key, 
//j    //    dsss::mpi::environment env = dsss::mpi::environment()) {
//j
//j    //  writeToStream(buffer, key, "avgTime", avgTime(key.first, key.second));
//j    //  writeToStream(buffer, key, "maxTime", maxTime(key.first, key.second));
//j    //  writeToStream(buffer, key, "minTime", minTime(key.first, key.second));
//j    //}
//j
//j    void writeToStream(std::stringstream& buffer,
//j        dsss::mpi::environment env = dsss::mpi::environment()) {
//j
//j      // time related values
//j      //std::vector<std::string> descriptions;
//j      std::vector<DescriptionIteration> descriptionIterations;
//j      //for (auto [description, not_used] : descriptionToStart)
//j      //  writeToStream(buffer, description); 
//j      for (auto [descriptionIteration, not_used] : descriptionIterationToStart)
//j        writeToStream(buffer, descriptionIteration);
//j      //for (auto [description, not_used] : barrierDescriptionToStart)
//j      //  writeToStreamNoLoss(buffer, description); 
//j      //for (auto [descriptionIteration, not_used] : barrierDescriptionIterationToStart)
//j      //  writeToStreamNoLoss(buffer, descriptionIteration);
//j
//j      //for (const std::string& description : descriptions)
//j      //  writeToStream(buffer, description);
//j      //for (const auto& descriptionIteration : descriptionIterations)
//j
//j
//j
//j      // write remaining values 
//j      //for (auto [description, value] : descriptionToValue)
//j      //  collectAndWriteToStream(buffer, description);
//j      for (auto [descriptionIteration, value] : descriptionIterationToValue)
//j        collectAndWriteToStream(buffer, descriptionIteration);
//j    }
//j
//j    private:
//j    bool measurementEnabled = true;
//j    size_t curInternIteration = 0;
//j    dsss::mpi::environment env;
//j
//j    std::string prefix;
//j    std::map<DescriptionIteration, size_t> descriptionIterationToValue;
//j    std::map<DescriptionIteration, PointInTime> descriptionIterationToStart;
//j    std::map<DescriptionIteration, TimeIntervalDataType> descriptionIterationToBarrierTime;
//j    std::map<DescriptionIteration, TimeIntervalDataType> descriptionIterationToActiveTime;
//j    std::map<DescriptionIteration, TimeIntervalDataType> descriptionIterationToTotalTime;
//j  };
  
  template<typename Key, typename Value>
  class Timer {
    using Clock = std::chrono::high_resolution_clock;
    using PointInTime = std::chrono::time_point<Clock>;
    using TimeUnit = std::chrono::nanoseconds;
    using TimeIntervalDataType = size_t;
    using PseudoKey = typename Key::PseudoKey;

    public: 
    
    Timer() = default;

        void start(const Key& key, const Value& value) {
      if (!measurementEnabled)
        return;
      if (keyToValue.find(key) != keyToValue.end()) {
        std::cout << "Key: " << key << " already added" << std::endl;
        std::abort();
      }
      const auto& itToValue = keyToValue.emplace(key, value);
      itToValue.first->second.setPseudoKeyCounter(incrementCounterPerPseudoKey(key));
      const auto&[itToStart, blabla] = keyToStart.emplace(key, PointInTime());

      // Start measurement
      const PointInTime start = Clock::now();
      itToStart->second = start;
    }

    void stop(const Key& key) {
      if (!measurementEnabled)
        return;

      const PointInTime endPoint = Clock::now();
      env.barrier();
      const PointInTime endPointAfterBarrier = Clock::now();
      //measurement stopped

      // check whether key is present
      auto it = keyToStart.find(key);
      if (it == keyToStart.end()){
        std::cout << "Key: " << key << " has no corresponding start" << std::endl;
        std::abort();
      }
      const PointInTime startPoint = it->second;

      TimeIntervalDataType elapsedActiveTime =
        std::chrono::duration_cast<std::chrono::nanoseconds>(endPoint - startPoint).count();
      TimeIntervalDataType elapsedTotalTime =
        std::chrono::duration_cast<TimeUnit>(endPointAfterBarrier - startPoint).count();

      keyToActiveTime.emplace(key, elapsedActiveTime);
      keyToTotalTime.emplace(key, elapsedTotalTime);
    }
    
    template <typename OutputIterator>
      void collect(OutputIterator out) {
        for (auto& [key, value] : keyToValue) {
          value.setType("avgTime");  
          value.setValue(avgTime(key));
          out = std::make_pair(key, value);

          value.setType("maxTime");  
          value.setValue(maxTime(key));
          out = std::make_pair(key, value);

          value.setType("minTime");  
          value.setValue(minTime(key));
          out = std::make_pair(key, value);
 
          value.setType("avgLoss");  
          value.setValue(avgLoss(key));
          out = std::make_pair(key, value);
 
          value.setType("maxLoss");  
          value.setValue(maxLoss(key));
          out = std::make_pair(key, value);
 
          value.setType("minLoss");  
          value.setValue(minLoss(key));
          out = std::make_pair(key, value);
        }
      }

    private:
    bool measurementEnabled = true;
    dsss::mpi::environment env;

    size_t incrementCounterPerPseudoKey(const Key& key) {
      auto itCounter = pseudoKeyToCounter.find(key.pseudoKey());
      if (itCounter == pseudoKeyToCounter.end()) {
        pseudoKeyToCounter.emplace(key.pseudoKey(), 1u);
        return 0;
      } else {
        return itCounter->second++;
      }
    }

    TimeIntervalDataType getLoss(const Key& key) {
      auto itKeyActiveTime = keyToActiveTime.find(key);
      auto itKeyTotalTime = keyToTotalTime.find(key);
      if (itKeyActiveTime == keyToActiveTime.end() || itKeyTotalTime == keyToTotalTime.end()) {
        std::cout << key << std::endl;
        std::cout << "Key " << key << " not present" << std::endl;
        std::abort();
      }
      return  (*itKeyTotalTime).second - (*itKeyActiveTime).second;
    }

    TimeIntervalDataType avgLoss(const Key& key,
        dsss::mpi::environment env = dsss::mpi::environment()) {

      TimeIntervalDataType localLoss = getLoss(key);
      TimeIntervalDataType sum = dsss::mpi::allreduce_sum(localLoss);
      return sum / env.size();
    }

    TimeIntervalDataType maxLoss(const Key& key) {
      TimeIntervalDataType localLoss = getLoss(key);
      return dsss::mpi::allreduce_max(localLoss);
    }

    TimeIntervalDataType minLoss(const Key& key) {
      TimeIntervalDataType localLoss = getLoss(key);
      return dsss::mpi::allreduce_min(localLoss);
    }

    TimeIntervalDataType avgTime(const Key& key) {
      auto itKeyTime = keyToActiveTime.find(key);
      if (itKeyTime == keyToActiveTime.end()) {
        std::cout << "Key " << key << "not present" << std::endl;
        std::abort();
      }
      TimeIntervalDataType sum = dsss::mpi::allreduce_sum((*itKeyTime).second);
      return sum / env.size();
    }

    TimeIntervalDataType maxTime(const Key& key) {
      auto itKeyTime = keyToActiveTime.find(key);
      if (itKeyTime == keyToActiveTime.end()) {
        std::cout << key << std::endl;
        std::cout << "Key " << key << " not present" << std::endl;
        std::abort();
      }
      return dsss::mpi::allreduce_max((*itKeyTime).second);
    }

    TimeIntervalDataType minTime(const Key& key) {
      auto itKeyTime = keyToActiveTime.find(key);
      if (itKeyTime == keyToActiveTime.end()) {
        std::cout << key << std::endl;
        std::cout << "Key " << key << " not present" << std::endl;
        std::abort();
      }
      return dsss::mpi::allreduce_min((*itKeyTime).second);
    }

    std::map<Key, Value> keyToValue;
    std::map<PseudoKey, size_t> pseudoKeyToCounter;
    std::map<Key, PointInTime> keyToStart;
    std::map<Key, TimeIntervalDataType> keyToActiveTime;
    std::map<Key, TimeIntervalDataType> keyToTotalTime;
  };

}
