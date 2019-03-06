#pragma once

#include <ostream>
#include <vector>
#include "util/valueTracker.hpp"
#include "util/timer.hpp"

namespace dss_schimek::measurement {

struct PhaseCounterPerPhaseRoundDescriptionTypeValue {
  using PseudoKey = std::string;

  std::string phase;
  size_t counterPerPhase;
  size_t round;
  std::string description;
  std::string type;
  size_t value;

  PhaseCounterPerPhaseRoundDescriptionTypeValue(const std::string& phase,
      const size_t counterPerPhase,
      const size_t round,
      const std::string& description,
      const std::string& type,
      const size_t value) 
    : phase(phase), counterPerPhase(counterPerPhase), description(description), round(round), type(type), value(value) {}

  const std::string& pseudoKey() const {
    return phase;
  }

  void setPseudoKeyCounter(size_t counter) {
    counterPerPhase = counter; 
  }

  void setType(const std::string& type_) {
    type = type_;
  }

  void setValue(size_t value_) {
    value = value_;
  }

  size_t getValue() const {
    return value;
  }

  friend std::ostream& operator<<(std::ostream& stream, PhaseCounterPerPhaseRoundDescriptionTypeValue& elem) {
    return stream << "[" << elem.phase << ", " << elem.counterPerPhase << ", " << elem.round << ", " << elem.description << ", " << elem.value << "]";
  }
};



struct PhaseRoundDescription {
  using PseudoKey = std::string;

  std::string phase;
  size_t round;
  std::string description;

  PhaseRoundDescription(const std::string& phase, const size_t round, const std::string& description) 
    : phase(phase), round(round), description(description) {}

  const std::string& pseudoKey() const {
    return phase;
  }

  bool operator< (const PhaseRoundDescription& rhs) const {
    return std::tie(phase, round, description) < std::tie(rhs.phase, rhs.round, rhs.description); 
  }

  friend std::ostream& operator<<(std::ostream& stream, const PhaseRoundDescription& elem) {
    return stream << "[" << elem.phase << ", " << elem.round << ", " << elem.description << "]";
  }
};

struct CounterPerPhaseTypeValue {
  size_t counterPerPhase;
  std::string type;
  size_t value;
  CounterPerPhaseTypeValue(const size_t counterPerPhase, const std::string& type, const size_t value)
    : counterPerPhase(counterPerPhase), type(type), value(value) {}
  
  void setPseudoKeyCounter(size_t counter) {
    counterPerPhase = counter;
  }

  void setType(const std::string& type_) {
    type = type_;
  }

  void setValue(size_t value_) {
    value = value_;
  }
};

  class MeasuringTool {
    // Columns:
    // Prefix | Phase | CounterPerPhase | Round | Description | Type | Value
    //
    // NonTimer: Duplicates are allowed (but since counterPerPhase is incremented there are no duplicates if all fields are taken into account)
    // Timer: Key = (Phase, Round, Description) Value = (CounterPerPhase, Type, Value)

    using NonTimerRecord = PhaseCounterPerPhaseRoundDescriptionTypeValue;
    // NonTimerRecord must contain type PseudoKey and functions pseudoKey() and setPseudoKeyCounter, setValue(), getValue()
    using TimerKey = PhaseRoundDescription;
    // TimerKey must contain type PseudoKey and function pseudoKey() 
    using TimerValue = CounterPerPhaseTypeValue;
    // TimerValue must contain functions setType() and setPseudoKeyCounter 

    struct OutputFormat {
      std::string prefix;
      std::string phase;
      size_t counterPerPhase;
      size_t round;
      std::string description;
      std::string type;
      size_t value;

      friend std::ostream& operator<<(std::ostream& stream, const OutputFormat& outputFormat) {
        return stream << outputFormat.prefix 
          << " phase=" << outputFormat.phase 
          << " counterPerPhase=" << outputFormat.counterPerPhase
          << " round=" << outputFormat.round 
          << " operation=" << outputFormat.description
          << " type=" << outputFormat.type
          << " value=" << outputFormat.value;
      }
    };
    
    public:
    static MeasuringTool& measuringTool() {
      static MeasuringTool measuringTool;
      return measuringTool;
    }
    
    void reset() {
      timer = Timer<TimerKey, TimerValue>();
      nonTimer = NonTimer<NonTimerRecord>();
    }
    void add(size_t value) {
      nonTimer.add(NonTimerRecord(curPhase, 0u, curRound, "unkown", "type", value));
    }

    void add(size_t value, const std::string& description) {
      nonTimer.add(NonTimerRecord(curPhase, 0u, curRound, description, "type", value));
    }

    void start(const std::string& description) {
      timer.start(TimerKey(curPhase, curRound, description), TimerValue(0u, "", 0u));
    }

    void stop(const std::string& description) {
      timer.stop(TimerKey(curPhase, curRound, description));
    }

    std::vector<OutputFormat> collect() {
      const size_t numMeasurements = 1000;

      std::vector<OutputFormat> data;
      data.reserve(10 * numMeasurements);
      std::vector<NonTimerRecord> nonTimerRecords;
      std::vector<std::pair<TimerKey, TimerValue>> timerRecords;
      nonTimerRecords.reserve(numMeasurements);
      timerRecords.reserve(6 * numMeasurements);

      nonTimer.collect(std::back_inserter(nonTimerRecords));
      timer.collect(std::back_inserter(timerRecords));

      for (const auto& nonTimerRecord : nonTimerRecords)
        data.push_back({prefix, nonTimerRecord.phase, nonTimerRecord.counterPerPhase, nonTimerRecord.round, nonTimerRecord.description, nonTimerRecord.type, nonTimerRecord.value});

      for (const auto& [timerKey, timerValue] : timerRecords)
        data.push_back({prefix, timerKey.phase, timerValue.counterPerPhase, timerKey.round, timerKey.description, timerValue.type, timerValue.value});
      return data;
    }

    void writeToStream(std::ostream& stream) {
      for (const auto& data : collect()) 
        stream << data << std::endl;
    }
    void setPhase(const std::string& phase) {
      curPhase = phase;
    }
    void setRound(size_t round) {
      curRound = round;
    }

    private:
      std::string prefix;
      std::string curPhase;
      size_t curRound;
      NonTimer<NonTimerRecord> nonTimer;
      Timer<TimerKey, TimerValue> timer;
  };
}
