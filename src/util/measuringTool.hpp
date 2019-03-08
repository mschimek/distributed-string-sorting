#pragma once

#include <ostream>
#include <vector>
#include "util/nonTimer.hpp"
#include "util/timer.hpp"

namespace dss_schimek {
  namespace measurement {

    struct PhaseCounterPerPhaseRoundDescriptionTypeRawCommunicationValue {
      using PseudoKey = std::string;

      std::string phase;
      size_t counterPerPhase;
      size_t round;
      std::string description;
      std::string type;
      bool rawCommunication;
      size_t value;

      PhaseCounterPerPhaseRoundDescriptionTypeRawCommunicationValue(const std::string& phase,
          const size_t counterPerPhase,
          const size_t round,
          const std::string& description,
          const std::string& type,
          const bool rawCommunication,
          const size_t value) 
        : phase(phase), counterPerPhase(counterPerPhase), description(description), round(round), 
          type(type), rawCommunication(rawCommunication), value(value) {}

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

      friend std::ostream& operator<<(std::ostream& stream, PhaseCounterPerPhaseRoundDescriptionTypeRawCommunicationValue& elem) {
        return stream << "[" << elem.phase << ", " << elem.counterPerPhase << ", " << elem.round << ", " << elem.description << ", " << elem.rawCommunication <<  ", " << elem.value << "]";
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

    struct CounterPerPhaseTypeRawCommunicationValue {
      size_t counterPerPhase;
      std::string type;
      bool rawCommunication;
      size_t value;
      CounterPerPhaseTypeRawCommunicationValue(const size_t counterPerPhase, const std::string& type, const bool rawCommunication, const size_t value)
        : counterPerPhase(counterPerPhase), type(type), rawCommunication(rawCommunication), value(value) {}

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
      // Prefix | Phase | CounterPerPhase | Round | Description | Type | RawCommunication | Value
      //
      // NonTimer: Duplicates are allowed (but since counterPerPhase is incremented there are no duplicates if all fields are taken into account)
      // Timer: Key = (Phase, Round, Description) Value = (CounterPerPhase, Type, RawCommunication, Value)

      using NonTimerRecord = PhaseCounterPerPhaseRoundDescriptionTypeRawCommunicationValue;
      // NonTimerRecord must contain type PseudoKey and functions pseudoKey() and setPseudoKeyCounter, setValue(), getValue()
      using TimerKey = PhaseRoundDescription;
      // TimerKey must contain type PseudoKey and function pseudoKey() 
      using TimerValue = CounterPerPhaseTypeRawCommunicationValue;
      // TimerValue must contain functions setType() and setPseudoKeyCounter 

      struct OutputFormat {
        std::string prefix;
        std::string phase;
        size_t counterPerPhase;
        size_t round;
        std::string description;
        std::string type;
        bool rawCommunication;
        size_t value;

        friend std::ostream& operator<<(std::ostream& stream, const OutputFormat& outputFormat) {
          return stream << outputFormat.prefix 
            << " phase=" << outputFormat.phase 
            << " counterPerPhase=" << outputFormat.counterPerPhase
            << " round=" << outputFormat.round 
            << " operation=" << outputFormat.description
            << " type=" << outputFormat.type
            << " rawCommunication=" << outputFormat.rawCommunication
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
        disabled = false;
        prefix = "";
        curPhase = "none";
        curRound = 0;
      }

      void add(size_t value) {
        if (disabled)
          return;
        add(value, "unkown");
      }

      void add(size_t value, const std::string& description) {
        if (disabled)
          return;
        nonTimer.add(NonTimerRecord(curPhase, 0u, curRound, description, "number", false, value));
      }

      void addRawCommunication(size_t value, const std::string& description) {
        if (disabled)
          return;
        nonTimer.add(NonTimerRecord(curPhase, 0u, curRound, description, "number", true, value));
      }

      void start(const std::string& description) {
        if (disabled)
          return;
        timer.start(TimerKey(curPhase, curRound, description), TimerValue(0u, "", false, 0u));
      }

      void stop(const std::string& description) {
        if (disabled)
          return;
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
          data.push_back({prefix, 
              nonTimerRecord.phase, 
              nonTimerRecord.counterPerPhase, 
              nonTimerRecord.round, 
              nonTimerRecord.description, 
              nonTimerRecord.type, 
              nonTimerRecord.rawCommunication, 
              nonTimerRecord.value});

        for (const auto& [timerKey, timerValue] : timerRecords)
          data.push_back({prefix, 
              timerKey.phase, 
              timerValue.counterPerPhase, 
              timerKey.round, 
              timerKey.description, 
              timerValue.type, 
              timerValue.rawCommunication, 
              timerValue.value});
        return data;
      }

      void writeToStream(std::ostream& stream) {
        for (const auto& data : collect()) 
          stream << data << std::endl;
      }

      void eanble() {
        disabled = false;
      }
      
      void disable() {
        disabled = true;
      }

      void setPrefix(const std::string& prefix_) {
        prefix = prefix_;
      }

      void setPhase(const std::string& phase) {
        curPhase = phase;
      }

      void setRound(size_t round) {
        curRound = round;
      }

      private:
      bool disabled = false;
      std::string prefix = "";
      std::string curPhase = "none";
      size_t curRound = 0;
      NonTimer<NonTimerRecord> nonTimer;
      Timer<TimerKey, TimerValue> timer;
    };
  }
}
