#pragma once

#include "mpi/environment.hpp"
#include "util/nonTimer.hpp"
#include "util/timer.hpp"
#include <ostream>
#include <vector>

namespace dss_schimek {
namespace measurement {

struct PhaseCounterPerPhaseRoundDescriptionTypeRawCommunicationSumUpValue {
    using PseudoKey = std::string;

    std::string phase;
    size_t counterPerPhase;
    size_t round;
    std::string description;
    std::string type;
    bool rawCommunication;
    bool sumUp;
    size_t value;

    PhaseCounterPerPhaseRoundDescriptionTypeRawCommunicationSumUpValue(
        const std::string& phase, const size_t counterPerPhase,
        const size_t round, const std::string& description,
        const std::string& type, const bool rawCommunication, const bool sumUp,
        const size_t value)
        : phase(phase), counterPerPhase(counterPerPhase), round(round),
          description(description), type(type),
          rawCommunication(rawCommunication), sumUp(sumUp), value(value) {}

    const std::string& pseudoKey() const { return phase; }

    void setPseudoKeyCounter(size_t counter) { counterPerPhase = counter; }

    void setType(const std::string& type_) { type = type_; }

    void setValue(size_t value_) { value = value_; }

    size_t getValue() const { return value; }

    bool getSumUp() { return sumUp; }

    friend std::ostream& operator<<(std::ostream& stream,
        PhaseCounterPerPhaseRoundDescriptionTypeRawCommunicationSumUpValue&
            elem) {
        return stream << "[" << elem.phase << ", " << elem.counterPerPhase
                      << ", " << elem.round << ", " << elem.description << ", "
                      << elem.rawCommunication << ", " << elem.sumUp
                      << elem.value << "]";
    }
};

struct PhaseRoundDescription {
    using PseudoKey = std::string;

    std::string phase;
    size_t round;
    std::string description;

    PhaseRoundDescription(const std::string& phase, const size_t round,
        const std::string& description)
        : phase(phase), round(round), description(description) {}

    const std::string& pseudoKey() const { return phase; }

    bool operator<(const PhaseRoundDescription& rhs) const {
        return std::tie(phase, round, description) <
               std::tie(rhs.phase, rhs.round, rhs.description);
    }

    friend std::ostream& operator<<(
        std::ostream& stream, const PhaseRoundDescription& elem) {
        return stream << "[" << elem.phase << ", " << elem.round << ", "
                      << elem.description << "]";
    }
};

struct CounterPerPhaseTypeRawCommunicationSumUpValue {
    size_t counterPerPhase;
    std::string type;
    bool rawCommunication;
    bool sumUp;
    size_t value;
    CounterPerPhaseTypeRawCommunicationSumUpValue(const size_t counterPerPhase,
        const std::string& type, const bool rawCommunication, const bool sumUp,
        const size_t value)
        : counterPerPhase(counterPerPhase), type(type),
          rawCommunication(rawCommunication), sumUp(sumUp), value(value) {}

    void setPseudoKeyCounter(size_t counter) { counterPerPhase = counter; }

    void setType(const std::string& type_) { type = type_; }

    void setValue(size_t value_) { value = value_; }
};

class MeasuringTool {
    // Columns:
    // Prefix | Phase | CounterPerPhase | Round | Description | Type |
    // RawCommunication | SumUp |  Value
    //
    // NonTimer: Duplicates are allowed (but since counterPerPhase is
    // incremented there are no duplicates if all fields are taken into account)
    // Timer: Key = (Phase, Round, Description) Value = (CounterPerPhase, Type,
    // RawCommunication, SumUp, Value)

    using NonTimerRecord =
        PhaseCounterPerPhaseRoundDescriptionTypeRawCommunicationSumUpValue;
    // NonTimerRecord must contain type PseudoKey and functions pseudoKey() and
    // setPseudoKeyCounter, setValue(), getValue()
    using TimerKey = PhaseRoundDescription;
    // TimerKey must contain type PseudoKey and function pseudoKey()
    using TimerValue = CounterPerPhaseTypeRawCommunicationSumUpValue;
    // TimerValue must contain functions setType() and setPseudoKeyCounter
    //
    dsss::mpi::environment env;

    struct OutputFormat {
        std::string prefix;
        std::string phase;
        size_t counterPerPhase;
        size_t round;
        std::string description;
        std::string type;
        bool rawCommunication;
        bool sumUp;
        size_t value;

        friend std::ostream& operator<<(
            std::ostream& stream, const OutputFormat& outputFormat) {
            return stream << outputFormat.prefix
                          << " phase=" << outputFormat.phase
                          << " counterPerPhase=" << outputFormat.counterPerPhase
                          << " round=" << outputFormat.round
                          << " operation=" << outputFormat.description
                          << " type=" << outputFormat.type
                          << " rawCommunication="
                          << outputFormat.rawCommunication
                          << " sumUp=" << outputFormat.sumUp
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
        verbose = false;
        prefix = "";
        curPhase = "none";
        curRound = 0;
    }

    void add(size_t value) {
        if (disabled) return;
        add(value, "unkown");
    }

    void add(
        size_t value, const std::string& description, const bool sumUp = true) {
        if (disabled) return;
        if (verbose && env.rank() == 0) std::cout << description << std::endl;
        nonTimer.add(NonTimerRecord(curPhase, 0u, curRound, description,
            "number", false, sumUp, value));
    }

    void addRawCommunication(size_t value, const std::string& description) {
        if (disabled) return;
        if (verbose && env.rank() == 0) std::cout << description << std::endl;
        nonTimer.add(NonTimerRecord(curPhase, 0u, curRound, description,
            "number", true, "true", value));
    }

    void start(const std::string& description) {
        if (disabled) return;
        if (verbose && env.rank() == 0) std::cout << description << std::endl;
        timer.start(TimerKey(curPhase, curRound, description),
            TimerValue(0u, "", false, false, 0u));
    }

    void stop(const std::string& description) {
        if (disabled) return;
        if (verbose && env.rank() == 0) std::cout << description << std::endl;
        timer.stop(TimerKey(curPhase, curRound, description));
    }

    std::vector<OutputFormat> collect() {
        disable();
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
            data.push_back(
                {prefix, nonTimerRecord.phase, nonTimerRecord.counterPerPhase,
                    nonTimerRecord.round, nonTimerRecord.description,
                    nonTimerRecord.type, nonTimerRecord.rawCommunication,
                    nonTimerRecord.sumUp, nonTimerRecord.value});

        for (const auto& [timerKey, timerValue] : timerRecords)
            data.push_back({prefix, timerKey.phase, timerValue.counterPerPhase,
                timerKey.round, timerKey.description, timerValue.type,
                timerValue.rawCommunication, timerValue.sumUp,
                timerValue.value});
        enable();
        return data;
    }

    void writeToStream(std::ostream& stream) {
        for (const auto& data : collect())
            stream << data << std::endl;
    }

    void enable() { disabled = false; }

    void disable() { disabled = true; }

    void setVerbose(const bool value) { verbose = value; }

    void setPrefix(const std::string& prefix_) { prefix = prefix_; }

    void setPhase(const std::string& phase) { curPhase = phase; }

    void setRound(size_t round) { curRound = round; }

private:
    bool disabled = false;
    bool verbose = false;
    std::string prefix = "";
    std::string curPhase = "none";
    size_t curRound = 0;
    NonTimer<NonTimerRecord> nonTimer;
    Timer<TimerKey, TimerValue> timer;
};
} // namespace measurement
} // namespace dss_schimek
