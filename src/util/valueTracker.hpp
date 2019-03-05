#pragma once

#include <vector>
#include <map>
#include <iostream>

#include "mpi/allgather.hpp"

struct ValueIterationName {
  size_t value;
  size_t iteration;
  std::string name;
  ValueIterationName(const size_t value, const size_t iteration, const std::string& name)
    : value(value), iteration(iteration), name(name) {}
};

struct PhaseCounterPerPhaseDescriptionRoundValue {
  std::string phase;
  size_t counterPerPhase;
  size_t round;
  std::string description;
  size_t value;
  PhaseCounterPerPhaseDescriptionRoundValue(const std::string& phase, const size_t counterPerPhase, const std::string& description, const size_t round, const size_t value) 
    : phase(phase), counterPerPhase(counterPerPhase), description(description), round(round), value(value) {}

  friend std::ostream& operator<<(std::ostream& stream, PhaseCounterPerPhaseDescriptionRoundValue& elem) {
    return stream << "[" << elem.phase << ", " << elem.counterPerPhase << ", " << elem.round << ", " << elem.description << ", " << elem.value << "]";
  }
};

struct DescriptionCounterPerPhase {
  std::string description;
  size_t counterPerPhase;

  DescriptionCounterPerPhase(const std::string& description, const size_t counterPerPhase)
    : description(description), counterPerPhase(counterPerPhase) {}

  bool operator< (const DescriptionCounterPerPhase& rhs) const {
    if (description != rhs.description)
      return description < rhs.description;
    return counterPerPhase < rhs.counterPerPhase;
  }

  friend std::ostream& operator<<(std::ostream& stream, const DescriptionCounterPerPhase& elem) {
    return stream << "[" << elem.description <<  ", " << elem.counterPerPhase << "]";
  }
};

struct RoundValue {
  size_t round;
  size_t value;
  RoundValue() = default;
  RoundValue(const size_t round, const size_t value) :
    round(round), value(value) {}
};


//Output:
//Phase | counterPerPhase | round (set by user) | description | value 
//
//Part of Key, i.e. parts that must be unique: (Phase, counterPerPhase)


class ValueTracker {
  using Phase = std::string;
  using Key = DescriptionCounterPerPhase;
  using Value = RoundValue;
  using InternContainer = std::map<Key, Value>;

  public:
  using Result = PhaseCounterPerPhaseDescriptionRoundValue;

  ValueTracker() : curPhase("") {}

  static ValueTracker& valueTracker() {
    static ValueTracker valueTracker;
    return valueTracker;
  }

  void setPhase(const Phase& phase) { curPhase = phase; }

  void add(size_t value_) {
    size_t counterPerPhase = incrementCounterPerPhase(curPhase); 
    Key key("none", counterPerPhase);
    Value value(0u, value_);
    add(key, value);
  }

  void add(size_t value_, size_t round) {
    size_t counterPerPhase = incrementCounterPerPhase(curPhase); 
    Key key("none", counterPerPhase);
    Value value(round, value_);
    add(key, value);
  }

  void add(size_t value_, size_t round, const std::string& description) {
    size_t counterPerPhase = incrementCounterPerPhase(curPhase); 
    Key key(description, counterPerPhase);
    Value value(round, value_);
    add(key, value);
  }
  
  template<typename OutIt>
    void collect(OutIt out) {
      for (const auto& [key, value] : phaseToData) {
        collect(key, out);
      }
    }

  private:
  Phase curPhase;
  std::map<Phase, size_t> phaseToCounter;
  std::map<Phase, InternContainer> phaseToData;
  
  size_t incrementCounterPerPhase(const Phase& phase) {
    auto itCounter = phaseToCounter.find(curPhase);
    if (itCounter == phaseToCounter.end()) {
      phaseToCounter.emplace(curPhase, 1u);
      return 0;
    } else {
      return itCounter->second++;
    }
  }

  void add(const Key& key, const Value& value) {
    InternContainer& container = phaseToData[curPhase];
    if (container.find(key) != container.end()) {
      std::cout << "key: " << key << " already added" << std::endl;
      std::abort();
    } 
    container.emplace(key, value);
  }
  
  

  void writeToStream(std::ostream& stream,
                     const std::string& prefix, 
                     const Phase& phase,
                     const Key& key,
                     const Value& value) {
    stream << prefix  
      << " phase=" << phase 
      << " operation=" << key.description
      << " counterPerPhase=" << std::to_string(key.counterPerPhase)
      << " round=" << std::to_string(value.round)
      << " type=number"
      << " value=" << std::to_string(value.value)
      << std::endl;
  }

  template<typename OutIt>
    void collect(const Phase& phase, const Key& key, OutIt out) {
      auto itContainer = phaseToData.find(phase);
      if ( itContainer == phaseToData.end()) {
        std::cout << "Phase: " << phase << " not found" << std::endl;
        std::abort();
      }
      const InternContainer& container = (*itContainer).second;
      auto itValue = container.find(key);
      if (itValue == container.end()) {
        std::cout << "Key: " << key << " not found" << std::endl;
        std::abort();
      }

      Value value = (*itValue).second;
      std::vector<Value> values = dsss::mpi::allgather(value);
      for (const Value& value : values) {
        out = Result(phase, key.counterPerPhase, key.description, value.round, value.value);
      }
    }

  template<typename OutIt>
    void collect(const Phase& phase, OutIt out) {
      for (const auto& [key, value] : phaseToData[phase])
        collect(phase, key, out);
    }

  

};
