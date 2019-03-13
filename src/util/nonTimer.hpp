#pragma once

#include "mpi/environment.hpp"
#include <vector>
#include <map>
#include <iostream>

#include "mpi/allgather.hpp"


template<typename Record>
class NonTimer {
  using PseudoKey = typename Record::PseudoKey;
  static const size_t allocationSize = 10000;
  public:

  NonTimer() {
    records.reserve(allocationSize);
  }

  template<typename OutIt>
    void collect(OutIt out) {
      for (Record& record : records) {
        collect(out, record);
      }
    }

  void add(const Record& record) {
    records.push_back(record);
    Record& storedRecord = records.back();
    storedRecord.setPseudoKeyCounter(incrementCounterPerPseudoKey(storedRecord));
  }

  private:
  dsss::mpi::environment env;
  std::vector<Record> records;
  std::map<PseudoKey, size_t> pseudoKeyToCounter;
  
  size_t incrementCounterPerPseudoKey(const Record& record) {
    auto itCounter = pseudoKeyToCounter.find(record.pseudoKey());
    if (itCounter == pseudoKeyToCounter.end()) {
      pseudoKeyToCounter.emplace(record.pseudoKey(), 1u);
      return 0;
    } else {
      return itCounter->second++;
    }
  }

    template<typename OutIt>
    inline void collect(OutIt out, Record& record) {
      // Cannot allgather Record as a whole since it is not trivially_copyable
      // -> workaround send only values
      auto ownValue = record.getValue();
      auto globalValues = dsss::mpi::allgather(ownValue);
      const size_t sum = std::accumulate(globalValues.begin(), globalValues.end(), 0);
      //for (auto& globalValue: globalValues) {
        record.setValue(sum);
        out = record;
      //}
    }
};
