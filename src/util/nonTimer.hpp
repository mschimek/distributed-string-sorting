#pragma once

#include "mpi/environment.hpp"
#include "mpi/type_mapper.hpp"
#include <algorithm>
#include <iostream>
#include <map>
#include <numeric>
#include <vector>

template <typename DataType>
inline std::vector<DataType> allgather_(DataType& send_data,
    dsss::mpi::environment env = dsss::mpi::environment()) {

    dsss::mpi::data_type_mapper<DataType> dtm;
    std::vector<DataType> receive_data(env.size());
    MPI_Allgather(&send_data, 1, dtm.get_mpi_type(), receive_data.data(), 1,
        dtm.get_mpi_type(), env.communicator());
    return receive_data;
}

template <typename Record>
class NonTimer {
    using PseudoKey = typename Record::PseudoKey;
    static const size_t allocationSize = 10000;

public:
    NonTimer() { records.reserve(allocationSize); }

    template <typename OutIt>
    void collect(OutIt out) {
        for (Record& record : records) {
            collect(out, record);
        }
    }

    void add(const Record& record) {
        records.push_back(record);
        Record& storedRecord = records.back();
        storedRecord.setPseudoKeyCounter(
            incrementCounterPerPseudoKey(storedRecord));
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
        }
        else {
            return itCounter->second++;
        }
    }

    template <typename OutIt>
    inline void collect(OutIt out, Record& record) {
        // Cannot allgather Record as a whole since it is not trivially_copyable
        // -> workaround send only values
        auto ownValue = record.getValue();
        auto globalValues = allgather_(ownValue);
        const size_t sum = std::accumulate(
            globalValues.begin(), globalValues.end(), static_cast<size_t>(0u));
        // for (auto& globalValue: globalValues) {
        record.setValue(sum);
        out = record;
        //}
    }
};
