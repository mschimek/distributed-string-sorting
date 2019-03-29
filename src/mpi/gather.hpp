#pragma once

#include "mpi/big_type.hpp"
#include "mpi/environment.hpp"
#include "mpi/type_mapper.hpp"
#include "util/measuringTool.hpp"
#include "util/string_set.hpp"

namespace dss_schimek {
namespace mpi {
template <typename DataType>
inline std::vector<DataType> gather(const DataType& send_data, int32_t root,
    dsss::mpi::environment env = dsss::mpi::environment()) {
    using namespace dsss::mpi;

    using dss_schimek::measurement::MeasuringTool;

    MeasuringTool& measuringTool = MeasuringTool::measuringTool();
    measuringTool.addRawCommunication(sizeof(DataType), "gather");

    data_type_mapper<DataType> dtm;
    std::vector<DataType> receive_data;
    if (env.rank() == root) receive_data.resize(env.size());
    MPI_Gather(&send_data, 1, dtm.get_mpi_type(), receive_data.data(), 1,
        dtm.get_mpi_type(), root, env.communicator());
    return receive_data;
}

template <typename DataType>
inline std::vector<DataType> gatherv(std::vector<DataType>& send_data,
    int32_t root, dsss::mpi::environment env = dsss::mpi::environment()) {
    using namespace dsss::mpi;
    using namespace dss_schimek;

    using dss_schimek::measurement::MeasuringTool;

    MeasuringTool& measuringTool = MeasuringTool::measuringTool();
    measuringTool.addRawCommunication(sizeof(DataType), "gatherv");

    std::vector<size_t> receiveCounts = gather(send_data.size(), root);
    std::vector<size_t> offsets;
    offsets.reserve(receiveCounts.size() + 1);
    offsets.push_back(0);
    std::partial_sum(receiveCounts.begin(), receiveCounts.end(),
        std::back_inserter(offsets));

    data_type_mapper<DataType> dtm;
    std::vector<DataType> receive_data;
    if (env.rank() == root) {
        receive_data.resize(offsets.back());
        std::copy(send_data.begin(), send_data.end(),
            receive_data.begin() + offsets[root]);
        for (int32_t i = root + 1; i < root + env.size(); ++i) {
            int32_t partner = i % env.size();
            MPI_Recv(receive_data.data() + offsets[partner],
                receiveCounts[partner], dtm.get_mpi_type(), partner, 42,
                env.communicator(), MPI_STATUSES_IGNORE);
        }
    } else {
      MPI_Send(send_data.data(), send_data.size(), dtm.get_mpi_type(), root, 42, env.communicator());
    }
    return receive_data;
}

} // namespace mpi
} // namespace dss_schimek

