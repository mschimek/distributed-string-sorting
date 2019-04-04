#pragma once 

#include <fstream>
#include <iostream>
#include <mpi.h>
#include <numeric>
#include <vector>

#include "mpi/big_type.hpp"
#include "mpi/environment.hpp"
#include "mpi/type_mapper.hpp"

namespace dss_schimek {

size_t getFileSize(const std::string& path) {
    std::ifstream in(path, std::ifstream::ate | std::ifstream::binary);
    if (!in.good()) {
        std::cout << " file not good" << std::endl;
        std::abort();
    }
    const size_t fileSize = in.tellg();
    in.close();
    return fileSize;
}

struct RawStringsLines {
    std::vector<unsigned char> rawStrings;
    size_t lines;
};

dss_schimek::RawStringsLines readFile(const std::string& path) {
    using dss_schimek::RawStringsLines;
    RawStringsLines data;
    const size_t fileSize = getFileSize(path);
    std::ifstream in(path);
    std::vector<unsigned char>& rawStrings = data.rawStrings;
    rawStrings.reserve(1.5 * fileSize);

    std::string line;
    data.lines = 0u;
    while (std::getline(in, line)) {
        ++data.lines;
        for (unsigned char curChar : line)
            rawStrings.push_back(curChar);
        rawStrings.push_back(0);
    }
    std::cout << "lines: " << data.lines << std::endl;
    in.close();
    return data;
}

std::vector<unsigned char> readFileAndDistribute(const std::string& path) {
    dss_schimek::mpi::environment env;
    std::vector<unsigned char> rawStrings;
    size_t recvCount = 0;
    dss_schimek::mpi::data_type_mapper<size_t> dtm;
    std::vector<size_t> sendCounts;
    if (env.rank() == 0) {
        dss_schimek::RawStringsLines data = readFile(path);
        rawStrings = std::move(data.rawStrings);
        std::cout << "lines in readFileAndDistribute: " << data.lines
                  << std::endl;
        const size_t smallPacketSize = data.lines / env.size();
        const size_t bigPacketSize =
            smallPacketSize + (data.lines % env.size());
        size_t curPos = 0;
        size_t lastPos = 0;
        for (size_t i = 0; i + 1 < env.size(); ++i) {
            for (size_t numStrings = 0; numStrings < smallPacketSize;
                 ++numStrings) {
                while (rawStrings[curPos] != 0)
                    ++curPos;
                ++curPos;
            }
            sendCounts.push_back(curPos - lastPos);
            lastPos = curPos;
        }
        std::cout << "rawStringsize: " << rawStrings.size()
                  << " lastPos: " << lastPos << std::endl;
        sendCounts.push_back(rawStrings.size() - lastPos);
        std::cout << "smallPacketSize: " << smallPacketSize
                  << " bigPacketSize: " << bigPacketSize << std::endl;
    }
    MPI_Scatter(sendCounts.data(), 1, dtm.get_mpi_type(), &recvCount, 1,
        dtm.get_mpi_type(), 0, env.communicator());

    std::vector<unsigned char> localRawStrings(recvCount);
    if (env.rank() == 0) {
        std::vector<size_t> offsets;
        offsets.push_back(0);
        std::partial_sum(
            sendCounts.begin(), sendCounts.end(), std::back_inserter(offsets));

        std::copy(rawStrings.begin(), rawStrings.begin() + sendCounts[0],
            localRawStrings.begin());
        for (int32_t i = 1; i < static_cast<int32_t>(env.size()); ++i) {
            auto receive_type =
                dss_schimek::mpi::get_big_type<unsigned char>(sendCounts[i]);
            MPI_Send(rawStrings.data() + offsets[i], 1, receive_type, i, 42,
                env.communicator());
        }
    }
    else {
        auto receive_type = dss_schimek::mpi::get_big_type<unsigned char>(recvCount);
        MPI_Recv(localRawStrings.data(), 1, receive_type, 0, 42,
            env.communicator(), MPI_STATUSES_IGNORE);
    }

    std::cout << "rank: " << env.rank() << " recvCount: " << recvCount
              << std::endl;
    return localRawStrings;
}

} // namespace dss_schimek
