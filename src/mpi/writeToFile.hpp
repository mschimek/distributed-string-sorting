#pragma once

#include "mpi/environment.hpp"
#include <string>
#include <vector>

namespace dss_schimek::mpi {
void writeToOwnFile(std::string filename, std::vector<unsigned char>& content,
    dss_schimek::mpi::environment env = dss_schimek::mpi::environment()) {
    MPI_File fh;
    filename += std::to_string(env.rank());
    MPI_File_open(MPI_COMM_SELF, filename.c_str(),
        MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
    MPI_File_write(
        fh, content.data(), content.size(), MPI_BYTE, MPI_STATUS_IGNORE);
    MPI_File_close(&fh);
}
} // namespace dss_schimek::mpi
