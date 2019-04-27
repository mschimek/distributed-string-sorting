/*****************************************************************************
 * This file is part of the Project JanusSortRBC
 *
 * Copyright (c) 2016-2017, Armin Wiebigke <armin.wiebigke@gmail.com>
 * Copyright (c) 2016-2019, Michael Axtmann <michael.axtmann@kit.edu>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#include <random>
#include <vector>

#include "JanusSort.hpp"
#include "RQuick.hpp"
#include "mpi/synchron.hpp"
#include "mpi/environment.hpp"
#include "strings/stringset.hpp"
#include "util/random_string_generator.hpp"

#define PRINT_ROOT(msg)                                                        \
    if (rank == 0) std::cout << msg << std::endl;

struct StringComparator {
    using String = dss_schimek::UCharLengthStringSet::String;
    bool operator()(String lhs, String rhs) {
        const unsigned char* lhsChars = lhs.string;
        const unsigned char* rhsChars = rhs.string;
        size_t counter = 0;
        //std::cout << "lhs: " << lhsChars << " rhs: " << rhsChars << std::endl;
        while (*lhsChars == *rhsChars && *lhsChars != 0) {
            ++lhsChars; ++rhsChars;
            counter++;
        }
        if (counter > 40) {
          std::cout << "attention!" << std::endl;
          std::abort();
        }
        return *lhsChars < *rhsChars;
    }
};

int main(int argc, char** argv) {
    using namespace dss_schimek;
    using namespace dss_schimek::mpi;
    using StringSet = UCharLengthStringSet;
    using Generator = DNRatioGenerator<StringSet>;
    using Container = StringLcpContainer<StringSet>;

    // Initialize the MPI environment
    dss_schimek::mpi::environment env;
    MPI_Comm comm = env.communicator();
    const uint64_t numStrings = 100;
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    // Create random input elements
    PRINT_ROOT("Create random input elements");
    Container container = Generator(numStrings);
    dss_schimek::mpi::execute_in_order([&](){

        std::cout << "rank: " << rank << std::endl;
    container.make_string_set().print();
        });
    std::mt19937_64 generator;
    int data_seed = 3469931 + rank;
    generator.seed(data_seed);
    std::uniform_real_distribution<double> dist(-100.0, 100.0);
    std::vector<double> data;
    for (int i = 0; i < 10; ++i)
        data.push_back(dist(generator));

    {
        /* RQuick */
        int tag = 11111;

        // Sort data descending

        auto data1 = data;
        StringComparator comp;
        PRINT_ROOT("Start sorting algorithm RQuick with MPI_Comm. "
                   << "RBC::Communicators are used internally.");
        container = RQuick::sort(generator, container.raw_strings(), MPI_BYTE,
            tag, comm, comp);
        env.barrier();
        dss_schimek::mpi::execute_in_order([&]() {
          std::cout << "rank: " << rank << " container size: " << container.size() << std::endl;
         container.make_string_set().print();
            });
        PRINT_ROOT("Elements have been sorted");

        // PRINT_ROOT("Start sorting algorithm RQuick with RBC::Comm.");
        // RBC::Comm rcomm;
        // RBC::Create_Comm_from_MPI(comm, &rcomm);
        // auto data2 = data;
        // RQuick::sort(generator, data2, MPI_DOUBLE, tag, rcomm,
        // std::greater<double>()); PRINT_ROOT("Elements have been sorted");

        // PRINT_ROOT("Start sorting algorithm RQuick with RBC::Comm. " <<
        //           "MPI communicators and MPI collectives are used.");
        // RBC::Comm rcomm1;
        // RBC::Create_Comm_from_MPI(comm, &rcomm1, true, true);
        // auto data3 = data;
        // RQuick::sort(generator, data3, MPI_DOUBLE, tag, rcomm1);
        // PRINT_ROOT("Elements have been sorted");
    }

    // Finalize the MPI environment
    // MPI_Finalize();
    env.finalize();
}
