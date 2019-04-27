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
#include "mpi/environment.hpp"
#include "mpi/synchron.hpp"
#include "strings/stringset.hpp"
#include "util/measuringTool.hpp"
#include "util/random_string_generator.hpp"
#include <tlx/sort/strings/radix_sort.hpp>

#define PRINT_ROOT(msg)                                                        \
    if (rank == 0) std::cout << msg << std::endl;

struct StringComparator {
    using String = dss_schimek::UCharLengthStringSet::String;
    bool operator()(String lhs, String rhs) {
        const unsigned char* lhsChars = lhs.string;
        const unsigned char* rhsChars = rhs.string;
        size_t counter = 0;
        // std::cout << "lhs: " << lhsChars << " rhs: " << rhsChars <<
        // std::endl;
        while (*lhsChars == *rhsChars && *lhsChars != 0) {
            ++lhsChars;
            ++rhsChars;
            counter++;
        }
        return *lhsChars < *rhsChars;
    }
};

int main(int argc, char** argv) {
    using namespace dss_schimek;
    using namespace dss_schimek::mpi;
    using StringSet = UCharLengthStringSet;
    using Generator = FileDistributer<StringSet>;
    // using Generator = DNRatioGenerator<StringSet>;
    using Container = StringLcpContainer<StringSet>;

    using measurement::MeasuringTool;
    MeasuringTool& measuringTool = MeasuringTool::measuringTool();

    // Initialize the MPI environment
    dss_schimek::mpi::environment env;
    MPI_Comm comm = env.communicator();
    const uint64_t numStrings = 100000;
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
    for (size_t i = 0; i < 5; ++i)
    {
      measuringTool.setPrefix("iteration=" + std::to_string(i));
    // Create random input elements
    PRINT_ROOT("Create random input elements");
    Container container = Generator("testData.dat");
    // Container container = Generator(numStrings);
    if (!container.isConsistent()) {
        std::cout << "initial input is corrupt" << std::endl;
        std::abort();
    }
    auto rawStrings = container.raw_strings();
    measuringTool.start("allgatherv");
    auto globalUnsortedRawStrings = allgatherv(rawStrings);
    measuringTool.stop("allgatherv");
    std::mt19937_64 generator;
    int data_seed = 3469931 + rank;
    generator.seed(data_seed);
     
        /* RQuick */
        int tag = 11111;

        // Sort data descending

        StringComparator comp;
        PRINT_ROOT("Start sorting algorithm RQuick with MPI_Comm. "
                   << "RBC::Communicators are used internally.");
        measuringTool.start("distributed_sort");
        auto sortedContainer = RQuick::sort(
            generator, container.raw_strings(), MPI_BYTE, tag, comm, comp);
        measuringTool.stop("distributed_sort");
        env.barrier();
        std::vector<unsigned char> sortedRawStrings(container.char_size());

        env.barrier();
        sortedContainer.orderRawStrings();
        std::vector<unsigned char> localSortedRawStrings =
            sortedContainer.raw_strings();
        const std::vector<unsigned char> globalSortedRawStrings =
            allgatherv(localSortedRawStrings);
        StringLcpContainer<StringSet> globalUnsortedStrings(
            std::move(globalUnsortedRawStrings));
        auto stringPtr = globalUnsortedStrings.make_string_ptr();
        measuringTool.start("local_sort");
        tlx::sort_strings_detail::radixsort_CI3(stringPtr, 0, 0);
        measuringTool.stop("local_sort");
        globalUnsortedStrings.orderRawStrings();

        auto finalInitRawStrings = globalUnsortedStrings.raw_strings();
        if (finalInitRawStrings.size() != globalSortedRawStrings.size()) {
            std::cout << "we have lost chars: " << std::endl;
            std::abort();
        }
        if (finalInitRawStrings != globalSortedRawStrings) {
            std::cout << "something is wrong " << std::endl;
            std::abort();
        }

        PRINT_ROOT("Elements have been sorted");
        std::stringstream buffer;
        measuringTool.writeToStream(buffer);
        if (env.rank() == 0) {
            std::cout << buffer.str() << std::endl;
        }
        measuringTool.reset();

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
