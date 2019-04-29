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
#include <tlx/cmdline_parser.hpp>
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
    // using Generator = FileDistributer<StringSet>;
    using Generator = DNRatioGenerator<StringSet>;
    using Container = StringLcpContainer<StringSet>;

    using measurement::MeasuringTool;
    MeasuringTool& measuringTool = MeasuringTool::measuringTool();

    tlx::CmdlineParser cp;
    unsigned int numberOfIterations = 1;
    unsigned int numberOfStrings = 1000;
    double dToNRatio = 1;
    bool strongScaling = false;
    bool check = false;
    bool splitterMode = false;
    unsigned int samplingFactor = 2;
    unsigned int stringLength = 100;
    cp.add_flag('b', "splitterMode", splitterMode, " ");
    cp.add_unsigned('s', "size", numberOfStrings, "");
    cp.add_unsigned('i', "numberOfIterations", numberOfIterations, "");
    cp.add_unsigned('c', "samplingFactor", samplingFactor, "");
    cp.add_flag('x', "strongScaling", strongScaling, " ");
    cp.add_flag('y', "check", check, " ");
    cp.add_double('r', "dToNRatio", dToNRatio, "D/N ratio");
    cp.add_unsigned('a', "stringLength", stringLength, " string Length ");
    if (!cp.process(argc, argv)) {
        return -1;
    }
    // Initialize the MPI environment
    dss_schimek::mpi::environment env;
    MPI_Comm comm = env.communicator();
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
    for (size_t i = 0; i < numberOfIterations; ++i) {
        if (strongScaling) numberOfStrings *= size;
        if (splitterMode) numberOfStrings = env.size() * (env.size() - 1) * samplingFactor;
        measuringTool.setPrefix(
            "RESULT numberProcessors=" + std::to_string(size) + " iteration=" +
            std::to_string(i) + " size=" + std::to_string(numberOfStrings) +
            " strongScaling=" + std::to_string(strongScaling) +
            " dToNRatio=" + std::to_string(dToNRatio) +
            " stringLength=" + std::to_string(stringLength));

        // Create random input elements
        PRINT_ROOT("Create random input elements");
        // Container container = Generator("testData.dat");

        Container container = Generator(numberOfStrings, stringLength, dToNRatio);
        if (!container.isConsistent()) {
            std::cout << "initial input is corrupt" << std::endl;
            std::abort();
        }
        auto rawStrings = container.raw_strings();
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

        if (check) {
            measuringTool.start("allgatherv");
            auto globalUnsortedRawStrings = allgatherv(rawStrings);
            measuringTool.stop("allgatherv");
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
