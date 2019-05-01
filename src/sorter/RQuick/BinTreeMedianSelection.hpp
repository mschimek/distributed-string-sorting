/*****************************************************************************
 * This file is part of the Project Karlsruhe Distributed Sorting Library
 * (KaDiS).
 *
 * Copyright (c) 2019, Michael Axtmann <michael.axtmann@kit.edu>
 *
 * This program is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program.  If not, see <http://www.gnu.org/licenses/>.
 *******************************************************************************/

#pragma once

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <random>
#include <utility>
#include <vector>

#include "strings/stringcontainer.hpp"
#include "strings/stringset.hpp"

#include "mpi/synchron.hpp"

#include <tlx/math.hpp>

#include "./RandomBitStore.hpp"

int globalIteration;
int globalRank;
namespace BinTreeMedianSelection {
template <class Iterator, class Comp, class Communicator>
std::vector<unsigned char> select(Iterator begin, Iterator end, size_t n,
    Comp&& comp, MPI_Datatype mpi_type, std::mt19937_64& async_gen,
    RandomBitStore& bit_gen, int tag, Communicator& comm);

namespace _internal {
template <class T, class Comp>
void selectMedians(std::vector<T>& recv_v, std::vector<T>& tmp_v,
    std::vector<T>& vs, size_t n, Comp&& comp, std::mt19937_64& async_gen,
    RandomBitStore& bit_gen) {
    using StringSet = dss_schimek::UCharLengthStringSet;
    using StringContainer = dss_schimek::StringContainer<StringSet>;
    using String = StringSet::String;
    // std::copy(recv_v.begin(), recv_v.end(), std::back_inserter(vs));
    // recv_v.clear();
    // return;
    // assert(recv_v.size() <= n);
    // assert(vs.size() <= n);

    tmp_v.resize(recv_v.size() +
                 vs.size()); // recv_v and vs contain 0 terminated raw strings
    StringContainer recv_v_container(std::move(recv_v));
    StringContainer vs_container(std::move(vs));
    // recv_v_container.make_string_set().print();
    // vs_container.make_string_set().print();
    std::vector<String> tmp_v_strings(
        recv_v_container.size() + vs_container.size());

    std::merge(recv_v_container.getStrings().begin(),
        recv_v_container.getStrings().end(), vs_container.getStrings().begin(),
        vs_container.getStrings().end(), tmp_v_strings.begin(),
        std::forward<Comp>(comp));

    {
        uint64_t curPos = 0;
        for (const auto& str : tmp_v_strings) {
            std::copy(
                str.string, str.string + str.length + 1, tmp_v.data() + curPos);
            curPos += str.length + 1;
        }
    }
    auto copyTmp = tmp_v;
    StringContainer tmpContainer(std::move(copyTmp));
    // std::cout << "iteration: " << globalIteration << " on rank: " <<
    // globalRank << std::endl;
    // tmpContainer.make_string_set().print("iteration" +
    // std::to_string(globalIteration) + " on rank: " +
    // std::to_string(globalRank));

    if (tmp_v_strings.size() <= n) {
        vs = vs_container.releaseRawStrings();
        vs.swap(tmp_v);
        // assert(std::is_sorted(vs.begin(), vs.end(),
        // std::forward<Comp>(comp)));
        return;
    }
    else {
        auto ss = tmpContainer.make_string_set();
        if ((tmp_v_strings.size() - n) % 2 == 0) {
            const auto offset = (tmp_v_strings.size() - n) / 2;
            assert(offset + n < tmp_v_strings.size());
            vs.clear();

            for (size_t i = offset; i < offset + n; ++i) {
                const auto str = ss[ss.begin() + i];
                const auto length = ss.get_length(str) + 1;
                auto chars = ss.get_chars(str, 0);
                std::copy_n(chars, length, std::back_inserter(vs));
            }
            return;
        }
        else {
            // We cannot remove the same number of elements at
            // the right and left end.
            const auto offset = (tmp_v_strings.size() - n) / 2;
            const auto padding_cnt = bit_gen.getNextBit(async_gen);
            assert(padding_cnt <= 1);
            assert(offset + padding_cnt + n <= tmp_v_strings.size());

            vs.clear();
            for (size_t i = offset + padding_cnt; i < offset + padding_cnt + n;
                 ++i) {
                const auto str = ss[ss.begin() + i];
                const auto length = ss.get_length(str) + 1;
                auto chars = ss.get_chars(str, 0);
                std::copy_n(chars, length, std::back_inserter(vs));
            }
            return;
        }
    }
}

template <class T>
int64_t selectMedian(const std::vector<T>& v, std::mt19937_64& async_gen,
    RandomBitStore& bit_gen) {
    if (v.size() == 0) {
        return -1;
    }

    assert(v.size() > 0);
    if (v.size() % 2 == 0) {
        if (bit_gen.getNextBit(async_gen)) {
            return v.size() / 2;
        }
        else {
            return (v.size() / 2) - 1;
        }
    }
    else {
        return v.size() / 2;
    }
}
} // namespace _internal

template <class Iterator, class Comp, class Communicator>
std::vector<unsigned char> select(Iterator begin, Iterator end, size_t n,
    Comp&& comp, MPI_Datatype mpi_type, std::mt19937_64& async_gen,
    RandomBitStore& bit_gen, int tag, Communicator& comm) {
    using T = unsigned char;
    using StringSet = dss_schimek::UCharLengthStringSet;
    using StringContainer = dss_schimek::StringContainer<StringSet>;

    int32_t myrank, nprocs;
    MPI_Comm_rank(comm, &myrank);
    MPI_Comm_size(comm, &nprocs);

    // assert(static_cast<size_t>(end - begin) <= n);

    std::vector<T> v(begin, end);
    std::vector<T> recv_v;
    std::vector<T> tmp_v;

    // v.reserve(2 * n);
    // recv_v.reserve(2 * n);
    // tmp_v.reserve(2 * n);

    assert(std::is_sorted(begin, end, std::forward<Comp>(comp)));

    const auto tailing_zeros = static_cast<unsigned>(tlx::ffs(myrank)) - 1;
    const auto logp = tlx::integer_log2_floor(nprocs);
    const auto iterations = std::min(tailing_zeros, logp);

    for (size_t it = 0; it != iterations; ++it) {
        const auto source = myrank + (1 << it);

        MPI_Status status;
        MPI_Probe(source, tag, comm, &status);
        int count = 0;
        MPI_Get_count(&status, mpi_type, &count);
        assert(static_cast<size_t>(count) <= n);
        recv_v.resize(count);

        MPI_Recv(recv_v.data(), count, mpi_type, source, tag, comm,
            MPI_STATUS_IGNORE);

        globalIteration = it;
        globalRank = myrank;
        _internal::selectMedians(
            recv_v, tmp_v, v, n, std::forward<Comp>(comp), async_gen, bit_gen);
    }
    if (myrank == 0) {
        StringContainer container(std::move(v));
        auto medianIndex =
            _internal::selectMedian(container.getStrings(), async_gen, bit_gen);

        auto requestedRawString = container.getRawString(medianIndex);
        int32_t size = requestedRawString.size();
        MPI_Bcast(&size, 1, MPI_INT, 0, comm);
        MPI_Bcast(requestedRawString.data(), size, MPI_BYTE, 0, comm);
        return requestedRawString;
    }
    else {
        int target = myrank - (1 << tailing_zeros);
        assert(v.size() <= n);
        MPI_Send(v.data(), v.size(), MPI_BYTE, target, tag, comm);

        int32_t medianSize;
        MPI_Bcast(&medianSize, 1, MPI_INT, 0, comm);
        std::vector<unsigned char> median(medianSize);
        MPI_Bcast(median.data(), medianSize, MPI_CHAR, 0, comm);
        return median;
    }
}
} // namespace BinTreeMedianSelection
