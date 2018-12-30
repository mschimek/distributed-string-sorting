/*******************************************************************************
 * mpi/alltoall.hpp
 *
 * Copyright (C) 2018 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include <algorithm>
#include <cstdint>
#include <mpi.h>
#include <numeric>
#include <iterator>
#include <vector>
#include <cstring>

#include "mpi/allreduce.hpp"
#include "mpi/big_type.hpp"
#include "mpi/environment.hpp"
#include "mpi/type_mapper.hpp"
#include "mpi/scan.hpp"
#include "mpi/synchron.hpp"
#include "mpi/byte_encoder.hpp"

#include "util/indexed_string_set.hpp"
#include "util/string.hpp"
#include "util/string_set.hpp"

#include "strings/stringptr.hpp"
#include "strings/stringcontainer.hpp"
#include "strings/stringtools.hpp"

namespace dsss::mpi {

  static constexpr bool debug_alltoall = false;

  template <typename DataType>
    inline std::vector<DataType> alltoall(const std::vector<DataType>& send_data,
        environment env = environment()) {
      std::vector<DataType> receive_data(send_data.size(), 0);
      data_type_mapper<DataType> dtm;
      MPI_Alltoall(send_data.data(),
          send_data.size() / env.size(),
          dtm.get_mpi_type(),
          receive_data.data(),
          send_data.size() / env.size(),
          dtm.get_mpi_type(),
          env.communicator());
      return receive_data;
    }
  template <typename DataType>
    inline std::vector<DataType> alltoallv_small(
        const std::vector<DataType>& send_data, const std::vector<size_t>& send_counts,
        environment env = environment()) {

      std::vector<int32_t> real_send_counts(send_counts.size());
      for (size_t i = 0; i < send_counts.size(); ++i) {
        real_send_counts[i] = static_cast<int32_t>(send_counts[i]);
      }
      std::vector<int32_t> receive_counts = alltoall(real_send_counts, env);

      std::vector<int32_t> send_displacements(real_send_counts.size(), 0);
      std::vector<int32_t> receive_displacements(real_send_counts.size(), 0);
      for (size_t i = 1; i < real_send_counts.size(); ++i) {
        send_displacements[i] = send_displacements[i - 1] + real_send_counts[i - 1];
        receive_displacements[i] = receive_displacements[i - 1] +
          receive_counts[i - 1];
      }
      std::vector<DataType> receive_data(
          receive_counts.back() + receive_displacements.back());

      if constexpr (debug_alltoall) {
        for (int32_t i = 0; i < env.size(); ++i) {
          if (i == env.rank()) {
            std::cout << i << ": send_counts.size() " << send_counts.size()
              << std::endl;
            std::cout << i << ": send counts: ";
            for (const auto sc : real_send_counts) { std::cout << sc << ", "; }
            std::cout << std::endl << "receive counts: ";

            for (const auto rc : receive_counts) { std::cout << rc << ", "; }
            std::cout << std::endl;
          }
          env.barrier();
        }
      }

      data_type_mapper<DataType> dtm;
      MPI_Alltoallv(send_data.data(),
          real_send_counts.data(),
          send_displacements.data(),
          dtm.get_mpi_type(),
          receive_data.data(),
          receive_counts.data(),
          receive_displacements.data(),
          dtm.get_mpi_type(),
          env.communicator());
      return receive_data;
    }


  class AllToAllvSmall {
    public:
      static std::string getName()  {
        return "AlltoAllvSmall";
      }
      template<typename DataType>
        static std::vector<DataType> alltoallv(
            DataType* const send_data, const std::vector<size_t>& send_counts,
            environment env = environment()) {

          std::vector<int32_t> real_send_counts(send_counts.size());
          for (size_t i = 0; i < send_counts.size(); ++i) {
            real_send_counts[i] = static_cast<int32_t>(send_counts[i]);
          }
          std::vector<int32_t> receive_counts = alltoall(real_send_counts, env);

          std::vector<int32_t> send_displacements(real_send_counts.size(), 0);
          std::vector<int32_t> receive_displacements(real_send_counts.size(), 0);
          for (size_t i = 1; i < real_send_counts.size(); ++i) {
            send_displacements[i] = send_displacements[i - 1] + real_send_counts[i - 1];
            receive_displacements[i] = receive_displacements[i - 1] +
              receive_counts[i - 1];
          }
          std::vector<DataType> receive_data(
              receive_counts.back() + receive_displacements.back());

          if constexpr (debug_alltoall) {
            for (int32_t i = 0; i < env.size(); ++i) {
              if (i == env.rank()) {
                std::cout << i << ": send_counts.size() " << send_counts.size()
                  << std::endl;
                std::cout << i << ": send counts: ";
                for (const auto sc : real_send_counts) { std::cout << sc << ", "; }
                std::cout << std::endl << "receive counts: ";

                for (const auto rc : receive_counts) { std::cout << rc << ", "; }
                std::cout << std::endl;
              }
              env.barrier();
            }
          }

          data_type_mapper<DataType> dtm;
          MPI_Alltoallv(send_data,
              real_send_counts.data(),
              send_displacements.data(),
              dtm.get_mpi_type(),
              receive_data.data(),
              receive_counts.data(),
              receive_displacements.data(),
              dtm.get_mpi_type(),
              env.communicator());
          return receive_data;
        }

  };

  class AllToAllvDirectMessages {
    public:
      static std::string getName() {
        return "AlltoAllvDirectMessages";
      }
      template<typename DataType>
        static std::vector<DataType> alltoallv(DataType* const send_data,
            const std::vector<size_t>& send_counts, environment env = environment()) {

          std::vector<size_t> receive_counts = alltoall(send_counts, env);
          size_t local_receive_count = std::accumulate(
              receive_counts.begin(), receive_counts.end(), 0);

          std::vector<size_t> send_displacements(env.size(), 0);
          for (size_t i = 1; i < send_counts.size(); ++i) {
            send_displacements[i] = send_displacements[i - 1] + send_counts[i - 1];
          }
          std::vector<size_t> receive_displacements(env.size(), 0);
          for (size_t i = 1; i < send_counts.size(); ++i) {
            receive_displacements[i] =
              receive_displacements[i - 1] + receive_counts[i - 1];
          }

          std::vector<MPI_Request> mpi_request(2 * env.size());
          std::vector<DataType> receive_data(receive_displacements.back() +
              receive_counts.back());
          for (int32_t i = 0; i < env.size(); ++i) {
            // start with self send/recv
            auto source = (env.rank() + (env.size() - i)) % env.size();
            auto receive_type = get_big_type<DataType>(receive_counts[source]);
            MPI_Irecv(receive_data.data() + receive_displacements[source],
                1,
                receive_type,
                source,
                44227,
                env.communicator(),
                &mpi_request[source]);
          }
          // dispatch sends
          for (int32_t i = 0; i < env.size(); ++i) {
            auto target = (env.rank() + i) % env.size();
            auto send_type = get_big_type<DataType>(send_counts[target]);
            MPI_Isend(send_data + send_displacements[target],
                1,
                send_type,
                target,
                44227,
                env.communicator(),
                &mpi_request[env.size() + target]);
          }
          MPI_Waitall(2 * env.size(), mpi_request.data(), MPI_STATUSES_IGNORE);
          return receive_data;
        }
  };

  template<typename AllToAllvSmallPolicy>
    class AllToAllvCombined {
      public:
        static std::string getName() {
          return "AllToAllvCombined";
        }
        template<typename DataType>
          static std::vector<DataType> alltoallv(DataType* const send_data,
              const std::vector<size_t>& send_counts, environment env = environment()) {

            size_t local_send_count = std::accumulate(
                send_counts.begin(), send_counts.end(), 0);

            std::vector<size_t> receive_counts = alltoall(send_counts, env);
            size_t local_receive_count = std::accumulate(
                receive_counts.begin(), receive_counts.end(), 0);

            size_t local_max = std::max(local_send_count, local_receive_count);
            size_t global_max = allreduce_max(local_max, env);

            if (global_max < env.mpi_max_int()) {
              return AllToAllvSmallPolicy::alltoallv(send_data, send_counts, env);
            } else {
              std::vector<size_t> send_displacements(env.size(), 0);
              for (size_t i = 1; i < send_counts.size(); ++i) {
                send_displacements[i] = send_displacements[i - 1] + send_counts[i - 1];
              }
              std::vector<size_t> receive_displacements(env.size(), 0);
              for (size_t i = 1; i < send_counts.size(); ++i) {
                receive_displacements[i] =
                  receive_displacements[i - 1] + receive_counts[i - 1];
              }

              std::vector<MPI_Request> mpi_request(2 * env.size());
              std::vector<DataType> receive_data(receive_displacements.back() +
                  receive_counts.back());
              for (int32_t i = 0; i < env.size(); ++i) {
                // start with self send/recv
                auto source = (env.rank() + (env.size() - i)) % env.size();
                auto receive_type = get_big_type<DataType>(receive_counts[source]);
                MPI_Irecv(receive_data.data() + receive_displacements[source],
                    1,
                    receive_type,
                    source,
                    44227,
                    env.communicator(),
                    &mpi_request[source]);
              }
              // dispatch sends
              for (int32_t i = 0; i < env.size(); ++i) {
                auto target = (env.rank() + i) % env.size();
                auto send_type = get_big_type<DataType>(send_counts[target]);
                MPI_Isend(send_data + send_displacements[target],
                    1,
                    send_type,
                    target,
                    44227,
                    env.communicator(),
                    &mpi_request[env.size() + target]);
              }
              MPI_Waitall(2 * env.size(), mpi_request.data(), MPI_STATUSES_IGNORE);
              return receive_data;
            }
          }
    };
  //template <typename DataType>
  //    inline std::vector<DataType> alltoallv_small(
  //        DataType* send_data, std::vector<size_t>& send_counts,
  //        environment env = environment()) {
  //
  //      static size_t counter = 0;
  //      std::cout << "my small send rank: " << env.rank() << " counter = " << counter++ << std::endl;
  //      std::vector<int32_t> real_send_counts(send_counts.size());
  //      for (size_t i = 0; i < send_counts.size(); ++i) {
  //        real_send_counts[i] = static_cast<int32_t>(send_counts[i]);
  //      }
  //      std::vector<int32_t> receive_counts = alltoall(real_send_counts, env);
  //
  //      std::vector<int32_t> send_displacements(real_send_counts.size(), 0);
  //      std::vector<int32_t> receive_displacements(real_send_counts.size(), 0);
  //      for (size_t i = 1; i < real_send_counts.size(); ++i) {
  //        send_displacements[i] = send_displacements[i - 1] + real_send_counts[i - 1];
  //        receive_displacements[i] = receive_displacements[i - 1] +
  //          receive_counts[i - 1];
  //      }
  //      std::vector<DataType> receive_data(
  //          receive_counts.back() + receive_displacements.back());
  //
  //      if constexpr (debug_alltoall) {
  //        for (int32_t i = 0; i < env.size(); ++i) {
  //          if (i == env.rank()) {
  //            std::cout << i << ": send_counts.size() " << send_counts.size()
  //              << std::endl;
  //            std::cout << i << ": send counts: ";
  //            for (const auto sc : real_send_counts) { std::cout << sc << ", "; }
  //            std::cout << std::endl << "receive counts: ";
  //
  //            for (const auto rc : receive_counts) { std::cout << rc << ", "; }
  //            std::cout << std::endl;
  //          }
  //          env.barrier();
  //        }
  //      }
  //
  //      data_type_mapper<DataType> dtm;
  //      MPI_Alltoallv(send_data,
  //          real_send_counts.data(),
  //          send_displacements.data(),
  //          dtm.get_mpi_type(),
  //          receive_data.data(),
  //          receive_counts.data(),
  //          receive_displacements.data(),
  //          dtm.get_mpi_type(),
  //          env.communicator());
  //      return receive_data;
  //    }
  template <typename DataType>
    inline std::vector<DataType> alltoallv(std::vector<DataType>& send_data,
        const std::vector<size_t>& send_counts, environment env = environment()) {

      size_t local_send_count = std::accumulate(
          send_counts.begin(), send_counts.end(), 0);

      std::vector<size_t> receive_counts = alltoall(send_counts, env);
      size_t local_receive_count = std::accumulate(
          receive_counts.begin(), receive_counts.end(), 0);

      size_t local_max = std::max(local_send_count, local_receive_count);
      size_t global_max = allreduce_max(local_max, env);

      if (global_max < env.mpi_max_int()) {
        return alltoallv_small(send_data, send_counts, env);
      } else {
        std::vector<size_t> send_displacements(env.size(), 0);
        for (size_t i = 1; i < send_counts.size(); ++i) {
          send_displacements[i] = send_displacements[i - 1] + send_counts[i - 1];
        }
        std::vector<size_t> receive_displacements(env.size(), 0);
        for (size_t i = 1; i < send_counts.size(); ++i) {
          receive_displacements[i] =
            receive_displacements[i - 1] + receive_counts[i - 1];
        }

        std::vector<MPI_Request> mpi_request(2 * env.size());
        std::vector<DataType> receive_data(receive_displacements.back() +
            receive_counts.back());
        for (int32_t i = 0; i < env.size(); ++i) {
          // start with self send/recv
          auto source = (env.rank() + (env.size() - i)) % env.size();
          auto receive_type = get_big_type<DataType>(receive_counts[source]);
          MPI_Irecv(receive_data.data() + receive_displacements[source],
              1,
              receive_type,
              source,
              44227,
              env.communicator(),
              &mpi_request[source]);
        }
        // dispatch sends
        for (int32_t i = 0; i < env.size(); ++i) {
          auto target = (env.rank() + i) % env.size();
          auto send_type = get_big_type<DataType>(send_counts[target]);
          MPI_Isend(send_data.data() + send_displacements[target],
              1,
              send_type,
              target,
              44227,
              env.communicator(),
              &mpi_request[env.size() + target]);
        }
        MPI_Waitall(2 * env.size(), mpi_request.data(), MPI_STATUSES_IGNORE);
        return receive_data;
      }
    }

  //template <typename DataType>
  //  inline std::vector<DataType> alltoallv(DataType* send_data,
  //      std::vector<size_t>& send_counts, environment env = environment()) {

  //    size_t local_send_count = std::accumulate(
  //        send_counts.begin(), send_counts.end(), 0);

  //    std::vector<size_t> receive_counts = alltoall(send_counts, env);
  //    size_t local_receive_count = std::accumulate(
  //        receive_counts.begin(), receive_counts.end(), 0);

  //    size_t local_max = std::max(local_send_count, local_receive_count);
  //    size_t global_max = allreduce_max(local_max, env);

  //    if (global_max < env.mpi_max_int() && false) {
  //      return alltoallv_small(send_data, send_counts, env);
  //    } else {
  //      std::vector<size_t> send_displacements(env.size(), 0);
  //      for (size_t i = 1; i < send_counts.size(); ++i) {
  //        send_displacements[i] = send_displacements[i - 1] + send_counts[i - 1];
  //      }
  //      std::vector<size_t> receive_displacements(env.size(), 0);
  //      for (size_t i = 1; i < send_counts.size(); ++i) {
  //        receive_displacements[i] =
  //          receive_displacements[i - 1] + receive_counts[i - 1];
  //      }

  //      std::vector<MPI_Request> mpi_request(2 * env.size());
  //      std::vector<DataType> receive_data(receive_displacements.back() +
  //          receive_counts.back());
  //      for (int32_t i = 0; i < env.size(); ++i) {
  //        // start with self send/recv
  //        auto source = (env.rank() + (env.size() - i)) % env.size();
  //        auto receive_type = get_big_type<DataType>(receive_counts[source]);
  //        MPI_Irecv(receive_data.data() + receive_displacements[source],
  //            1,
  //            receive_type,
  //            source,
  //            44227,
  //            env.communicator(),
  //            &mpi_request[source]);
  //      }
  //      // dispatch sends
  //      for (int32_t i = 0; i < env.size(); ++i) {
  //        auto target = (env.rank() + i) % env.size();
  //        auto send_type = get_big_type<DataType>(send_counts[target]);
  //        MPI_Isend(send_data + send_displacements[target],
  //            1,
  //            send_type,
  //            target,
  //            44227,
  //            env.communicator(),
  //            &mpi_request[env.size() + target]);
  //      }
  //      MPI_Waitall(2 * env.size(), mpi_request.data(), MPI_STATUSES_IGNORE);
  //      return receive_data;
  //    }
  //  }

  inline std::vector<dsss::char_type> alltoallv_strings(
      dsss::string_set& send_data, const std::vector<size_t>& send_counts,
      environment env = environment()) {

    const size_t size = send_counts.size();
    std::vector<size_t> send_counts_char(size, 0);
    std::vector<dsss::char_type> send_buffer;
    send_buffer.reserve(send_data.data_container().size());
    // Determine the number of character that must be sended to each node
    for (size_t interval = 0, offset = 0; interval < size; ++interval) {
      // We cannot be sure that the pointers are ordered anymore (i.e., the memory
      // positions are not monotone Increasing)
      for (size_t j = offset; j < send_counts[interval] + offset; ++j) {
        const size_t string_length = dsss::string_length(send_data[j]) + 1;
        send_counts_char[interval] += string_length;
        std::copy_n(send_data[j], string_length, std::back_inserter(send_buffer));
      }
      offset += send_counts[interval];
    }

    if constexpr (debug_alltoall) {
      const size_t total_chars_sent = std::accumulate(
          send_counts_char.begin(), send_counts_char.end(), 0);
      const size_t total_chars_count = send_data.data_container().size();

      for (int32_t rank = 0; rank < env.size(); ++rank) {
        if (env.rank() == rank) {
          std::cout << rank << ": total_chars_sent " << total_chars_sent
            << ", total_chars_count: " << total_chars_count << std::endl;
        }
        env.barrier();
      }
    }
    return alltoallv(send_buffer, send_counts_char, env);
  }

  template <typename IndexType>
    inline dsss::indexed_string_set<IndexType> alltoallv_indexed_strings(
        dsss::indexed_string_set<IndexType>& send_data,
        std::vector<size_t>& send_counts_strings,
        environment env = environment()) {

      // Send the strings
      assert(send_counts_strings.size() == env.size());
      const size_t size = send_counts_strings.size();
      std::vector<dsss::char_type> real_send_data;
      std::vector<IndexType> index_send_data;
      std::vector<size_t> send_counts_char(size, 0);
      std::vector<size_t> send_displacements(size, 0);
      // Determine the number of character that must be sended to each node
      for (size_t scs_pos = 0, string_pos = 0; scs_pos < size; ++scs_pos) {
        for (size_t i = 0; i < send_counts_strings[scs_pos]; ++i) {
          index_send_data.emplace_back(send_data[string_pos].index);
          // The "+1" is there to also send the terminating 0
          const size_t string_length =
            dsss::string_length(send_data[string_pos].string) + 1;
          send_counts_char[scs_pos] += string_length;
          std::copy_n(send_data[string_pos].string, string_length,
              std::back_inserter(real_send_data));
          ++string_pos;
        }
      }

      std::vector<dsss::char_type> receive_data = alltoallv(real_send_data,
          send_counts_char,
          env);
      std::vector<IndexType> receive_data_indices = alltoallv(index_send_data,
          send_counts_strings,
          env);
      return dsss::indexed_string_set<IndexType>(std::move(receive_data),
          std::move(receive_data_indices));
    }

  template<typename StringSet, typename ByteEncoder>
    static inline std::vector<size_t> computeSendCountsBytes(
        const StringSet& ss,
        const std::vector<size_t>& intervals,
        const ByteEncoder& byteEncoder) {

      const auto begin = ss.begin();
      std::vector<size_t> sendCountsBytes(intervals.size(), 0);
      for (size_t interval = 0, offset = 0; interval < intervals.size(); ++interval) {
        size_t sendCountChars = 0;
        for (size_t j = offset; j < intervals[interval] + offset; ++j) {
          size_t stringLength = ss.get_length(ss[begin + j]); 
          sendCountChars += stringLength + 1;
        }
        sendCountsBytes[interval] = byteEncoder.computeNumberOfSendBytes(sendCountChars, 
            intervals[interval]); 
        offset += intervals[interval];
      }
      return sendCountsBytes;
    }

  template<typename StringSet, typename AllToAllPolicy, typename ByteEncoderPolicy, typename Timer>
    struct AllToAllStringImpl;

  template <typename StringSet, typename AllToAllPolicy, typename ByteEncoderPolicy, typename Timer>
    dss_schimek::StringLcpContainer<StringSet> alltoallv(
        dss_schimek::StringLcpContainer<StringSet>& container,
        const std::vector<size_t>& sendCountsString,
        Timer& timer,
        environment env = environment()) {
      static AllToAllStringImpl<StringSet, AllToAllPolicy, ByteEncoderPolicy, Timer> sender;
      return sender.alltoallv(container, sendCountsString, timer, env);

    }
  template<typename StringSet, typename AllToAllPolicy, typename ByteEncoderPolicy, typename Timer>
    struct AllToAllStringImpl : private ByteEncoderPolicy {
      dss_schimek::StringLcpContainer<StringSet> alltoallv(
          dss_schimek::StringLcpContainer<StringSet>& container,
          const std::vector<size_t>& sendCountsString,
          Timer& timer,
          environment env = environment()) {

        
        timer.start("all_to_all_strings_intern_copy", env);
        using namespace dss_schimek;
        using String = typename StringSet::String;
        using CharIterator = typename StringSet::CharIterator;

        const ByteEncoderPolicy byteEncoder;
        const StringSet& sendSet = container.make_string_set();

        std::vector<size_t> sendCountsTotal = 
          computeSendCountsBytes<StringSet,ByteEncoderPolicy>(sendSet, sendCountsString, byteEncoder);

        size_t totalNumberSendBytes = std::accumulate(sendCountsTotal.begin(), sendCountsTotal.end(), 0);

        std::vector<unsigned char> buffer(totalNumberSendBytes);
        unsigned char* curPos = buffer.data();
        for (size_t interval = 0, stringsWritten = 0; interval < sendCountsString.size(); ++interval) {
          auto begin = sendSet.begin() + stringsWritten;
          StringSet subSet = sendSet.sub(begin, begin + sendCountsString[interval]);
          curPos = byteEncoder.write(curPos, 
              subSet,
              container.lcp_array() + stringsWritten);
          stringsWritten += sendCountsString[interval];
        }
        timer.end("all_to_all_strings_intern_copy", env);

        timer.start("all_to_all_strings_mpi", env);
        std::vector<unsigned char> recv = AllToAllPolicy::alltoallv(buffer.data(), sendCountsTotal);
        timer.end("all_to_all_strings_mpi", env);
        timer.start("all_to_all_strings_read", env);
        auto [rawStrings, rawLcps] = byteEncoder.read(recv.data(), recv.size());
        timer.end("all_to_all_strings_read", env);
        return StringLcpContainer<StringSet>(std::move(rawStrings), std::move(rawLcps));
      }

    };

  template<typename StringSet, typename AllToAllPolicy, typename Timer> 
    struct AllToAllStringImpl<StringSet, AllToAllPolicy, dss_schimek::EmptyByteEncoderCopy, Timer> {
      dss_schimek::StringLcpContainer<StringSet> alltoallv(
          dss_schimek::StringLcpContainer<StringSet>& send_data,
          const std::vector<size_t>& send_counts,
          Timer& timer,
          environment env = environment()){

        using namespace dss_schimek;
        using String = typename StringSet::String;
        using CharIt = typename StringSet::CharIterator;
        timer.start("all_to_all_strings_intern_copy", env);
        const StringSet ss = send_data.make_string_set();

        if (send_data.size() == 0)
          return dss_schimek::StringLcpContainer<StringSet>();

        std::vector<unsigned char> receive_buffer_char;
        std::vector<size_t> receive_buffer_lcp;
        std::vector<size_t> send_counts_lcp(send_counts);
        std::vector<size_t> send_counts_char(send_counts.size(), 0);

        std::vector<unsigned char> send_buffer;
        send_buffer.reserve(send_data.char_size());
        for (size_t interval = 0, offset = 0; interval < send_counts.size(); ++interval) {
          for (size_t j = offset; j < send_counts[interval] + offset; ++j) {
            String str = ss[ss.begin() + j];
            size_t string_length = ss.get_length(str) + 1; 

            send_counts_char[interval] += string_length;
            std::copy_n(ss.get_chars(str, 0), string_length, std::back_inserter(send_buffer));
          }
          offset += send_counts[interval];
        }

        timer.end("all_to_all_strings_intern_copy", env);
        timer.start("all_to_all_strings_mpi", env);
        receive_buffer_char = AllToAllPolicy::alltoallv(send_buffer.data(), send_counts_char, env);
        receive_buffer_lcp = AllToAllPolicy::alltoallv(send_data.lcps().data(), send_counts, env);
        timer.end("all_to_all_strings_mpi", env);

        // no bytes are read in this version only for evaluation layout
        timer.start("all_to_all_strings_read", env);
        timer.end("all_to_all_strings_read", env);
        return dss_schimek::StringLcpContainer<StringSet>(
            std::move(receive_buffer_char), std::move(receive_buffer_lcp));
      }
    };
  
  template<typename StringSet, typename AllToAllPolicy, typename Timer> 
    struct AllToAllStringImpl<StringSet, AllToAllPolicy, dss_schimek::EmptyByteEncoderMemCpy, Timer> {
      dss_schimek::StringLcpContainer<StringSet> alltoallv(
          dss_schimek::StringLcpContainer<StringSet>& send_data,
          const std::vector<size_t>& sendCountsString,
          Timer& timer,
          environment env = environment()){

        using namespace dss_schimek;
        using String = typename StringSet::String;
        using CharIt = typename StringSet::CharIterator;
        timer.start("all_to_all_strings_intern_copy", env);
        const EmptyByteEncoderMemCpy byteEncoder;
        const StringSet ss = send_data.make_string_set();

        if (send_data.size() == 0)
          return dss_schimek::StringLcpContainer<StringSet>();

        std::vector<unsigned char> receive_buffer_char;
        std::vector<size_t> receive_buffer_lcp;
        std::vector<size_t> send_counts_lcp(sendCountsString);
        std::vector<size_t> send_counts_char(sendCountsString.size());

        std::vector<unsigned char> buffer(send_data.char_size());
        unsigned char* curPos = buffer.data();
        for (size_t interval = 0, stringsWritten = 0; interval < sendCountsString.size(); ++interval) {
          auto begin = ss.begin() + stringsWritten;
          StringSet subSet = ss.sub(begin, begin + sendCountsString[interval]);
          size_t numWrittenChars = 0;
          std::tie(curPos,  numWrittenChars)= byteEncoder.write(curPos, 
              subSet);

          send_counts_char[interval] = numWrittenChars;
          stringsWritten += sendCountsString[interval];
        }

        timer.end("all_to_all_strings_intern_copy", env);
        timer.start("all_to_all_strings_mpi", env);
        receive_buffer_char = AllToAllPolicy::alltoallv(buffer.data(), send_counts_char, env);
        receive_buffer_lcp = AllToAllPolicy::alltoallv(send_data.lcps().data(), sendCountsString, env);
        timer.end("all_to_all_strings_mpi", env);

        // no bytes are read in this version only for evaluation layout
        timer.start("all_to_all_strings_read", env);
        timer.end("all_to_all_strings_read", env);
        return dss_schimek::StringLcpContainer<StringSet>(
            std::move(receive_buffer_char), std::move(receive_buffer_lcp));
      }
    };

  template<typename StringSet, typename AllToAllPolicy, typename Timer> 
    struct AllToAllStringImpl<StringSet,
      AllToAllPolicy,
      dss_schimek::SequentialDelayedByteEncoder,
      Timer> {

      dss_schimek::StringLcpContainer<StringSet> alltoallv(
          dss_schimek::StringLcpContainer<StringSet>& container,
          const std::vector<size_t>& sendCountsString,
          Timer& timer,
          environment env = environment()){

        using namespace dss_schimek;
        using String = typename StringSet::String;
        using CharIterator = typename StringSet::CharIterator;

        const SequentialDelayedByteEncoder byteEncoder;

        timer.start("all_to_all_strings_intern_copy", env);
        const StringSet& sendSet = container.make_string_set();
        std::vector<unsigned char> contiguousStrings = 
          dss_schimek::getContiguousStrings(sendSet, container.char_size());

        const std::vector<size_t>& sendCountsLcp(sendCountsString);
        std::vector<size_t> sendCountsChar(sendCountsString.size(), 0);
        std::vector<size_t> sendCountsTotal(sendCountsString.size(), 0);

        for (size_t interval = 0, offset = 0; interval < sendCountsString.size(); ++interval) {
          for (size_t j = offset; j < sendCountsString[interval] + offset; ++j) {
            size_t stringLength = sendSet.get_length(sendSet[sendSet.begin() + j]); 
            sendCountsChar[interval] += stringLength + 1;
          }
          sendCountsTotal[interval] = byteEncoder.computeNumberOfSendBytes(sendCountsChar[interval], 
              sendCountsLcp[interval]); 
          offset += sendCountsString[interval];
        }

        size_t totalNumberSendBytes = 
          byteEncoder.computeNumberOfSendBytes(sendCountsChar, sendCountsLcp);

        std::vector<unsigned char> buffer(totalNumberSendBytes);
        unsigned char* curPos = buffer.data();

        for (size_t interval = 0, offset = 0, stringsWritten = 0; 
            interval < sendCountsString.size(); ++interval) {

          curPos = byteEncoder.write(curPos, 
              contiguousStrings.data() + offset, 
              sendCountsChar[interval], 
              container.lcp_array() + stringsWritten, 
              sendCountsLcp[interval]);
          offset += sendCountsChar[interval];
          stringsWritten += sendCountsLcp[interval];
        }

        timer.end("all_to_all_strings_intern_copy", env);
        timer.start("all_to_all_strings_mpi", env);
        std::vector<unsigned char> recv = AllToAllPolicy::alltoallv(buffer.data(), sendCountsTotal);
        timer.end("all_to_all_strings_mpi", env);
        timer.start("all_to_all_strings_read", env);
        auto [rawStrings, rawLcps] = byteEncoder.read(recv.data(), recv.size());
        timer.end("all_to_all_strings_read", env);
        return StringLcpContainer<StringSet>(std::move(rawStrings), std::move(rawLcps));
      }  
    };
} // namespace dsss::mpi
/******************************************************************************/
