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

#include "util/indexed_string_set.hpp"
#include "util/string.hpp"
#include "util/string_set.hpp"

#include "strings/stringptr.hpp"
#include "strings/stringcontainer.hpp"
#include "strings/stringtools.hpp"

namespace dsss::mpi {

static constexpr bool debug_alltoall = false;

template <typename DataType>
inline std::vector<DataType> alltoall(std::vector<DataType>& send_data,
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
  std::vector<DataType>& send_data, std::vector<size_t>& send_counts,
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

template <typename DataType>
inline std::vector<DataType> alltoallv_small(
  DataType* send_data, std::vector<size_t>& send_counts,
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

template <typename DataType>
inline std::vector<DataType> alltoallv(std::vector<DataType>& send_data,
    std::vector<size_t>& send_counts, environment env = environment()) {

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
    std::vector<size_t> send_displacements(0, env.size());
    for (size_t i = 1; i < send_counts.size(); ++i) {
      send_displacements[i] = send_displacements[i - 1] + send_counts[i - 1];
    }
    std::vector<size_t> receive_displacements(0, env.size());
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

template <typename StringSet>
dss_schimek::StringLcpContainer<StringSet> alltoallv(
    dss_schimek::StringLcpContainer<StringSet>& send_data,
    const std::vector<size_t>& send_counts,
    environment env = environment()){

  using namespace dss_schimek;
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
      size_t string_length; 
      if constexpr(std::is_same<StringSet, UCharLengthStringSet>::value)
       string_length = dss_schimek::string_length(send_data[j].string) + 1;
      else
        string_length = dss_schimek::string_length(send_data[j]) + 1;
      send_counts_char[interval] += string_length;
      if constexpr(std::is_same<StringSet, UCharLengthStringSet>::value)
      std::copy_n(send_data[j].string, string_length, std::back_inserter(send_buffer));
      else
      std::copy_n(send_data[j], string_length, std::back_inserter(send_buffer));
    }
    offset += send_counts[interval];
  }
  //for (size_t i = 0; i < send_counts.size(); ++i)
  //  std::cout << i << " counts: " << send_counts_char[i] << std::endl;

  receive_buffer_char = alltoallv(send_buffer, send_counts_char, env);
  receive_buffer_lcp = alltoallv(send_data.lcps(), send_counts_lcp, env);

  return dss_schimek::StringLcpContainer<StringSet>(
      std::move(receive_buffer_char), std::move(receive_buffer_lcp));
}

class ByteEncoder {
  /*
   * Encoding: 
   * (|Number Chars : size_t | Number Strings: size_t | [char : unsigned char] | [number : size_t] |)*
   *
   */
  public:
    size_t computeNumberOfSendBytes(size_t charsToSend, size_t numbersToSend) const {
      return charsToSend + sizeof(size_t) * numbersToSend + 2 * sizeof(size_t);
    }

    size_t computeNumberOfSendBytes(const std::vector<size_t>& charsToSend,
        const std::vector<size_t>& numbersToSend) const {
      assert(charsToSend.size() == numbersToSend.size());
      size_t numberOfSendBytes = 0;
      for (size_t i = 0; i < charsToSend.size(); ++i) {
        numberOfSendBytes += computeNumberOfSendBytes(charsToSend[i], numbersToSend[i]); 
      }
      return numberOfSendBytes;
    }

    std::pair<size_t, size_t> computeNumberOfRecvData(const char unsigned* buffer, size_t size) const {
     size_t recvChars = 0, recvNumbers  = 0;
     size_t i = 0; 
     while(i < size) {
       size_t curRecvChars = 0;
       size_t curRecvNumbers = 0;
       memcpy(&curRecvChars, buffer + i, sizeof(size_t));
       recvChars += curRecvChars;
       i += sizeof(size_t);
       memcpy(&curRecvNumbers, buffer + i, sizeof(size_t));
       recvNumbers += curRecvNumbers;
       i += sizeof(size_t);
       i += curRecvChars + sizeof(size_t) * curRecvNumbers;
     }
     return std::make_pair(recvChars, recvNumbers); 
    }

    // start work here
    template<typename StringSet>
    unsigned char* write(unsigned char* buffer,
        const StringSet ss,
        const size_t* numbersToWrite) const {

      using String = typename StringSet::String;
      using CharIt = typename StringSet::CharIterator;

      unsigned char* const startOfBuffer = buffer;
      size_t numChars = 0;
      const size_t size = ss.size();

      std::cout << "size: " << size << std::endl;
      buffer += sizeof(size_t);
      memcpy(buffer, &size, sizeof(size_t));
      buffer += sizeof(size_t);
      auto beginOfSet = ss.begin();
      for (size_t i = 0; i < size; ++i) {
        String str = ss[beginOfSet + i];
        size_t stringLength = ss.get_length(str) + 1;
        numChars += stringLength;
        memcpy(buffer, ss.get_chars(str, 0), stringLength);
        buffer += stringLength;
      }
      memcpy(startOfBuffer, &numChars, sizeof(size_t));
      memcpy(buffer, numbersToWrite, size * sizeof(size_t));
      buffer += size * sizeof(size_t);
      dsss::mpi::environment env;
      size_t output1 = 0;
      size_t output2 = 0; 
      memcpy(&output1, startOfBuffer, sizeof(size_t));
      memcpy(&output2, startOfBuffer + sizeof(size_t), sizeof(size_t));

      std::cout << "rank: " << env.rank() << " numChars:" << output1 << " numStrings " << output2 << std::endl;
      return buffer;
    }

    unsigned char*  write(unsigned char* buffer,
        const unsigned char* charsToWrite,
        size_t numChars,
        const size_t* numbersToWrite,
        size_t numNumbers) const {
      const size_t alignmentSizeT = alignof(size_t);
      const size_t sizeOfSizeT = sizeof(size_t);

      memcpy(buffer, &numChars, sizeOfSizeT);
      buffer += sizeOfSizeT;
      memcpy(buffer, &numNumbers, sizeOfSizeT);
      buffer += sizeOfSizeT;
      memcpy(buffer, charsToWrite, numChars);
      //for (size_t i = 0; i < numChars; ++i)
      //  std::cout << buffer[i];
      //std::cout << std::endl;
      buffer += numChars;
      const size_t bytesToWrite = sizeof(size_t) * numNumbers;
      memcpy(buffer, numbersToWrite, bytesToWrite);
      //for(size_t i = 0; i < bytesToWrite; ++i)
      //  std::cout << (int) buffer[i];
      //std::cout << std::endl;
      buffer += bytesToWrite;
      return buffer;
    }

    void read_(unsigned char* buffer,
        std::vector<unsigned char>& charsToRead,
        std::vector<size_t>& numbersToRead) const {

      size_t numCharsToRead = 0;
      size_t numNumbersToRead = 0;
      memcpy(&numCharsToRead, buffer, sizeof(size_t));
      buffer += sizeof(size_t);
      memcpy(&numNumbersToRead, buffer, sizeof(size_t));
      buffer += sizeof(size_t);
      //std::cout << "READ: numberChar " << numCharsToRead << " numberNumbers " << numNumbersToRead << std::endl;
      charsToRead.clear();
      charsToRead.reserve(numCharsToRead);
      numbersToRead.clear();
      numbersToRead.reserve(numNumbersToRead);

      for (size_t i = 0; i < numCharsToRead; ++i) 
        charsToRead.emplace_back(buffer[i]);
      for (size_t i = 0; i < numNumbersToRead; ++i) {
        size_t tmp;
        std::memcpy(&tmp, buffer + numCharsToRead + i * sizeof(size_t), sizeof(size_t));
        numbersToRead.emplace_back(tmp);
      }  
    }
    std::pair<std::vector<unsigned char>, std::vector<size_t>> read(unsigned char* buffer, size_t size) const {

      auto recvData = computeNumberOfRecvData(buffer, size);
      size_t numCharsToRead = recvData.first;
      size_t numNumbersToRead = recvData.second;
      dsss::mpi::environment env;
      std::cout << "env: " << env.rank() << " numChars: " << numCharsToRead << " numNumbers: " <<numNumbersToRead << std::endl; 
      std::vector<unsigned char> charsToRead;
      std::vector<size_t> numbersToRead;
      charsToRead.reserve(numCharsToRead);
      numbersToRead.reserve(numNumbersToRead);

      size_t curPos = 0;
      while(curPos < size) {
        memcpy(&numCharsToRead, buffer + curPos, sizeof(size_t));
        curPos += sizeof(size_t);
        memcpy(&numNumbersToRead, buffer + curPos, sizeof(size_t));
        curPos += sizeof(size_t);

        for(size_t i = 0; i < numCharsToRead; ++i) 
          charsToRead.emplace_back(*(buffer + curPos + i));
        curPos += numCharsToRead;
        for (size_t i = 0; i < numNumbersToRead; ++i) {
          size_t tmp;
          std::memcpy(&tmp, buffer + curPos + i * sizeof(size_t), sizeof(size_t));
          numbersToRead.emplace_back(tmp);
        }
        curPos += numNumbersToRead * sizeof(size_t);
      }
      return make_pair(std::move(charsToRead), std::move(numbersToRead));
    }
};


template <typename StringSet>
dss_schimek::StringLcpContainer<StringSet> alltoallv(
    dss_schimek::StringLcpPtr<StringSet>& strptr,
    const std::vector<size_t>& send_counts,
    environment env = environment()){

  using namespace dss_schimek;
  using String = typename StringSet::String;
  using CharIterator = typename StringSet::CharIterator;
  const StringSet& sendSet = strptr.active();
  if (sendSet.size() == 0)
    return StringLcpContainer<StringSet>();

  std::vector<unsigned char> receive_buffer_char;
  std::vector<size_t> receive_buffer_lcp;
  std::vector<size_t> send_counts_lcp(send_counts);
  std::vector<size_t> send_counts_char(send_counts.size(), 0);

  std::vector<unsigned char> send_buffer;

  for (size_t interval = 0, offset = 0; interval < send_counts.size(); ++interval) {
    for (size_t j = offset; j < send_counts[interval] + offset; ++j) {
      size_t string_length = sendSet.get_length(sendSet[sendSet.begin() + j]); 
      send_counts_char[interval] += string_length;
    }
    offset += send_counts[interval];
  }
  size_t totalNumSendChars = 
    std::accumulate(send_counts_char.begin(), send_counts_char.end(), 0);
  size_t totalNumStrings = 
    std::accumulate(send_counts_lcp.begin(), send_counts_lcp.end(), 0);

  size_t totalBytesToSend = totalNumSendChars + sizeof(size_t) * totalNumStrings;

  std::cout << "totalBytesToSend: " << totalBytesToSend << std::endl;
  std::cout << "totalNumSendChar: " << totalNumSendChars << std::endl;
  std::cout << "totalNumStrings: " << totalNumStrings << std::endl;
  

  unsigned char* buffer = new unsigned char[totalBytesToSend]; 

  
  for (size_t interval = 0, bufferOffset = 0, stringOffset = 0;
      interval < send_counts.size(); ++interval) {

    for(size_t j = 0; j < send_counts[interval]; ++j) {
      String str = sendSet[sendSet.begin() + stringOffset + j];
      CharIterator charIt = sendSet.get_chars(str, 0);
      size_t stringLength = sendSet.get_length(str) + 1;
      std::cout << charIt << std::endl;
      memcpy(buffer + bufferOffset, charIt, stringLength);
      std::cout << buffer + bufferOffset << std::endl;
      bufferOffset += stringLength;
    }
    size_t sizeOfLcpChunk = send_counts[interval] + sizeof(size_t);
    memcpy(buffer + bufferOffset, strptr.lcp_array() + stringOffset, sizeOfLcpChunk);
    bufferOffset += sizeOfLcpChunk; 
    stringOffset += send_counts[interval];
  }
  

  for (size_t i = 0; i < totalBytesToSend; ++i) {
    std::cout << (int) buffer[i] << ";";
  }
  std::cout << std::endl;
  return StringLcpContainer<StringSet>();
}

template <typename StringSet>
dss_schimek::StringLcpContainer<StringSet> alltoallv_2(
    dss_schimek::StringLcpContainer<StringSet>& container,
    const std::vector<size_t>& sendCountsString,
    environment env = environment()){

  using namespace dss_schimek;
  using String = typename StringSet::String;
  using CharIterator = typename StringSet::CharIterator;

  const ByteEncoder byteEncoder;
  const StringSet& sendSet = container.make_string_set();
  std::vector<unsigned char> contiguousStrings = 
    dss_schimek::getContiguousStrings(sendSet, container.char_size());

  //for (size_t i = 0, offset = 0; i < sendSet.size(); ++i) {
  //  size_t length = dss_schimek::string_length(contiguousStrings.data() + offset);
  //  std::cout << contiguousStrings.data() + offset << std::endl;
  //  offset += length + 1;
  //}

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
  
  //for (size_t i = 0; i < totalNumberSendBytes; ++i) {
  //  std::cout << (int) buffer[i] << ";";
  //}
  //std::cout << std::endl;

  std::vector<unsigned char> recv = alltoallv_small(buffer, sendCountsTotal);
   
  std::vector<unsigned char> rawStrings;
  std::vector<size_t> rawLcps;
  std::tie(rawStrings, rawLcps) = byteEncoder.read(recv.data(), recv.size());

  return StringLcpContainer<StringSet>(std::move(rawStrings), std::move(rawLcps));
}

template <typename StringSet>
dss_schimek::StringLcpContainer<StringSet> alltoallv_3(
    dss_schimek::StringLcpContainer<StringSet>& container,
    const std::vector<size_t>& sendCountsString,
    environment env = environment()){

  using namespace dss_schimek;
  using String = typename StringSet::String;
  using CharIterator = typename StringSet::CharIterator;

  const ByteEncoder byteEncoder;
  const StringSet& sendSet = container.make_string_set();

  //for (size_t i = 0, offset = 0; i < sendSet.size(); ++i) {
  //  size_t length = dss_schimek::string_length(contiguousStrings.data() + offset);
  //  std::cout << contiguousStrings.data() + offset << std::endl;
  //  offset += length + 1;
  //}

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

    std::cout << "sendSet size: " << sendSet.size() << " sendCountsLcp[interval]: " << sendCountsLcp[interval] << std::endl;
    auto begin = sendSet.begin() + stringsWritten;
    StringSet subSet = sendSet.sub(begin, begin + sendCountsLcp[interval]);
    //subSet.print();
    curPos = byteEncoder.write(curPos, 
        subSet,
        container.lcp_array() + stringsWritten);
    stringsWritten += sendCountsLcp[interval];
  }
  
  //std::cout << "rank: " << env.rank() << " bytes written" << std::endl;
  //for (size_t i = 0; i < totalNumberSendBytes; ++i) {
  //  std::cout << (int) buffer[i] << ";";
  //}
  //std::cout << std::endl;

  std::vector<unsigned char> recv = alltoallv_small(buffer, sendCountsTotal);
  // 
  std::vector<unsigned char> rawStrings;
  std::vector<size_t> rawLcps;
  std::tie(rawStrings, rawLcps) = byteEncoder.read(recv.data(), recv.size());

  return StringLcpContainer<StringSet>(std::move(rawStrings), std::move(rawLcps));
}
} // namespace dsss::mpi

/******************************************************************************/
