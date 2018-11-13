/*******************************************************************************
 * mpi/shift.hpp
 *
 * Copyright (C) 2018 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include <algorithm>
#include <cstdint>
#include <mpi.h>

#include "mpi/environment.hpp"
#include "mpi/type_mapper.hpp"
#include "mpi/allreduce.hpp"
#include "util/string.hpp"
#include "strings/stringtools.hpp"

namespace dsss::mpi {

template <typename DataType>
static inline DataType shift_left(DataType& send_data,
  environment env = environment()) {

  if (env.size() == 1) {
    return send_data;
  }

  std::int32_t destination = env.size() - 1;
  if (env.rank() > 0) {
    destination = env.rank() - 1;
  }
  std::int32_t source = 0;
  if (env.rank() + 1 < env.size()) {
    source = env.rank() + 1;
  }

  data_type_mapper<DataType> dtm;
  DataType receive_data;
  MPI_Sendrecv(&send_data,
               1,
               dtm.get_mpi_type(),
               destination,
               0, // chose arbitrary tag
               &receive_data,
               1,
               dtm.get_mpi_type(),
               source,
               MPI_ANY_TAG,
               env.communicator(),
               MPI_STATUS_IGNORE);
  return receive_data;
}

template <typename DataType>
static inline DataType shift_right(DataType& send_data,
  environment env = environment()) {

  if (env.size() == 1) {
    return send_data;
  }

  std::int32_t destination = 0;
  if (env.rank() + 1 < env.size()) {
    destination = env.rank() + 1;
  }
  std::int32_t source = env.size() - 1;
  if (env.rank() > 0) {
    source = env.rank() - 1;
  }

  data_type_mapper<DataType> dtm;
  DataType receive_data;
  MPI_Sendrecv(&send_data,
               1,
               dtm.get_mpi_type(),
               destination,
               0, // chose arbitrary tag
               &receive_data,
               1,
               dtm.get_mpi_type(),
               source,
               MPI_ANY_TAG,
               env.communicator(),
               MPI_STATUS_IGNORE);
  return receive_data;
}

template <typename DataType>
static inline std::vector<DataType> shift_left(DataType* send_data,
  std::size_t count, environment env = environment()) {

  std::int32_t destination = env.size() - 1;
  if (env.rank() > 0) {
    destination = env.rank() - 1;
  }
  std::int32_t source = 0;
  if (env.rank() + 1 < env.size()) {
    source = env.rank() + 1;
  }

  std::size_t receive_count = shift_left(count, env);

  std::vector<DataType> receive_data(receive_count);
  data_type_mapper<DataType> dtm;
  MPI_Sendrecv(send_data,
               count,
               dtm.get_mpi_type(),
               destination,
               0, // chose arbitrary tag
               receive_data.data(),
               receive_count,
               dtm.get_mpi_type(),
               source,
               MPI_ANY_TAG,
               env.communicator(),
               MPI_STATUS_IGNORE);
  return receive_data;
}

template <typename DataType>
static inline std::vector<DataType> shift_right(DataType* send_data,
  std::size_t count, environment env = environment()) {

  std::int32_t destination = 0;
  if (env.rank() + 1 < env.size()) {
    destination = env.rank() + 1;
  }
  std::int32_t source = env.size() - 1;
  if (env.rank() > 0) {
    source = env.rank() - 1;
  }

  std::size_t receive_count = shift_right(count, env);

  std::vector<DataType> receive_data(count);
  data_type_mapper<DataType> dtm;
  MPI_Sendrecv(send_data,
               count,
               dtm.get_mpi_type(),
               destination,
               0, // chose arbitrary tag
               receive_data.data(),
               receive_count,
               dtm.get_mpi_type(),
               source,
               MPI_ANY_TAG,
               env.communicator(),
               MPI_STATUS_IGNORE);
  return receive_data;
}

template<typename DataType>
struct ShiftReturn {
  DataType recv_data;
  bool all_PE_empty;
};

template <typename DataType>
static inline ShiftReturn<DataType> get_shift_left_receiv_counts(const DataType& send_data, 
    environment env = environment(), 
    bool is_empty = false) {

  if (env.size() == 1)
    return {send_data, is_empty};

  const bool is_overall_empty = dsss::mpi::allreduce_and(is_empty, env);
  if (is_overall_empty)
    return {send_data, true};

  std::cout << "scalar is_overall_empty " << is_overall_empty << std::endl;
  data_type_mapper<DataType> dtm;
  
  std::int32_t destination = env.size() - 1;
  if (env.rank() > 0) {
    destination = env.rank() - 1;
  }
  std::int32_t source = 0;
  if (env.rank() + 1 < env.size()) {
    source = env.rank() + 1;
  }

  std::cout << "scalar rank: " << env.rank() << " is_empty: " << is_empty << std::endl;
  DataType recv_data;
  if (is_empty) {
    MPI_Recv(
        &recv_data,
        1,
        dtm.get_mpi_type(),
        source,
        MPI_ANY_TAG,
        env.communicator(),
        MPI_STATUS_IGNORE);
    std::cout << "scalar rank: " << env.rank() << " recv_data  " << recv_data << std::endl;
    MPI_Send(
        &recv_data,
        1,
        dtm.get_mpi_type(),
        destination,
        42,
        env.communicator()); 
  } else {
    MPI_Request send_request, recv_request;
    MPI_Irecv(
        &recv_data,
        1,
        dtm.get_mpi_type(),
        source,
        MPI_ANY_TAG,
        env.communicator(),
        &recv_request);
    MPI_Isend(
        &send_data,
        1,
        dtm.get_mpi_type(),
        destination,
        42,
        env.communicator(),
        &send_request);
    MPI_Wait(&recv_request, MPI_STATUS_IGNORE);
    MPI_Wait(&send_request, MPI_STATUS_IGNORE);
  }
  return {recv_data, false};
}

static inline std::vector<unsigned char> shift_string_left_(
  unsigned char* send_data, environment env = environment()) {

  std::size_t string_length = dss_schimek::string_length(send_data);
  std::size_t send_length = string_length + 1;
  std::vector<char_type> real_send_data;
  std::copy_n(send_data, send_length, std::back_inserter(real_send_data));

  if (env.size() == 1) {
    return real_send_data;
  }

  ShiftReturn result = get_shift_left_receiv_counts(string_length + 1, env, string_length == 0);
  if (result.all_PE_empty)
    return real_send_data;

  size_t recv_length = result.recv_data;

  std::int32_t destination = env.size() - 1;
  if (env.rank() > 0) {
    destination = env.rank() - 1;
  }
  std::int32_t source = 0;
  if (env.rank() + 1 < env.size()) {
    source = env.rank() + 1;
  }

  std::vector<unsigned char> receive_data(recv_length); 
    data_type_mapper<unsigned char> dtm;
    MPI_Request send_request;
    MPI_Request recv_request;

  if (string_length > 0)
  {
      MPI_Isend(
      send_data,
      send_length,
      dtm.get_mpi_type(),
      destination,
      42, 
      env.communicator(),
      &send_request);
    MPI_Irecv(
      receive_data.data(),
      recv_length,
      dtm.get_mpi_type(),
      source,
      MPI_ANY_TAG,
      env.communicator(),
      &recv_request);
    MPI_Wait(&send_request, MPI_STATUS_IGNORE);
    MPI_Wait(&recv_request, MPI_STATUS_IGNORE);
  } else {
    MPI_Recv(
        receive_data.data(),
        recv_length,
        dtm.get_mpi_type(),
        source,
        MPI_ANY_TAG,
        env.communicator(),
        MPI_STATUS_IGNORE);
    MPI_Send(
        receive_data.data(),
        recv_length,
        dtm.get_mpi_type(),
        destination,
        42,
        env.communicator());
    return std::vector<unsigned char>(1, 0);
  }
  return receive_data;
  }

static inline std::vector<dsss::char_type> shift_string_left(
  dsss::string send_data, environment env = environment()) {

  std::size_t send_length = dsss::string_length(send_data) + 1;
  std::vector<char_type> real_send_data;
  std::copy_n(send_data, send_length, std::back_inserter(real_send_data));

  if (env.size() == 1) {
    return real_send_data;
  }

  std::size_t receive_length = shift_left(send_length, env);

  std::int32_t destination = env.size() - 1;
  if (env.rank() > 0) {
    destination = env.rank() - 1;
  }
  std::int32_t source = 0;
  if (env.rank() + 1 < env.size()) {
    source = env.rank() + 1;
  }

  std::vector<dsss::char_type> receive_data(receive_length);
  data_type_mapper<dsss::char_type> dtm;
  MPI_Sendrecv(real_send_data.data(),
               send_length,
               dtm.get_mpi_type(),
               destination,
               0, // chose arbitrary tag
               receive_data.data(),
               receive_length,
               dtm.get_mpi_type(),
               source,
               MPI_ANY_TAG,
               env.communicator(),
               MPI_STATUS_IGNORE);
  return receive_data;
}

static inline std::vector<dsss::char_type> shift_string_right(
  dsss::string send_data, environment env = environment()) {

  std::size_t send_length = dsss::string_length(send_data) + 1;
  std::vector<char_type> real_send_data;
  std::copy_n(send_data, send_length, std::back_inserter(real_send_data));

  if (env.size() == 1) {
    return real_send_data;
  }

  std::size_t receive_length = shift_right(send_length, env);

  std::int32_t destination = 0;
  if (env.rank() + 1 < env.size()) {
    destination = env.rank() + 1;
  }
  std::int32_t source = env.size() - 1;
  if (env.rank() > 0) {
    source = env.rank() - 1;
  }

  std::vector<dsss::char_type> receive_data(receive_length);
  data_type_mapper<dsss::char_type> dtm;
  MPI_Sendrecv(real_send_data.data(),
               send_length,
               dtm.get_mpi_type(),
               destination,
               0, // chose arbitrary tag
               receive_data.data(),
               receive_length,
               dtm.get_mpi_type(),
               source,
               MPI_ANY_TAG,
               env.communicator(),
               MPI_STATUS_IGNORE);
  return receive_data;
}

} // namespace dsss::mpi

/******************************************************************************/
