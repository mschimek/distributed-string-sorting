#include <iostream>
#include <vector>
#include <random>

#include "mpi/environment.hpp"
#include "mpi/allgather.hpp"
#include "util/timer.hpp"
#include "mpi/synchron.hpp"

#include <tlx/cmdline_parser.hpp>

void allgather(size_t size) {
  std::random_device rd;
  std::mt19937 generator(rd);
  std::uniform_int_distribution dist(65, 90);
  std::vector<unsigned char> sendData;
  sendData.reserve(size);
  
  for (size_t i = 0; i < size; ++i) 
    sendData.push_back(dist(generator));

  dsss::mpi::data_type_mapper<unsigned char> dtm;
  std::vector<unsigned char> recvData(env.size());
  MPI_Allgather(
      sendData.data(),
      size,
      dtm.get_mpi_type(),
      recvData.data(),
      1,
      dtm.get_mpi_type(),
      env.communicator());


  dsss::mpi::allgather(sendData.data(), ) 
}

int main (int argc, char *argv[]) {

  dsss::mpi::environment env;
  tlx::CmdlineParser cp;

  unsigned int iterations = 1;
  unsigned int sizeInBytes = 10000;

  cp.add_unsigned('s', "size", sizeInBytes, " number of bytes to send");
  cp.add_unsigned('i', "numberOfIterations", iterations, "");

  for (size_t i = 0; i < iterations; ++i) {

  }
  
  env.finalize();
}
