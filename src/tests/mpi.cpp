#include <iostream>
#include <vector>
#include <random>
#include <numeric>

#include "mpi/environment.hpp"
#include "mpi/allgather.hpp"
#include "mpi/alltoall.hpp"
#include "util/timer.hpp"
#include "mpi/synchron.hpp"

#include <tlx/cmdline_parser.hpp>

class InitializedUniformIntDistribution {
  public:
    InitializedUniformIntDistribution(const size_t min, const size_t max) 
    : min(min), max(max), generator(initialSeed()), dist(min, max) {};

  size_t operator() () {
    return dist(generator);
  }

  private: 
  size_t min;
  size_t max;
  std::mt19937 generator;
  std::uniform_int_distribution<unsigned char> dist;

  size_t initialSeed() {
    static std::random_device rd;
    static size_t seed(rd());
    return seed;
  }


};


std::vector<unsigned char> allgatherv(const size_t sizeInBytes, dss_schimek::Timer& timer) {
  dsss::mpi::environment env;
  InitializedUniformIntDistribution dist(65, 90);
  std::vector<unsigned char> sendData;
  sendData.reserve(sizeInBytes);

    for (size_t i = 0; i < sizeInBytes - 1; ++i) 
      sendData.push_back(dist());
    sendData.push_back(0);
  timer.start("allgatherv");
  std::vector<unsigned char> recvData = dsss::mpi::allgatherv(sendData);
  timer.end("allgatherv");
  return recvData; 
}

std::vector<unsigned char> alltoall(const size_t sizeInBytesPerPE, dss_schimek::Timer& timer) {
  dsss::mpi::environment env;
  InitializedUniformIntDistribution dist(65, 90);
  std::vector<size_t> sendCounts(env.size(), sizeInBytesPerPE);
  std::vector<unsigned char> sendData;
  sendData.reserve(sizeInBytesPerPE * env.size());

  for (size_t j = 0; j < env.size(); ++j) {
    for (size_t i = 0; i < sizeInBytesPerPE - 1; ++i) 
      sendData.push_back(dist());
    sendData.push_back(0);
  }
  timer.start("alltoall");
  std::vector<unsigned char> recvData = dsss::mpi::AllToAllvSmall::alltoallv(sendData.data(), sendCounts);
  timer.end("alltoall");
  return recvData; 
}

void doNotOptimizeAway(const std::vector<unsigned char> data) {
  volatile size_t sum = 0;
  sum = std::accumulate(data.begin(), data.end(), 0);
}

void runIteration(const size_t sizeInBytesAllgather, const size_t sizeInBytesAllToAll, const size_t curIteration) {
  dsss::mpi::environment env;
  std::string prefix = std::string("RESULT") +
    " numberProcessors=" + std::to_string(env.size()) +
    " iteration=" + std::to_string(curIteration);
  " sizeInBytesAllgather=" + std::to_string(sizeInBytesAllgather) +
    " sizeInBytesAllToAll=" + std::to_string(sizeInBytesAllToAll);

  dss_schimek::Timer timer(prefix);
  doNotOptimizeAway(allgatherv(sizeInBytesAllgather, timer));
  doNotOptimizeAway(alltoall(sizeInBytesAllToAll, timer));

  std::stringstream buffer;
  timer.writeToStream(buffer);
  if (env.rank() == 0) {
    std::cout << buffer.str() << std::endl;
  }
}

void run(const size_t sizeInBytesAllgather, const size_t sizeInBytesAllToAll, const size_t numberOfIterations) {
  for(size_t i = 0; i < numberOfIterations; ++i)
    runIteration(sizeInBytesAllgather, sizeInBytesAllToAll, i);
}


int main (int argc, char *argv[]) {

  dsss::mpi::environment env;
  tlx::CmdlineParser cp;

  unsigned int iterations = 1;
  unsigned int sizeInBytesAllgather = 10000;
  unsigned int sizeInBytesAllToAll = 10000;

  cp.add_unsigned('_', "sizeAllGather", sizeInBytesAllgather, " number of bytes to send");
  cp.add_unsigned('_', "sizeAllToAll", sizeInBytesAllToAll, " number of bytes to send");
  cp.add_unsigned('i', "numberOfIterations", iterations, "");
  

  if (!cp.process(argc, argv)) {
    return -1;
  }

  if (env.rank() == 0) {
    std::cout << "startup info:" << std::endl;
    cp.print_result();
  }
 
  run(sizeInBytesAllgather, sizeInBytesAllToAll, iterations);
 
  env.finalize();
}
