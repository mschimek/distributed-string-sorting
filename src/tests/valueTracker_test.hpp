#pragma once

#include <tlx/die/core.hpp>
#include "util/valueTracker.hpp"
#include "mpi/environment.hpp"
#include "mpi/allgather.hpp"
#include "mpi/alltoall.hpp"

void valueTracker_test_1() {
  dsss::mpi::environment env;

  ValueTracker tracker;
  tracker.setPhase("init");
  tracker.add(42);
  tracker.add(42);
  tracker.setPhase("work");
  tracker.add(42, 10);
  tracker.add(42, 10, "addingRandomNumbers");

  std::vector<ValueTracker::Result> results;
  tracker.collect(std::back_inserter(results));

  for (size_t i = 0; i < env.size(); ++i) {
    auto it = std::find_if(results.begin(), results.end(), [&] (const ValueTracker::Result& result) {
        return result.phase == "init" && result.counterPerPhase == 0 && result.description == "none" && result.value == 42;
        });
    tlx_die_if((it == results.end()));
    results.erase(it);
  }
  
  for (size_t i = 0; i < env.size(); ++i) {
    auto it = std::find_if(results.begin(), results.end(), [&] (const ValueTracker::Result& result) {
        return result.phase == "init" && result.counterPerPhase == 1 && result.description == "none" && result.value == 42;
        });
    tlx_die_if((it == results.end()));
    results.erase(it);
  }
  
  for (size_t i = 0; i < env.size(); ++i) {
    auto it = std::find_if(results.begin(), results.end(), [&] (const ValueTracker::Result& result) {
        return result.phase == "work" && result.counterPerPhase == 0 && result.description == "none" && result.value == 42 && result.round == 10;
        });
    tlx_die_if((it == results.end()));
    results.erase(it);
  }
  
  for (size_t i = 0; i < env.size(); ++i) {
    auto it = std::find_if(results.begin(), results.end(), [&] (const ValueTracker::Result& result) {
        return result.phase == "work" && result.counterPerPhase == 1 && result.description == "addingRandomNumbers" && result.value == 42 && result.round == 10;
        });
    tlx_die_if((it == results.end()));
    results.erase(it);
  }

  tlx_die_unless((results.empty()));
}

void valueTracker_test_2() {
  dsss::mpi::environment env;

  ValueTracker& tracker = ValueTracker::valueTracker();
  tracker.setPhase("init");
  std::vector<size_t> testData(env.size(), 42);
  dsss::mpi::alltoall(testData);

  std::vector<ValueTracker::Result> results;
  tracker.collect(std::back_inserter(results));

  for (size_t i = 0; i < env.size(); ++i) {
    auto it = std::find_if(results.begin(), results.end(), [&] (const ValueTracker::Result& result) {
        return result.phase == "init" && result.counterPerPhase == 0 && result.description == "none" && result.value == sizeof(size_t) * env.size();
        });
    tlx_die_if((it == results.end()));
    results.erase(it);
  }
  tlx_die_unless((results.empty()));
}
