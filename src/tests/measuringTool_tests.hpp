#pragma once

#include <tlx/die/core.hpp>
#include "util/measuringTool.hpp"
#include "mpi/environment.hpp"

void measuringTool_basicTest1() {
  dsss::mpi::environment env;
  using namespace dss_schimek::measurement;

  MeasuringTool& measuringTool = MeasuringTool::measuringTool();

  measuringTool.setPhase("init");
  measuringTool.add(42);
  measuringTool.add(43);
  measuringTool.setPhase("work");
  measuringTool.start("myOp");
  measuringTool.stop("myOp");

  auto data = measuringTool.collect();

  std::vector<std::string> types = {"avgTime", "minTime", "maxTime", "avgLoss", "minLoss", "maxLoss"};

  tlx_die_unless(data.size() == (env.size() * 2 + 6));

  for (size_t i = 0; i < env.size(); ++i) {
    auto it = std::find_if(data.begin(), data.end(), [&] (const auto& elem) {
        return elem.phase == "init" && elem.counterPerPhase == 0 && elem.description == "unkown" && elem.value == 42;
        });
    tlx_die_if((it == data.end()));
    data.erase(it);
  }
  
  for (size_t i = 0; i < env.size(); ++i) {
    auto it = std::find_if(data.begin(), data.end(), [&] (const auto& elem) {
        return elem.phase == "init" && elem.counterPerPhase == 1 && elem.description == "unkown" && elem.value == 43;
        });
    tlx_die_if((it == data.end()));
    data.erase(it);
  }
  
  for (const std::string& type : types) {
    auto it = std::find_if(data.begin(), data.end(), [&] (const auto& elem) {
        return elem.phase == "work" && elem.counterPerPhase == 0 && elem.description == "myOp" && elem.type == type;
        });
    tlx_die_if((it == data.end()));
    data.erase(it);
  }
  //
  //for (size_t i = 0; i < env.size(); ++i) {
  //  auto it = std::find_if(results.begin(), results.end(), [&] (const ValueTracker::Result& result) {
  //      return result.phase == "work" && result.counterPerPhase == 0 && result.description == "none" && result.value == 42 && result.round == 10;
  //      });
  //  tlx_die_if((it == results.end()));
  //  results.erase(it);
  //}
  //
  //for (size_t i = 0; i < env.size(); ++i) {
  //  auto it = std::find_if(results.begin(), results.end(), [&] (const ValueTracker::Result& result) {
  //      return result.phase == "work" && result.counterPerPhase == 1 && result.description == "addingRandomNumbers" && result.value == 42 && result.round == 10;
  //      });
  //  tlx_die_if((it == results.end()));
  //  results.erase(it);
  //}

  tlx_die_unless((data.empty()));
}

//void valueTracker_test_2() {
//  dsss::mpi::environment env;
//
//  ValueTracker& tracker = ValueTracker::valueTracker();
//  tracker.setPhase("init");
//  std::vector<size_t> testData(env.size(), 42);
//  dsss::mpi::alltoall(testData);
//
//  std::vector<ValueTracker::Result> results;
//  tracker.collect(std::back_inserter(results));
//
//  for (size_t i = 0; i < env.size(); ++i) {
//    auto it = std::find_if(results.begin(), results.end(), [&] (const ValueTracker::Result& result) {
//        return result.phase == "init" && result.counterPerPhase == 0 && result.description == "none" && result.value == sizeof(size_t) * env.size();
//        });
//    tlx_die_if((it == results.end()));
//    results.erase(it);
//  }
//  tlx_die_unless((results.empty()));
//}
