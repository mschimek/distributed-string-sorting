#pragma once

#include "mpi/environment.hpp"

namespace dss_schimek::mpi {
  template <typename Functor>
    void execute_in_order(const Functor& functor, 
        dss_schimek::mpi::environment env = dss_schimek::mpi::environment()){
      for (size_t i = 0; i < env.size(); ++i) {
        env.barrier();
        if (env.rank() == i)
          functor();
      volatile size_t j = 0;
      for(j = 0; j < 100000; ++j) j++;
        env.barrier();
      }

      env.barrier();
    }
}
