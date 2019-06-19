#!/bin/bash -x

rank=${OMPI_COMM_WORLD_RANK}
if [ -z "$rank"  ]; then
    rank=${PMI_RANK}
    if [ -z "$rank" ]; then
        rank=${SLURM_NODEID}
    fi
fi

exec perf record -o perf-${rank}.data $@
