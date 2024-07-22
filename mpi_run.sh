#!/bin/bash
n_nodes=12

mpi_options="--use-hwthread-cpus --allow-run-as-root --bind-to hwthread -np $n_nodes"
bench_path=./devutils/benchs/bench.sh
run_path=./run.sh
OMP_NUM_THREADS=1 mpiexec  $mpi_options $bench_path $run_path