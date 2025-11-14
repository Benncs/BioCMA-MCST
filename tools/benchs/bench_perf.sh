export OMP_PROC_BIND=spread
export OMP_PLACES=threads

#CLI=$(echo $(./tools/runner.py debug -n 1 --dry-run) | cut -d' ' -f2-)
CLI=$(echo $(./tools/runner.py bench -n 1 --dry-run))

FLAG="-g -q --call-graph dwarf "

EVENT="-e cache-misses,cycles,instructions"

OMP_NUM_THREADS=6 VLAMBDA=693147.18056 perf record $FLAG $EVENT $CLI

chown $(whoami) perf.data
mv perf.data /tmp/perf_$(date +"%Y%m%d_%H%M%S").data
