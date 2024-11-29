export OMP_PROC_BIND=spread
export OMP_PLACES=threads

CLI=$(echo $(./tools/runner.py bench -n 4 --dry-run) | cut -d' ' -f2-)

FLAG="-g -q --call-graph dwarf "

EVENT="-e cache-misses,cycles,instructions"

perf record $FLAG $EVENT  $CLI
chown benjamin perf.data
