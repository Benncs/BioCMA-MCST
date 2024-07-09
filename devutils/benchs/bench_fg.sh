#!/bin/bash

process_name="BioCMA-MCST_cli"

flamegraph_path="/mnt/c/Users/casale/Documents/code/FlameGraph/"

PID=7124   
#$(pgrep $process_name)

# Check if PID is empty
if [ -z "$PID" ]; then
    echo "Error: Process not found or not running."
    exit 1
fi

# Check if PID and duration are provided as arguments
if [ "$#" -ne 1 ]; then
    echo "Usage: $0  <duration>"
    exit 1
fi

DURATION="$1"

# # Record performance data using perf
echo "Recording performance data for PID $PID for $DURATION seconds..."
sudo perf record -F 99 -g -p "$PID" -a -- sleep "$DURATION"

# Convert perf data to folded stack format
echo "Converting perf data to folded stack format..."
sudo perf script | $flamegraph_path"stackcollapse-perf.pl" --all > out.folded

# Generate flame graph
echo "Generating flame graph..."
$flamegraph_path"flamegraph.pl" out.folded > flamegraph.svg

# # Open flame graph in web browser
# echo "Opening flame graph in web browser..."
# xdg-open flamegraph.svg &> /dev/null

echo "Done."
