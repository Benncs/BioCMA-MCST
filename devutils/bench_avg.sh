#!/bin/bash

# Check if the correct number of arguments is provided
if [ "$#" -lt 1 ]; then
    echo "Usage: $0 <command> [args...]"
    exit 1
fi

num_runs=5

# Command to be benchmarked
command_to_run="$@"

# Initialize total time
total_time=0

# Run the command multiple times and measure the time
echo "Benchmarking: $command_to_run"
for ((i=1; i<=$num_runs; i++))
do
    echo "Run $i"
    # Measure time for each run
    time_output=$( { /usr/bin/time -f "%e" $command_to_run 2>&1; } 2>&1 )
    # Extract the elapsed time
    elapsed_time=$(echo "$time_output" | tail -n 1)
    # Add to total time
    total_time=$(echo "$total_time + $elapsed_time" | bc)
done

# Calculate average time
average_time=$(echo "scale=2; $total_time / $num_runs" | bc)

# Print average time
echo -e "\nAverage Time elapsed over $num_runs runs: $average_time seconds"

exit $?