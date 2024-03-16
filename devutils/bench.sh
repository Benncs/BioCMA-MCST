#!/bin/bash

# Check if the correct number of arguments is provided
if [ "$#" -lt 1 ]; then
    echo "Usage: $0 <command> [args...]"
    exit 1
fi

# Command to be benchmarked
command_to_run="$@"

# Run the command and measure the time
echo "Benchmarking: $command_to_run"
/usr/bin/time -f "\nTime elapsed: %e seconds\nCPU usage: %P" $command_to_run

exit $?