#!/usr/bin/bash 

PID=$(pgrep biocma)
echo "sending SIGUSR1 to $PID"
kill -SIGUSR1 $PID 
