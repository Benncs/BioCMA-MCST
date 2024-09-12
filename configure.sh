#!/bin/bash
# -----------------------------------------------------------------------------
# Author: CASALE Benjamin
# Date: 09/09/2024
# Description: Script to configure build directory
# -----------------------------------------------------------------------------

if [ $# -eq 0 ]; then
    echo "No arguments provided. Exiting."
    exit 1
elif [ $# -eq 1 ]; then
    clang_compiler=""
    gcc_compiler=$1
elif [ $# -eq 2 ]; then
    clang_compiler=$2
    gcc_compiler=$1
else
    echo "Too many arguments provided. Exiting."
    exit 1
fi

compilers=($clang_compiler $gcc_compiler)

compilers=($(echo "${compilers[@]}" | tr ' ' '\n' | grep -v '^$' | tr '\n' ' '))


compilers_name=("gcc" "clang")
types=("debugoptimized" "release")
types_name=("debug" "release")

function exec() {
    local current_compiler=$1
    local current_name=$2
    local current_type=$3
    CXX=$current_compiler meson setup __builddir/$current_name --buildtype=$current_type >>/dev/null
}

index_compiler=0
for compiler in "${compilers[@]}"; do
    folder=${compilers_name[$index_compiler]}
    for i in "${!types[@]}"; do
        current_type=${types[$i]}
        current_name=${types_name[$i]}_$folder
        exec $compiler $current_name $current_type
        
        if [ $? -eq 0 ]; then
            echo "Success: $compiler for $current_name with type $current_type"
      
        else
            echo "Failure: $compiler for $current_name with type $current_type"
        fi  
    done
    ((index_compiler++))  
done