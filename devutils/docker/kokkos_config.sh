#!/bin/bash

set -e  # Exit immediately if a command exits with a non-zero status

current_pwd=$(pwd)
kokkos_version=4.4.00
folder_name="kokkos-4.4.00"
tar_name="${folder_name}.tar.gz"
tar_url="https://github.com/kokkos/kokkos/releases/download/$kokkos_version/kokkos-$kokkos_version.tar.gz"

back_end_omp="1"
back_end_cuda="0"


cleanup() {
    cd "$current_pwd" || exit 1
}

trap cleanup EXIT

get_kokkos_source() {
    echo "Getting Kokkos.."
    wget -q "$1" -O "$2"

    if [ $? -ne 0 ]; then
        echo "Error: Failed to download Kokkos source."
        exit 1
    fi

    tar -xvf "$2" > /dev/null

    if [ $? -ne 0 ]; then
        echo "Error: Failed to extract Kokkos archive."
        exit 1
    fi
}

cd /tmp || { echo "Error: Failed to change directory to /tmp"; exit 1; }

get_kokkos_source "$tar_url" "$tar_name"

cd "$folder_name" || { echo "Error: Failed to change directory to $folder_name"; exit 1; }

mkdir -p kokkos_build
cd kokkos_build || { echo "Error: Failed to change directory to kokkos_build"; exit 1; }

flag_cmake="-DCMAKE_POSITION_INDEPENDENT_CODE=ON -DCMAKE_CXX_STANDARD=20 -B . -S .. -DCMAKE_BUILD_TYPE=Release"

if [[ "$back_end_omp" == "1" ]]; then
    flag_cmake="${flag_cmake} -DKokkos_ENABLE_OPENMP=ON"
fi

if [[ "$back_end_cuda" == "1" ]]; then
    flag_cmake="${flag_cmake} -DKokkos_ENABLE_CUDA=ON"
fi

cmake $flag_cmake
cmake --build .
cmake --install .

cd /tmp
rm -rf /tmp/*
