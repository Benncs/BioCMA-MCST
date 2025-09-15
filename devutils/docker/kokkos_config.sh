#!/bin/bash

check_sudo() {
  command -v sudo &>/dev/null
}
cleanup() {
  cd "$current_pwd" || exit 1
}
get_kokkos_source() {
  echo "Getting Kokkos.."
  echo $tar_url
  wget -q "$1" -O "$2"

  if [ $? -ne 0 ]; then
    echo "Error: Failed to download Kokkos source."
    exit 1
  fi

  tar -xvf "$2" >/dev/null

  if [ $? -ne 0 ]; then
    echo "Error: Failed to extract Kokkos archive."
    exit 1
  fi
}
set -e # Exit immediately if a command exits with a non-zero status

current_pwd=$(pwd)
back_end_omp=0
back_end_cuda=0
clang_version=-1
kokkos_version="4.7.00"
while [[ "$#" -gt 0 ]]; do
  case $1 in
  --omp) back_end_omp=1 ;;
  --cuda) back_end_cuda=1 ;;
  --clang)
    if [[ -n "$2" && "$2" != --* ]]; then
      clang_version="$2"
      shift
    else
      echo "Error: --clang option requires a version number."
      exit 1
    fi
    ;;
  *)
    if [[ "$1" != --* ]]; then
      kokkos_version="$1"
    else
      echo "Unknown option: $1"
      exit 1
    fi
    ;;
  esac
  shift
done

folder_name="kokkos-$kokkos_version"
tar_name="${folder_name}.tar.gz"
tar_url="https://github.com/kokkos/kokkos/releases/download/$kokkos_version/$tar_name"

trap cleanup EXIT

cd /tmp || {
  echo "Error: Failed to change directory to /tmp"
  exit 1
}

get_kokkos_source "$tar_url" "$tar_name"

cd "$folder_name" || {
  echo "Error: Failed to change directory to $folder_name"
  exit 1
}

mkdir -p kokkos_build
cd kokkos_build || {
  echo "Error: Failed to change directory to kokkos_build"
  exit 1
}

flag_cmake="-DCMAKE_POSITION_INDEPENDENT_CODE=ON  -DCMAKE_CXX_STANDARD=20 -B . -S .. -DCMAKE_BUILD_TYPE=Release"

if [[ "$clang_version" != "-1" ]]; then
  flag_cmake="${flag_cmake} -DCMAKE_CXX_COMPILER=clang++-$clang_version"
fi

#flag_cmake="${flag_cmake} -DKokkos_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE=ON"
# flag_cmake="${flag_cmake} -DCMAKE_CXX_COMPILER=clang++-$clang_version"

if [[ "$back_end_omp" == "1" ]]; then
  flag_cmake="${flag_cmake} -DKokkos_ENABLE_OPENMP=ON"
fi

if [[ "$back_end_cuda" == "1" ]]; then
  # flag_cmake="${flag_cmake} -DCUDA_ROOT=/usr/local/cuda-12.6/"
  flag_cmake="${flag_cmake} -DKokkos_ARCH_TURING75=ON"
  flag_cmake="${flag_cmake} -DKokkos_ENABLE_CUDA=ON"
  flag_cmake="${flag_cmake} -DKokkos_ENABLE_CUDA_CONSTEXPR=ON"
fi

flag_cmake="${flag_cmake} -DKokkos_ENABLE_HWLOC=ON"
cmake $flag_cmake
cmake --build .

if check_sudo; then
  sudo cmake --install .
else
  cmake --install .
fi

cd /tmp
rm -rf /tmp/$tar_name /tmp/$folder_name
