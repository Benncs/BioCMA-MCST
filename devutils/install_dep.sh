#!/bin/bash

LLVM_VERSION=18
INSTALL_HDF5=false

# Function to print error messages
error() {
  echo "$1" >&2
}

# Function to check if a command exists
command_exists() {
  command -v "$1" &> /dev/null
}


# Check if apt-get is installed
if ! command_exists apt-get; then
  error "apt-get could not be found. This script is for Ubuntu systems."
  exit 1
fi

# Check if sudo is installed
if ! command_exists sudo; then
  error "sudo could not be found. Please install sudo or run this script as root."
  exit 1
fi

# Update package list
echo "Updating package list..."
if ! sudo apt-get update; then
  error "Failed to update package list."
  exit 1
fi

# Install necessary packages
echo "Installing necessary packages..."
if ! sudo apt-get install -y python3-dev python3 python3-pip wget software-properties-common gnupg libomp-dev libopenmpi-dev libtbb-dev libeigen3-dev pkg-config ninja-build; then
  error "Failed to install necessary packages."
  exit 1
fi

# Check for --hdf5 flag
for arg in "$@"; do
  if [ "$arg" == "--hdf5" ]; then
    INSTALL_HDF5=true
    break
  fi
done

# Install optional HDF5 if flag is set
if $INSTALL_HDF5; then
  echo "Installing HDF5..."
  if ! sudo apt-get install -y libhdf5-dev; then
    error "Failed to install HDF5."
    exit 1
  fi
fi

# Install meson via pip
echo "Installing meson..."
if ! pip3 install meson --break-system-packages; then
  error "Failed to install meson."
  exit 1
fi

echo "Installing omp dev..."
  if ! sudo apt-get install -y libomp-18-dev; then
    error "Failed to install HDF5."
    exit 1
  fi

# Run clang configuration script
echo "Running clang configuration script..."
if ! sh ./devutils/docker/clang_config.sh $LLVM_VERSION; then
  error "Failed to run clang configuration script."
  exit 1
fi

# Install specific version of clang
echo "Installing clang-$LLVM_VERSION..."
if ! sudo apt-get install -y clang-$LLVM_VERSION; then
  error "Failed to install clang-$LLVM_VERSION."
  exit 1
fi

echo "All dependencies have been installed successfully."
