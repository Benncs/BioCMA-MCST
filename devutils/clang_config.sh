#!/bin/bash
# -----------------------------------------------------------------------------
# Author: CASALE Benjamin
# Date: 16/04/2024
# Description: Script to get the right version of clang-tidy and clang-format 
#              for CI and container execution. Not recommended to run this 
#              on a local machine.
# -----------------------------------------------------------------------------


# Define LLVM version
LLVM_VERSION=17

# Download LLVM installation script
wget https://apt.llvm.org/llvm.sh && \

# Make the script executable
chmod +x llvm.sh && \

# Run the LLVM installation script with the specified version
./llvm.sh $LLVM_VERSION && \

# Install clang-format and clang-tidy
apt-get install -y clang-format-$LLVM_VERSION clang-tidy-$LLVM_VERSION && \

# Create symbolic links for convenience
ln -s /usr/bin/clang-tidy-$LLVM_VERSION /usr/bin/clang-tidy && \
ln -s /usr/bin/run-clang-tidy-$LLVM_VERSION /usr/bin/run-clang-tidy
