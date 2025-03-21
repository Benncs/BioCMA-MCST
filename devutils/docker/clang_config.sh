#!/bin/bash
# -----------------------------------------------------------------------------
# Author: CASALE Benjamin
# Date: 16/04/2024
# Description: Script to get the right version of clang-tidy and clang-format 
#              for CI and container execution. Not recommended to run this 
#              on a local machine.
# -----------------------------------------------------------------------------


if [ "$#" -ne 1 ]; then
    echo "Usage: $0  <VERSION>"
    exit 1
fi

# Define LLVM version
LLVM_VERSION="$@"

# Download LLVM installation script
wget -q https://apt.llvm.org/llvm.sh 

# Check if download was successful
if [ $? -ne 0 ]; then
    echo "Error: Failed to download LLVM installation script."
    exit 1
fi
mv llvm.sh  /tmp/llvm.sh 
# Make the script executable
chmod +x /tmp/llvm.sh 

# Run the LLVM installation script with the specified version
/tmp/llvm.sh  $LLVM_VERSION all 

# Check if installation was successful
if [ $? -ne 0 ]; then
    echo "Error: LLVM installation failed."
    exit 1
fi

# Install clang-format and clang-tidy
# apt-get install -y clang-format clang-tidy 
apt-get install -y clang-format-$LLVM_VERSION clang-tidy-$LLVM_VERSION 

# Check if installation was successful
if [ $? -ne 0 ]; then
    echo "Error: Failed to install clang-format and clang-tidy."
    exit 1
fi

# # Create symbolic links for convenience
ln -sf /usr/bin/clang-tidy-$LLVM_VERSION /usr/bin/clang-tidy 
ln -sf /usr/bin/clang-format-$LLVM_VERSION /usr/bin/clang-format 
ln -sf /usr/bin/run-clang-tidy-$LLVM_VERSION /usr/bin/run-clang-tidy

# Check if symbolic links were created successfully
if [ $? -ne 0 ]; then
    echo "Error: Failed to create symbolic links."
    exit 1
fi
echo "LLVM, clang-format, and clang-tidy setup completed successfully."
