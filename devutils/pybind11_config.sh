#!/bin/bash
# -----------------------------------------------------------------------------
# Author: CASALE Benjamin
# Date: 16/04/2024
# Description: Script to install pybind11 and create pybind11.pc file.
# -----------------------------------------------------------------------------

# Install pybind11 and its development package
sudo apt-get install -y pybind11-dev python3-pybind11

# Create directory if it doesn't exist
sudo mkdir -p /usr/local/share/pkgconfig

# Create or overwrite pybind11.pc file
sudo tee /usr/local/share/pkgconfig/pybind11.pc > /dev/null <<EOF
prefix=\${pcfiledir}/../../
includedir=\${prefix}/include

Name: pybind11
Description: Seamless operability between C++11 and Python
Version: 2.12.0
Cflags: -I\${includedir}
EOF

