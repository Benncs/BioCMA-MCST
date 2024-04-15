#!/bin/bash

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
