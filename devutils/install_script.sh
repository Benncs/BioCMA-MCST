#!/bin/bash

# Directory where the application will be installed
INSTALL_DIR="/opt/biomc"
# Check if RESULT_DIR is provided as an argument
if [ -z "$1" ]; then
    echo "Error: RESULT_DIR argument is missing."
    exit 1
fi
RESULT_DIR="$1"

# Create the installation directory if it doesn't exist
mkdir -p "$INSTALL_DIR"
mkdir -p "$RESULT_DIR"

if [ -z "$MESON_SOURCE_ROOT" ]; then
    MESON_SOURCE_ROOT='./'
fi

# Copy additional files
cp $MESON_SOURCE_ROOT/tools/runner.py "$INSTALL_DIR/runner.py"
cp $MESON_SOURCE_ROOT/tools/cli_formater.py "$INSTALL_DIR/cli_formater.py"
cp $MESON_SOURCE_ROOT/devutils/exec.py "$INSTALL_DIR/exec.py"
cp -r $MESON_SOURCE_ROOT/devutils/datamodel "$INSTALL_DIR/datamodel"

# Create a symbolic link for the runner script
ln -sf "$INSTALL_DIR/runner.py" /usr/bin/biomc_runner

# Create the results directory
mkdir -p $RESULT_DIR
