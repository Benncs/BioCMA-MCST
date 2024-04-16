FROM ubuntu:latest
RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y \
    python3-dev \
    python3 \
    python3-pip \
    wget \
    software-properties-common \
    gnupg \ 
    libomp-dev libopenmpi-dev libeigen3-dev  pkg-config ninja-build

# Install Meson
RUN pip install meson

# Create and set the working directory
WORKDIR /app

# Copy all files and directories from the parent directory into /app
COPY . /app

# Install LLVM clang-tidy/clang-format
RUN sh ./devutils/clang_config.sh

# Configure project with Meson
RUN mkdir -p __builddir && \
    meson setup __builddir/ --buildtype=release

# Format the code using clang-format
RUN ninja -C __builddir/ clang-format 

# Run clang-tidy on the codebase
RUN ninja -C __builddir/ clang-tidy