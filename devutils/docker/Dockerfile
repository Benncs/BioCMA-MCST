FROM nvcr.io/nvidia/nvhpc:24.5-devel-cuda_multi-ubuntu22.04
ARG LLVM_VERSION=17


RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y \
    python3-dev \
    python3 \
    python3-pip \
    wget \
    software-properties-common \
    gnupg \ 
    libomp-dev libopenmpi-dev libeigen3-dev  pkg-config ninja-build

RUN apt-get install -y libblas-dev libhdf5-dev 
RUN apt-get install -y libblas-dev libhdf5-dev libvtk9-dev

# Install Meson
RUN pip install meson

# Create and set the working directory
WORKDIR /app
COPY clang_config.sh /app
# Install LLVM clang-tidy/clang-format
RUN sh ./clang_config.sh $LLVM_VERSION



RUN apt-get install -y clang-$LLVM_VERSION

